#!/usr/bin/env python

import os
import re
import json
import gffutils
import bed_reader
import numpy as np
import pandas as pd

from .csq import var_csq
from .vcf_reader import VcfReader
from .load_input import load_config

from flask_app import (
  get_lock,
  master_lock)

from app.dirs import data_dir, tool_dir

"""
Allele mining pipeline
v1.0 author: Martijn van Kaauwen
v2.0 author: Peter Poláček

Usage: To be called from within app.py
"""

def compress_sv(seq):
    """Shortens nucleotide sequence when longer than 5

    seq: str, nucleotide sequence to shorten

    returns: str, shortened sequence
    """
    # If amino acid sequence (ends with Xxx)
    if re.findall("[A-Z][a-z]{2}$", seq):
        # Convert to format XxxXxx..XxxXxx (n AA)
        return f"{seq[:6]}..{seq[-6:]} ({len(seq)//3} AA)"

    # If nucleotide sequence (ends with X)
    elif re.findall("[A-Z]$", seq):
        # Convert to format XX..XX (n bp)
        return f"{seq[:2]}..{seq[-2:]} ({len(seq)} bp)"
    else:
        return seq

def genes_list(genes, gene_attrs):
    # Convert genes and annotations from dict to list, add header
    if genes:
        genes = list(genes.values())
        gene_attrs = [
        "ID",
        "Chromosome",
        "Start",
        "End",
        "Length",
        "Strand",
        "Variants"] + gene_attrs
        genes = [gene_attrs] + genes
    return genes

def am_pipeline(genome, form):
    """Pipeline generating allele mining output"""

###############################################################################
### Initialise output structures
###############################################################################

    # [[sample_1, [pos1, qual1], [pos2, qual2], ...], [sample2, ...]]
    gts = list()

    # {(gt1, gt2, ...): [sample_1, sample_2, ...], ...}
    # Later converted to list of (key, value) tuples
    groups = dict()

    # {gene_1: [gene_attributes], ...}
    # Later converted to list of [gene_attributes]
    genes = dict()

    # [[chr1, pos1, [affected_genes], ref, [csq_list], [sev_list]], ...]
    anns = list()

    errors = list()

###############################################################################
### Initialise files and conversion dictionaries
###############################################################################

    gconfig = load_config(genome)
    query = form.data

    try:
        if query["impact"] != ["geneonly"]:
            if query["dataset"] == "":
                vcf_path = os.path.join(
                  data_dir,
                  genome,
                  gconfig["vcf"])
            else:
                markers_root = os.path.join(
                  data_dir,
                  genome,
                  "assoc",
                  query["dataset"])
                vcf_path = os.path.join(
                  markers_root,
                  "gt",
                  "gt.vcf.gz")
                params_path = os.path.join(
                  markers_root,
                  "stats",
                  "params")
                with open(params_path, "r") as params_in:
                    is_hap = json.load(params_in)[1][11] == "Yes"
            vcf_reader = VcfReader(vcf_path)

            sev_path = os.path.join(
              data_dir,
              genome,
              f"sev.{gconfig['vcf']}")
            if os.path.exists(sev_path):
                sev_reader = VcfReader(sev_path)
            else:
                sev_reader = None
    except Exception as e:
        errors.append("Problem loading VCF file.")
        
    try:
        db_path = os.path.join(
          data_dir,
          genome,
          f"{genome}.db")
        db = gffutils.FeatureDB(db_path)
    except:
        errors.append("Problem loading GFFutils DB file.")

    if errors:
        # Without any of these, the program aborts and returns error message
        # All structures returned to avoid errors with output unpacking
        return (
          gts,
          groups,
          genes,
          anns,
          errors)

    # Initialise dictionaries for print-formatting consequences and samples
    csq_dict = dict(var_csq(impact_tag=False))

    samp_path = os.path.join(
      data_dir,
      genome,
      gconfig["samp"])
    with master_lock:
        lock = get_lock(samp_path)
    with lock:
        samp_dict = (pd
          .read_csv(
            samp_path,
            sep="\t",
            header=0,
            index_col=0,
            usecols=[0,1])
          .set_axis(
            ["Name"],
            axis=1)["Name"]
          .to_dict())

    for group in query["samples_group"]:
        query["samples"] += json.loads(group)

    # Get attribute names of additional information for genes
    if gconfig["info"]:
        info_path = os.path.join(
          data_dir,
          genome,
          gconfig["info"])
        with open(info_path, "r") as info_in:
            gene_attrs = next(info_in).strip().split("\t")[1:]
    else:
        gene_attrs = []

    # Add "indel" to variant types choice for filtering later
    if query["types"]:
        if "ins" in query["types"] or "del" in query["types"]:
            query["types"].append("indel")

    target = query["target"]
    prom_size = query["prom_size"]

    regions = list()

###############################################################################
### Get list of genes and regions to search through - gene input
###############################################################################

    if query["fmt"] == "gene":
        if query["gene"]:
            query["gene"] = [
              gene.strip() for gene in
              query["gene"].split(", ")]
        else:
            query["gene"] = []

        for group in query["gene_groups"]:
            query["gene"] += json.loads(group)

        recs = list()
        for gene in query["gene"]:
            try:
                rec = db[gene]
                if rec not in recs:
                    recs.append(rec)
                else:
                    errors.append("Duplicate gene entered: " + gene)

            except gffutils.exceptions.FeatureNotFoundError as error:
                errors.append("Gene not found: " + error.feature_id)
        if not recs:
            return (
              gts,
              groups,
              genes,
              anns,
              errors)

        # Compose a list of records
        for rec in recs:
            if target == "prom":
                # Get promoter upstream of transcription start
                # (for genes with multiple transcripts, get each promoter)
                transcripts = db.children(
                  rec.id,
                  featuretype="mRNA")

                for t in transcripts:
                    if t.strand == "+":
                        start = t.start - prom_size
                        end = t.start - 1
                    else:
                        start = t.end + 1
                        end = t.end + prom_size
                    regions.append([rec.chrom, start, end])
            else:
                regions.append([rec.chrom, rec.start, rec.end])

            genes[rec.id] = [
              rec.id,
              rec.chrom,
              rec.start,
              rec.end,
              abs(rec.start - rec.end) + 1,
              rec.strand,
              0 if query["impact"] != ["geneonly"] else "-",
              *[
                rec.attributes.get(attr, [""])[0]
                for attr in gene_attrs]
            ]

###############################################################################
### Get list of genes and regions to search through - region/coord input
###############################################################################

    else:
        if query["fmt"] == "coord":
            for region in query["coord"].split(", "):
                if region.strip():
                    if region.count(":") != 1 or region.count("-") > 1:
                        errors.append("Invalid coordinates: " + region)
                        continue
                    
                    chrom, span = region.split(":")
                    start, end = span.split("-") if "-" in span else (span, span)
                    try:
                        start = int(start.replace(",", "").replace(" ", ""))
                        end   = int(end.replace(",", "").replace(" ", ""))
                        assert start <= end
                    except:
                        errors.append("Invalid coordinates: " + region)
                        continue

                    chrom = chrom.strip()
                    if chrom not in form.chrom.choices:
                        errors.append("Invalid coordinates: " + region)
                        continue

                    regions.append([chrom, start, end])

            if not regions:
                return (gts,
                        groups,
                        genes,
                        anns,
                        errors)

            query["region"] = regions[:]
            regions = list()

        else:
            query["region"] = [(query["chrom"], query["start"], query["end"])]

        for region in query["region"]:
            if target == "full":
                regions.append(region)

            recs = db.region(
              region=(region),
              featuretype="gene",
              completely_within=False)

            # Compose a list of records
            for rec in recs:
                if target == "prom":
                    # Get promoter upstream of transcription start
                    # (for genes with multiple transcripts, get each promoter)
                    transcripts = db.children(
                      rec.id,
                      featuretype="mRNA")

                    for t in transcripts:
                        if t.strand == "+":
                            start = t.start - prom_sise
                            end = t.start - 1
                        else:
                            start = t.end + 1
                            end = t.end + prom_size
                        # Clip region if partly outside query
                        regions.append(
                          [rec.chrom,
                           max(region[1], start),
                           min(region[2], end)])
                else:
                    # Clip region if partly outside query
                    regions.append(
                      [rec.chrom,
                       max(region[1], rec.start),
                       min(region[2], rec.end)])

                genes[rec.id] = [
                  rec.id,
                  rec.chrom,
                  rec.start,
                  rec.end,
                  abs(rec.start - rec.end) + 1,
                  rec.strand,
                  0 if query["impact"] != ["geneonly"] else "-",
                  *[
                    rec.attributes.get(attr, [""])[0]
                    for attr in gene_attrs]
                ]

###############################################################################
### Iterate over variant records in regions
###############################################################################

    if query["impact"] == ["geneonly"]:
        genes = genes_list(genes, gene_attrs)
        return (
          gts,
          groups,
          genes,
          anns,
          errors)

    ann_positions = set()

    # Filter by sample
    if query["samples"]:
        samples = [
          sample for sample
          in vcf_reader.header.samples
          if sample in query["samples"]]
    else:
        samples = vcf_reader.header.samples

    # Initiate sample entries beforehand
    print(query["samples"])
    for sample in samples:
        gts.append([])

    for region in regions:
        chrom = region[0]
        # VcfReader is 0-based
        vcf_reader.set_region(
          chrom,
          region[1] - 1,
          region[2])

        # Last variant position kept to catch duplicate variant entries
        prev_pos = 0
        dup_cnt = 0

###############################################################################
### Filter for type and severity
###############################################################################

        for var in vcf_reader:
            # Filter by type
            if query["dataset"] == "":
                if query["types"] and (var.VAR_TYPE not in query["types"]):
                    continue

                if "SEV" in query["impact"]:
                    if all([
                      sev.rsplit("|", 1)[1] == "TOLERATED"
                      for sev in var.parsed_info_fields().get("SEV", "")]):
                        continue

###############################################################################
### Iterate over annotations, extract those matching query
###############################################################################

            ann_list = list()
            pos = var.POS
            ref = var.REF
            typ = var.VAR_TYPE

            # Multiple vars per position are given an index to distinguish
            if pos == prev_pos:
                dup_cnt += 1
                pos = ".".join(map(str, (pos, dup_cnt)))
            else:
                dup_cnt = 0
                prev_pos = pos

            # For overlapping regions, update only affected genes/promoters
            # rest of the data is the same
            if (chrom, pos) in ann_positions:
                continue

            # Annotation may be missing (var.CSQ returns an error)
            # Annotation entry kept empty
            print_anns = list()
            try:
                csq = var.CSQ
            except:
                pass
            else:
                for ann in csq:
                    if ann["Annotation_Impact"] not in query["impact"]:
                        continue

                    if query["csqs"]:
                        if not [
                          c for c in ann["Annotation"].split("&")
                          if c in query["csqs"]]:
                            continue

                    ann_list.append(ann)

                if ann_list:
                    # Iterate again to have at least one annotation
                    # for each alternative variant
                    alts = {ann["Allele"] for ann in ann_list}
                    if len(alts) < len(var.ALT):
                        for ann in csq:
                            if ann["Allele"] not in alts:
                                ann_list.append(ann)
                                alts.add(ann["Allele"])

                    # Sort based on the order in var.INFO field
                    ann_list.sort(key=var.CSQ.index)
                else:
                    continue

###############################################################################
### Shorten long sequences in annotations (for html visualisation)
###############################################################################

                for ann in ann_list:
                    print_ann = list()
                    # Add variant index
                    try:
                        var_num = var.ALLELES.index(ann["Allele"])
                        print_ann.append(typ + str(var_num))
                    except:
                        # Sometimes given allele is not found in var.ALT
                        print_ann.append("Unknown")

                    # Shorten ALT allele
                    alt = ann["Allele"]
                    if len(alt) > 5:
                        alt = compress_sv(alt)
                    print_ann.append(alt)
                    print_ann.append(ann["Annotation_Impact"])

                    print_anns.append(" | ".join(print_ann))
                    print_ann = []

                    # Split consequence names
                    print_csq = [
                      csq_dict[c]
                      if c in csq_dict.keys() else c
                      for c in ann["Annotation"].split("&")]
                    print_ann.append(", ".join(print_csq))

                    # Shorten AA change
                    if ann["HGVS.p"]:
                        aa = ann["HGVS.p"]
                        if len(ref) > 5 or len(alt) > 5:
                            aa = re.sub(
                              "(([A-Z][a-z]{2})+)",
                              lambda x: (
                                compress_sv(x.group())
                                if len(x.group()) > 12
                                else x.group()),
                              aa)
                        aa = re.sub(r"^p\.", "", aa)
                        print_ann.append(aa)

                    # Format warnings
                    if ann["ERRORS / WARNINGS / INFO"]:
                        warnings = ann["ERRORS / WARNINGS / INFO"]
                        print_ann.append(", ".join(warnings.split("&")))

                    print_anns.append(" | ".join(print_ann))

###############################################################################
### Get and process severity score (if available)
###############################################################################

                if not print_anns:
                    continue

            sevs = list()
            if sev_reader and query["dataset"] == "":
                sev_reader.set_region(
                  var.CHROM,
                  var.POS - 1,
                  var.POS)
                for sev in sev_reader:
                    if sev.CHROM != var.CHROM or sev.POS != var.POS:
                        continue
                    for sev_ann in sev.parsed_info_fields().get("SEV", ""):
                        data = sev_ann.split("|")
                        try:
                            var_num = var.ALLELES.index(data[0])
                            str_var = typ + str(var_num)
                        except:
                            str_var = "Unknown"
                        score = str(round(float(data[7]), 2))
                        sevs.append(
                          " | ".join(
                            (str_var,
                             data[0],
                             score)))
                        sevs.append(
                             data[11])

###############################################################################
### Save variant entry
###############################################################################

            # Save annotation entry
            if query["dataset"] == "":
                if len(ref) > 5:
                    ref = compress_sv(ref)
            else:
                if is_hap:
                    ref = f"({len(ref)} bp)"

            anns.append(
              [chrom,
               pos,
               sorted({ann["Gene_ID"] for ann in ann_list}),
               ref,
               print_anns,
               sevs])
            ann_positions.add((chrom, pos))

###############################################################################
### Get genotype data
###############################################################################

            min_c = 0
            mis_c = 0
            het_c = 0
            for samp_idx, sample in enumerate(samples):
                data = var.get_sample_call(sample)
                gt = data["GT"]

                for delim in ("/", "|", "\\"):
                    if delim in gt:
                        a, b = gt.split(delim)
                        if a == "." or b == ".":
                            mis_c += 1
                        else:
                            min_c += sum(1 for v in (a, b) if v != "0")
                            if a != b:
                                het_c += 1
                        a = "·" if a == "0" else "—" if a == "." else typ + a
                        b = "·" if b == "0" else "—" if b == "." else typ + b
                        genotype = delim.join((a, b)) if a != b else a
                        break
                else:
                    genotype = "·" if gt == "0" else "—" if gt == "." else gt

###############################################################################
### Get allelic depths and qualities for a call
###############################################################################

                if query["dataset"] == "":
                    quals = list()
                    qual = ""
                    if "GQ" in data.keys():
                        try:
                            qual = int(data["GQ"])
                        except ValueError:
                            pass

                    if qual == "" and "PL" in data.keys():
                        for q in data["PL"].split(","):
                            try:
                                quals.append(int(q))
                            except ValueError:
                                pass
                        qual = sorted(quals)[1] if quals else ""

                    if qual == "" and "GL" in data.keys():
                        # Genotype likelihoods converted to Phred-scale
                        for q in data["GL"].split(","):
                            try:
                                quals.append(float(q))
                            except ValueError:
                                pass
                        qual = sorted(quals)[1] * -10 if quals else ""
                else:
                    qual = ""

                gts[samp_idx].append([genotype, qual])

            gen_c = (len(samples) - mis_c) * 2
            mis = mis_c / len(samples)
            if gen_c > 0:
                het = het_c / gen_c
                maf = min(min_c / gen_c, 1 - min_c / gen_c)
            if (gen_c == 0 and query["mis"]) or (
              (query["mis"] and mis > query["mis"]) or 
              (query["het"] and het > query["het"]) or
              (query["maf"] and maf < query["maf"])):
                _ = anns.pop()
                gts = [gt[:-1] for gt in gts]
            else:
                if gen_c == 0:
                    anns[-1] += [np.round(mis, 3), "—", "—"]
                else:
                    anns[-1] += [
                      np.round(mis, 3),
                      np.round(het, 3),
                      np.round(maf, 3)]

###############################################################################
### Filter out positions without variation
###############################################################################

    pos_rmv = list()
    for pos in range(len(gts[0])):
        # Check if pos has gts other than those to be filtered out (".", "—")
        calls = {samp_gts[pos][0] for samp_gts in gts}
        if calls.issubset(query["skip_variants"]):
            pos_rmv.append(pos)

        # If skipping variants with no variation
        elif query["skip_no_var"] and len(calls) == 1:
            pos_rmv.append(pos)

    anns = [
      ann for idx, ann in enumerate(anns)
      if idx not in pos_rmv]
    gts  = [
      [gt for idx, gt in enumerate(sample_gts)
       if idx not in pos_rmv]
      for sample_gts in gts]

    # If there are no valid genotypes
    if not gts[0]:
        gts = list()
        genes = genes_list(genes, gene_attrs)
        return (
          gts,
          groups,
          genes,
          anns,
          errors)

###############################################################################
### Group samples based on haplotype
###############################################################################

    var_num = len(gts[0])
    for sample, sample_gts in zip(samples, gts):
        haplotype = tuple(sample_gts[pos][0] for pos in range(var_num))

        name = samp_dict[sample] if samp_dict else sample
        groups.setdefault(haplotype, list()).append(name)

###############################################################################
### Filter out samples without variation
###############################################################################

    samp_rmv = list()
    # If looking for differences between samples,
    # all samples will be visible
    if not query["skip_no_var"]:
        for samp_idx, (samp, samp_gts) in enumerate(zip(samples, gts)):
            calls = {samp_gts[pos][0] for pos in range(var_num)}

            if calls.issubset(query["skip_samples"]):
                samp_rmv.append(samp_idx)

        gts = [
          gt for idx, gt in enumerate(gts)
          if idx not in samp_rmv]
        samples = [
          samp for idx, samp in enumerate(samples)
          if idx not in samp_rmv]

###############################################################################
### Keep track of variant count
###############################################################################

    for ann in anns:
        for gene in ann[2]:
            if gene in genes.keys():
                genes[gene][6] += 1

###############################################################################
### Rename sample names and append to genotypes
###############################################################################

    if samp_dict:
        samples = [samp_dict.get(samp, samp) for samp in samples]

    # Sort positions
    order = sorted(
              range(len(anns)),
              key=lambda x: (anns[x][0], int(str(anns[x][1]).split(".", 1)[0])))

    anns   = [anns[idx] for idx in order]
    for i in range(len(anns)):
        if "." in str(anns[i][1]):
            idx = str(anns[i][1]).split(".", 1)[1]
        else:
            idx = ""
        anns[i][1] = f"{int(str(anns[i][1]).split('.', 1)[0]):,}"
        if idx:
            anns[i][1] += "." + idx

    gts    = [[samp_data[idx] for idx in order] for samp_data in gts]
    groups = [([k[idx] for idx in order], v) for k, v in groups.items()]

    # Add sample names to gts
    gts = [[sample] + gt for sample, gt in zip(samples, gts)]

###############################################################################
### Return
###############################################################################

    genes = genes_list(genes, gene_attrs)
    return (
      gts,
      groups,
      genes,
      anns,
      errors)

if __name__ == "__main__":
    pass
