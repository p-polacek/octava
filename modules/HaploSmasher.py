#!/usr/bin/env python

import re
import json
import os.path
import gffutils
import numpy as np

from collections import OrderedDict
from sortedcontainers import SortedDict
from gffutils.exceptions import FeatureNotFoundError

from .parseVcf import VcfReader
from .getQueryOptions import (getDbDir,
                              getFile,
                              getVariantConsequences,
                              convertSampleNames)

"""
Haplosmasher v2.0 main pipeline
v1.0 author: Martijn van Kaauwen
v2.0 author: Peter Poláček

Usage: To be called from within HaploSmasherApp.py
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


def HaploSmasherPipeline(genome, form):
    """Main HaploSmasher pipeline"""

###############################################################################
### Initialise output structures
###############################################################################

    genotypes      = list()  # [<sample1>[pos1, pos2, ...], <sample2>[], ...]

    groups         = SortedDict()  # {group: [sample_1, ..., sample_n]}

    genes          = OrderedDict() # [{gene_id: {"ID":          gene_id,
                                   #             "Chromosome":  chromosome,
                                   #             "Start":       start,
                                   #             "End":         end,
                                   #             "Length":      length,
                                   #             "Strand":      strand,
                                   #             "Promoters":   [(start,
                                   #                              end,
                                   #                              id)],
                                   #             "Variants":    var_count}}]

    annotations    = list()  # [<pos1>[chromosome,
                                   #         position,
                                   #         short_ref,
                                   #         var_type,
                                   #         [affected_gene1, ...],
                                   #         [short_ann1, ...]], <pos2>[], ...]

                                   # Converted to lists of dictionaries later

    errors         = list()

###############################################################################
### Initialise input parameters
###############################################################################

    query = form.data
    # Format
    if query["fmt"] == "gene":
        if query["gene"]:
            query["gene"] = [gene.strip() for gene in
                             query["gene"].split(", ")]
        else:
            query["gene"] = []

        for group in query["gene_groups"]:
            query["gene"] += json.loads(group)

    elif query["fmt"] == "region":
        query["region"] = [[query["chrom"], query["start"], query["end"]]]
        _ = query.pop("gene")

    elif query["fmt"] == "coord":
        regions = list()
        for region in query["coord"].split(", "):
            if region.strip():
                if region.count(":") != 1 or region.count("-") > 1:
                    errors.append("Invalid coordinates: " + region)
                    continue
                
                chrom, span = region.split(":")
                start, end = span.split("-") if "-" in span else [span, span]
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
            return (genotypes,
                    groups,
                    genes,
                    annotations,
                    errors)

        query["region"] = regions
        _ = query.pop("gene")

    for group in query["sample_groups"]:
        query["samples"] += group.split(",")

###############################################################################
### Initialise files and conversion dictionaries
###############################################################################

    try:
        # Get target annotation and variant files
        gff_file = getFile(genome, "gff")
        if query["impact"] != ["geneonly"]:
            vcf_file = getFile(genome, "vcf")

        # Initialize gffutils database
        db_dir = getDbDir()
        db = gffutils.FeatureDB(f"{db_dir}/{genome}.db")

    except OSError as error:
        errors.append(error.args[0])

        # Without any of these, the program aborts and returns error message
        # All structures returned to avoid errors with output unpacking
        return (genotypes,
                groups,
                genes,
                annotations,
                errors)

    # Initialise dictionaries for print-formatting consequences and samples
    csq_dict = dict(getVariantConsequences(impact_tag=False))
    try:
        sample_dict = convertSampleNames(genome)
    except:
        sample_dict = None

    # Add "indel" to type_choice for easier filtering later
    if query["types"]:
        if "ins" in query["types"] or "del" in query["types"]:
            query["types"].append("indel")

###############################################################################
### Compose a list of genes within the query
###############################################################################

    target_region = query["target"]

    # If input is a gene identifier
    if "gene" in query.keys():
        records = list()
        for gene in query["gene"]:

            try:
                gene_record = db[gene]
                if gene_record not in records:
                    records.append(gene_record)
                else:
                    errors.append("Duplicate gene entered: " + gene)

            except FeatureNotFoundError as error:
                errors.append("Gene not found: " + error.feature_id)

        if not records:
            return (genotypes,
                    groups,
                    genes,
                    annotations,
                    errors)

        # Compose a list of records
        for record in records:
            if target_region == "prom":
                # Get promoter upstream of transcription start
                # (for genes with multiple transcripts, get each promoter)
                transcripts = db.children(record.id,
                                          featuretype="mRNA")
                prom_size = query["prom_size"]
                promoters = [(
                              (t.start - prom_size,
                               t.start - 1,
                               t["ID"][0])
                              if t.strand == "+" else
                              (t.end + 1,
                               t.end + prom_size,
                               t["ID"][0])
                              )
                             for t in transcripts]
            else:
                promoters = None

            genes[record.id] = {"ID":              record.id,
                                "Chromosome":      record.chrom,
                                "Start":           record.start,
                                "End":             record.end,
                                "Length":          abs(record.start -
                                                       record.end) + 1,
                                "Strand":          record.strand,
                                #"Promoters":       promoters,
                                "Variants":        0,
                                "Region":          None}


            for key, value in record.attributes.items():
                if key != "ID":
                    genes[record.id][key] = value[0]

    # If input is a chromosome region
    elif "region" in query.keys():
        for i, region in enumerate(query["region"]):
            records = db.region(region=(region),
                                featuretype="gene",
                                completely_within=False)

            # Compose a list of records
            for record in records:
                if target_region == "prom":
                    # Get promoter upstream of transcription start
                    # (for genes with multiple transcripts, get each promoter)
                    transcripts = db.children(record.id,
                                              featuretype="mRNA")
                    prom_size = query["prom_size"]
                    promoters = [(
                                  (t.start - prom_size,
                                   t.start - 1,
                                   t["ID"][0])
                                  if t.strand == "+" else
                                  (t.end + 1,
                                   t.end + prom_size,
                                   t["ID"][0])
                                  )
                                 for t in transcripts]
                else:
                    promoters = None

                genes[record.id] = {"ID":              record.id,
                                    "Chromosome":      record.chrom,
                                    "Start":           record.start,
                                    "End":             record.end,
                                    "Length":          abs(record.start -
                                                           record.end) + 1,
                                    "Strand":          record.strand,
                                    #"Promoters":       promoters,
                                    "Variants":        0,
                                    "Region":          i}

                for key, value in record.attributes.items():
                    if key != "ID":
                        genes[record.id][key] = value[0]

###############################################################################
### Set the regions to search through
###############################################################################

    if target_region == "full":

        if "gene" in query.keys():
            # Get each gene from start to end,
            # with no entry for affected promoter
            regions = [[gene["Chromosome"],
                        gene["Start"],
                        gene["End"],
                        gene["ID"],
                        None]
                       for gene in genes.values()]

        else:
            # Get full query region,
            # with no entry for affected gene or transcript
            regions = [[region[0],
                        region[1],
                        region[2],
                        None,
                        None] for region in query["region"]]

    else:
        if target_region == "cds":

            # Get each gene from start to end,
            # with no entry for affected transcript
            regions = [[gene["Chromosome"],
                        gene["Start"],
                        gene["End"],
                        gene["ID"],
                        None,
                        gene["Region"]]
                       for gene in genes.values()]

        elif target_region == "prom":

            # Get each promoter from start to end,
            # with entry for affected gene and affected transcript
            regions = list()
            for gene in genes.values():
                for promoter in gene["Promoters"]:
                    regions.append([gene["Chromosome"],
                                    promoter[0],
                                    promoter[1],
                                    gene["ID"],
                                    promoter[2]])

        # Clip regions partly outside query
        if "region" in query.keys():
            regions = [[region[0],
                        max(region[1], query["region"][region[5]][1]),
                        min(region[2], query["region"][region[5]][2]),
                        region[3],
                        region[4],
                        region[5]]
                       for region in regions]

        # Information on which gene came from which region kep for clipping
        # Removed after clipping finished
        regions = [region[:-1] for region in regions]

###############################################################################
### Iterate over variant records in regions
###############################################################################

    # Check if vcf file has been loaded
    if "vcf_file" not in locals():

        return (genotypes,
                groups,
                genes,
                annotations,
                errors)

    vcf_reader = VcfReader(vcf_file)

    # Filter by sample
    if query["samples"]:
        samples = [sample for sample
                   in vcf_reader.header.samples
                   if sample in query["samples"]]
    else:
        samples = vcf_reader.header.samples

    # Initiate sample entries beforehand, to avoid running setdefault()
    # for every call of every sample
    for sample in samples:
        genotypes.append([])

    ann_positions = set()
    # Iterate over regions
    for region in regions:
        chrom = region[0]
        # VcfReader does not consider first position of a set region
        vcf_reader.set_region(chrom,
                              region[1] - 1,
                              region[2])

        # Last variant position kept to catch duplicate variant entries
        last_pos = 0
        last_pos_offset = 0

        for record in vcf_reader:
            # Filter by type
            if query["types"] and (record.VAR_TYPE
                                   not in query["types"]):
                continue

###############################################################################
### Iterate over annotations, extract those matching query
###############################################################################

            ann_list = list()
            pos = record.POS
            ref = record.REF
            typ = record.VAR_TYPE

            # A position can have multiple entries in VCF file
            # Add an index to the position to distinguish them
            # e.g. 104829, 104829.1, 104829.2 etc.
            if pos == last_pos:
                last_pos_offset += 1
                pos = ".".join(map(str(last_pos, last_pos_offset)))
            else:
                last_pos_offset = 0
            last_pos = record.POS

            # In case regions to look through have overlap,
            # add the genes/promoters of that region to Affected genes
            # (the rest of the data is the same)
            if (chrom, pos) in ann_positions:
                for ann_index in range(len(annotations)):
                    ann = annotations[ann_index]
                    if chrom != ann[0] or pos != ann[1]:
                        continue
                    if region[3] not in ann[2]:
                        annotations[ann_index][2].append(region[3])
                    if target_region == "prom":
                        annotations[ann_index][2].append(region[4])                 ## CORRECT THIS LATER
                continue

            # Annotation may be missing (record.CSQ returns an error)
            # Handled below by keeping empty annotation entry
            try:
                for ann in record.CSQ:
                    if ann["Annotation_Impact"] not in query["impact"]:
                        continue

                    if query["csqs"]:
                        if not [c for c in ann["Annotation"].split("&")
                                if c in query["csqs"]]:
                            continue

                    ann_list.append(ann)

                if not ann_list:
                    continue

                # Iterate again to have at least one annotation
                # for each alternative variant
                alts = {ann["Allele"] for ann in ann_list}
                if len(alts) < len(record.ALT):
                    for ann in record.CSQ:
                        if ann["Allele"] not in alts:
                            ann_list.append(ann)
                            alts.add(ann["Allele"])

                # Sort based on the order in record.INFO field
                ann_list.sort(key=record.CSQ.index)

###############################################################################
### Shorten long sequences in annotations (for html visualisation)
###############################################################################

                print_anns = list()
                for ann in ann_list:
                    print_ann = list()
                    # Add variant index
                    try:
                        var_num = record.ALLELES.index(ann["Allele"])
                        print_ann.append(typ + str(var_num))
                    except:
                        # Sometimes given allele is not found in record.ALT
                        print_ann.append("Unknown")

                    # Shorten ALT allele
                    alt = ann["Allele"]
                    if len(alt) > 5:
                        alt = compress_sv(alt)
                    print_ann.append(alt)

                    # Split consequence names
                    print_csq = [csq_dict[c]
                                 if c in csq_dict.keys() else c
                                 for c in ann["Annotation"].split("&")]
                    print_ann.append(", ".join(print_csq))

                    # Shorten AA change
                    if ann["HGVS.p"]:
                        aa = ann["HGVS.p"]
                        if len(ref) > 5 or len(alt) > 5:
                            aa = re.sub("(([A-Z][a-z]{2})+)",
                                        lambda x: (compress_sv(x.group())
                                                   if len(x.group()) > 12
                                                   else x.group()),
                                        aa)
                        aa = re.sub("^p\.", "", aa)
                        print_ann.append(aa)

                    # Format warnings
                    if ann["ERRORS / WARNINGS / INFO"]:
                        warnings = ann["ERRORS / WARNINGS / INFO"]
                        print_ann.append(", ".join(warnings.split("&")))

                    print_anns.append(" | ".join(print_ann))

            except:
                print_anns = ""

###############################################################################
### Get and process severity score (if available)
###############################################################################

            prov_anns = list()
            for prov_ann in record.parsed_info_fields().get("PROV", ""):
                fields = prov_ann.split("|")
                try:
                    var_num = record.ALLELES.index(fields[0])
                    var = typ + str(var_num)
                except:
                    var = "Unknown"
                score = str(round(float(fields[7]), 2))
                prov_anns.append(" | ".join((var, fields[0], score, fields[11])))

###############################################################################
### Save variant entry
###############################################################################

            if target_region == "prom":
                affected_genes = [region[4]]
            else:
                if region[3]:
                    affected_genes = [region[3]]
                else:
                    affected_genes = sorted({a["Gene_ID"] for a in ann_list})

            # Save annotation entry
            if len(ref) > 5:
                ref = compress_sv(ref)

            annotations.append([chrom, record.POS, affected_genes, ref, print_anns, prov_anns])
            ann_positions.add((chrom, pos))

###############################################################################
### Get genotype data
###############################################################################

            for sample_index, sample in enumerate(samples):
                data = record.get_sample_call(sample)
                gt = data["GT"]

                for delim in ("/", "|", "\\"):
                    if delim in gt:
                        a, b = gt.split(delim)
                        a = "·" if a == "0" else "—" if a == "." else typ + a
                        b = "·" if b == "0" else "—" if b == "." else typ + b
                        genotype = delim.join((a, b)) if a != b else a
                        break
                else:
                    genotype = "·" if gt == "0" else "—" if gt == "." else gt

###############################################################################
### Get allelic depths and qualities for a call
###############################################################################

                try:
                    # freeBayes (and others) report genotype likelihoods
                    if "GL" in data.keys():
                        # Genotype likelihoods need to be Phred-scaled (GL * -10)
                        # Calculated from most likely genotype except the called one
                        # (which is typically 0)
                        if isinstance(data["GL"], str):
                            qual = sorted(list(map(int, data["GL"].split(","))))[1] * -10
                        else:
                            qual = sorted(data["GL"])[1] * -10

                    # bcftools  (and others) report genotype quality
    #                elif "GQ" in data.keys():
    #                    qual = data["GQ"]
                    elif "PL" in data.keys():
                        if isinstance(data["PL"], str):
                            qual = sorted(list(map(int, data["PL"].split(","))))[1]
                        else:
                            qual = sorted(data["PL"])[1]

                    else:
                        qual = ""
                except:
                    qual = ""

                genotypes[sample_index].append([genotype, qual])

###############################################################################
### Filter out positions without variation
###############################################################################

    pos_remove = list()
    for pos in range(len(genotypes[0])):
        # Check if position has genotypes apart from those to be filtered out
        calls = {sample_gts[pos][0] for sample_gts in genotypes}

        # query["skip_variants"] can be ".", "—", or both
        # Check if anything else is present
        if calls.issubset(query["skip_variants"]):
            pos_remove.append(pos)

        # If looking only for variants varying between samples
        # and at least one sample differs from others at the position
        elif query["skip_no_var"] and len(calls) == 1:
            pos_remove.append(pos)

    annotations = [a for i, a in enumerate(annotations) if i not in pos_remove]
    genotypes = [[gt for i, gt in enumerate(sample_gts) if i not in pos_remove]
                 for sample_gts in genotypes]

    # If there are no genotypes left after filtering
    if not genotypes[0]:
        # genotypes are checked if true when rendering template - set to empty
        genotypes = list()

        return (genotypes,
                groups,
                genes,
                annotations,
                errors)

###############################################################################
### Group samples based on haplotype
###############################################################################

    varnum = len(genotypes[0])
    for sample, sample_gts in zip(samples, genotypes):
        haplotype = tuple(sample_gts[pos][0] for pos in range(varnum))

        sample_name = sample_dict[sample] if sample_dict else sample
        groups.setdefault(haplotype, list()).append(sample_name)

    #groups = list(groups.values())

###############################################################################
### Filter out samples without variation
###############################################################################

    sample_remove = list()
    # If looking for differences between samples,
    # all samples will always be visible
    if not query["skip_no_var"]:
        for sample_index, (sample, sample_gts) in enumerate(zip(samples, genotypes)):
            calls = {sample_gts[pos][0] for pos in range(varnum)}

            if calls.issubset(query["skip_samples"]):
                sample_remove.append(sample_index)

        genotypes = [gt for i, gt in enumerate(genotypes)
                     if i not in sample_remove]
        samples = [sample for i, sample in enumerate(samples)
                   if i not in sample_remove]

###############################################################################
### Keep track of variant count
###############################################################################

    for ann in annotations:
        for gene in ann[2]:
            try:
                genes[gene]["Variants"] += 1
            except:
                pass

###############################################################################
### Rename sample names and append to genotypes
###############################################################################

    if sample_dict:
        samples = [sample_dict[s] if s in sample_dict.keys() else s
                   for s in samples]

    # Sort positions
    order = sorted(range(len(annotations)), key=lambda x: (annotations[x][0], 
                                                           annotations[x][1]))
    annotations = [annotations[i] for i in order]
    genotypes = [[sample_data[i] for i in order] for sample_data in genotypes]
    groups = [([k[i] for i in order], v) for k, v in groups.items()]

    genotypes = [[sample] + gt for sample, gt in zip(samples, genotypes)]

###############################################################################
### Return
###############################################################################

    return (genotypes,
            groups,
            genes,
            annotations,
            errors)

if __name__ == "__main__":
    pass
