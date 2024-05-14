#!/usr/bin/env python

import os
import re
import sys
import json
import pysam
import gffutils
import pandas as pd
from collections import Counter

from app.dirs import data_dir

from modules.assoc import run_qc
from modules.vcf_reader import VcfReader
from modules.load_input import load_config

def configure(genome):
    gconfig = load_config(genome)
    if not all([
      gconfig[typ]
      for typ in ("ann", "vcf", "seq")]):
        print("ANN, VCF and SEQ files are all required.")
        raise OSError

    if not os.path.exists(gconfig["gff_db"]):
        os.mkdirs(gconfig["gff_db"])
    # Create a gffutils database if none exists
    ann_path = os.path.join(
      data_dir,
      genome,
      gconfig["ann"])
    db_path = os.path.join(
      gconfig["gff_db"],
      f"{genome}.db")
    if not os.path.exists(db_path):
        gffutils.create_db(
          ann_path,
          dbfn=db_path,
          merge_strategy="merge")
    db = gffutils.FeatureDB(db_path)

    # Add additional attributes to GFFUtils database
    # Count occurrences per attribute value for functional enrichment
    if gconfig["gene_info"]:
        annotate_db(genome, gconfig)
    if gconfig["enr_cols"]:
        num_genes = db.count_features_of_type("gene")
        calc_enr(genome, gconfig, num_genes)

    if gconfig["phe"] == ".tsv":
        phe_path = os.path.join(
          data_dir,
          genome,
          ".tsv")
        with open(phe_path, "w") as phe_out:
            phe_out.write("\n")

    vcf_path = os.path.join(
      data_dir,
      genome,
      gconfig["vcf"])

    # Create a genome data file if none exists
    gdata_path = os.path.join(
      data_dir,
      genome,
      ".gdata")
    if not os.path.exists(gdata_path):
        gdata = dict()
        # Get chromosome names, sizes and offsets
        idx_path = os.path.join(
          data_dir,
          genome,
          gconfig["seq"] + ".fai")
        idx_df = pd.read_csv(
          idx_path,
          sep="\t",
          header=None)
        chroms = idx_df.iloc[:, 0].tolist()
        sizes = idx_df.iloc[:, 1].tolist()

        offset = 0
        offsets = list()
        for size in sizes:
            offsets.append(offset)
            offset += size

        gdata["chr_names"] = chroms
        gdata["chr_sizes"] = sizes
        gdata["chr_offsets"] = offsets

        # Get first gene ID in the file
        for rec in db.region(
          region=chroms[0],
          featuretype="gene"):
            gdata["id_fmt"] = rec.id
            break

        # Save sample list
        vcf_reader = VcfReader(vcf_path)
        samples = vcf_reader.header.samples

        vcf_header = str(vcf_reader.header)
        vcf_chrs = re.findall(
          "\#\#contig\=\<ID\=(.+?)\,l",
          vcf_header)
        vcf_sizes = re.findall(
          "length\=(\d+)",
          vcf_header)
        vcf_sizes = list(map(int, vcf_sizes))
        if vcf_chrs != chroms or vcf_sizes != sizes:
            print(
              "Chromosomes do not match in size or names between "
              "VCF and SEQ file.")
            raise OSError

        # Rename samples based on names specified
        if gconfig["samp"]:
            samp_path = os.path.join(
              data_dir,
              genome,
              gconfig["samp"])

            with open(samp_path, "r") as samp_in:
                samp_dict = dict(
                  line.strip().split("\t")[:2]
                  for line in samp_in)

            gdata["samples"] = [
              (sample, samp_dict.get(sample, sample))
              for sample in samples]
        else:
            gdata["samples"] = samples

        # Initialise empty list for groups
        gdata["groups"] = dict()

        # Save genome data to JSON file
        with open(gdata_path, "w") as gdata_out:
            json.dump(gdata, gdata_out)

    # Create dummy VCF for faster IGV visualisation
    try:
        subprocess.run([
          "bcftools", "--version"],
          check=True)
    except subprocess.CalledProcessError:
        pass
    else:
        var_only = os.path.join(
          data_dir,
          genome,
          "var_only.vcf.gz")
        info_rm = subprocess.Popen(
          ["bcftools", "annotate",
           "-x", "INFO",
           vcf_path],
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        gt_rm = subprocess.Popen(
          ["bcftools", "view",
           "-G",
           "-Oz", "-o",
           var_only],
          stdin=info_rm.stdout,
          stdout=subprocess.PIPE,
          stderr=subprocess.PIPE)
        out, err = gt_rm.communicate()
        if gt_rm.returncode != 0:
            raise OSError(err.decode())
        # Index VCF
        pysam.tabix_index(
          var_only,
          preset="vcf")

    # Run raw marker set QC
    stats_dir = os.path.join(
      data_dir,
      genome,
      "stats")
    if not os.path.exists(stats_dir):
        os.mkdir(stats_dir)
    run_qc(genome)

def annotate_db(genome, gconfig):
    info_path = os.path.join(
      data_dir,
      genome,
      gconfig["gene_info"])
    info = pd.read_csv(
      info_path,
      index_col=[0],
      sep=None,
      engine="python")
    info = info.fillna("")
    info = info.astype(str)
    cols = info.columns

    db_path = os.path.join(
      gconfig["gff_db"],
      f"{genome}.db")
    db = gffutils.FeatureDB(db_path)
    new_entries = list()
    for gene in db.features_of_type("gene"):
        # GFFutils splits attributes by comma
        # Rejoin them for uniformity
        for attr, val in gene.attributes.items():
            gene.attributes[attr] = ",".join(val)
        # Save new attributes
        if gene.id in info.index:
            data = info.loc[gene.id]
            for col in cols:
                val = data[col]
                gene.attributes[col] = val if not pd.isna(val) else ""
        new_entries.append(gene)
    db = db.update(
      new_entries,
      merge_strategy="replace",
      make_backup=False)

def calc_enr(genome, gconfig, num_genes):
    info_path = os.path.join(
      data_dir,
      genome,
      gconfig["gene_info"])
    info = pd.read_csv(
      info_path,
      header=[0],
      sep=None,
      engine="python")

    enr = dict()
    enr["num_genes"] = num_genes
    cols = gconfig["enr_cols"].split(",")
    for col in cols:
        col = col.strip()
        enr_vals = dict()
        attrs = [
          attr.strip() for row_attrs in info[col].dropna()
          for attr in row_attrs.split(",") if attr]
        enr[col] = Counter(attrs)

    enr_path = os.path.join(
      data_dir,
      genome,
      ".enr")
    with open(enr_path, "w") as enr_out:
        json.dump(enr, enr_out)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Requires a single argument: name of folder in data/ folder.") 
    else:
        try:
            configure(sys.argv[1])
            print("Finished initialising genome.")
        except:
            print("Problem initialising genome, check config.ini.")