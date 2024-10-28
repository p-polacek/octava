#!/usr/bin/env python

import os
import re
import sys
import json
import pysam
import shutil
import gffutils
import subprocess
import numpy as np
import configparser
import pandas as pd
from collections import Counter

from app.dirs import data_dir
from app.task_worker import progress
from app.gffutils_db import gffutils_db_path

from modules.assoc import run_qc

def genome_setup(
  genome,
  params,
  task_id):

###############################################################################
### Set up GFFUtils database
###############################################################################

    if not os.path.exists(gffutils_db_path):
        os.makedirs(gffutils_db_path)
    ann_path = os.path.join(
      data_dir,
      genome,
      params["ann"])
    db_path = os.path.join(
      gffutils_db_path,
      f"{genome}.db")
    if not os.path.exists(db_path):
        progress(task_id, "Creating an annotation database")
        gffutils.create_db(
          ann_path,
          dbfn=db_path,
          merge_strategy="merge")
    db = gffutils.FeatureDB(db_path)

###############################################################################
### Add gene attributes to GFFUtils database, calculate enrichment
###############################################################################

    if params["info"]:
        progress(task_id, "Updating database information")
        annotate_db(
          genome,
          db,
          params["info"],
          task_id)
    if params["enr"]:
        progress(task_id, "Calculating gene enrichment")
        num_genes = db.count_features_of_type("gene")
        calc_enr(
          genome,
          params["info"],
          params["enr"],
          num_genes,
          task_id)

###############################################################################
### BGZip and tabix files if necessary
###############################################################################

    for file in ("vcf", "ann", "seq", "exp"):
        file_path = os.path.join(
          data_dir,
          genome,
          params[file])
        if params[file] and not file_path.endswith(".gz"):
            if os.path.islink(file_path):
                progress(task_id, f"Copying {file} file")
                link_path = os.readlink(file_path)
                os.unlink(file_path)
                shutil.copy2(link_path, file_path)

            progress(task_id, f"Compressing {file} file")
            pysam.tabix_compress(file_path, file_path + ".gz")
            os.remove(file_path)
            params[file] += ".gz"

        file_path = os.path.join(
          data_dir,
          genome,
          params[file])
        if file == "seq":
            if not os.path.exists(file_path + ".fai") or not os.path.exists(file_path + ".gzi"):
                progress(task_id, f"Indexing {file} file")
                _ = pysam.FastaFile(file_path)
        else:
            if params[file] and not os.path.exists(file_path + ".tbi"):
                progress(task_id, f"Indexing {file} file")
                if file == "exp":
                    pysam.tabix_index(
                      file_path,
                      seq_col=1,
                      start_col=2,
                      end_col=3,
                      line_skip=1)
                else:
                    pysam.tabix_index(
                      file_path,
                      preset=("vcf" if file == "vcf" else "gff"))
    for file in params["bam"]:
        file_path = os.path.join(
          data_dir,
          genome,
          file)
        if not os.path.exists(file_path + ".bai"):
            progress(task_id, "Indexing a bam file")
            pysam.index(file)

###############################################################################
### Get genome and sample data
###############################################################################

    progress(task_id, "Parsing vcf file")
    vcf_path = os.path.join(
      data_dir,
      genome,
      params["vcf"])

    # Create a genome data file if none exists
    gdata_path = os.path.join(
      data_dir,
      genome,
      ".gdata")
#    if not os.path.exists(gdata_path):
    gdata = dict()
    # Get chromosome names, sizes and offsets
    vcf = pysam.VariantFile(vcf_path)
    chroms = list(vcf.header.contigs)
    sizes = [vcf.header.contigs[c].length for c in chroms]
    offsets = [sum(sizes[:i]) for i in range(len(sizes))]
    gdata["chr_names"] = chroms
    gdata["chr_sizes"] = sizes
    gdata["chr_offsets"] = offsets

    # Get first gene ID in the file
    for rec in db.region(
      region=str(chroms[0]),
      featuretype="gene"):
        gdata["id_fmt"] = rec.id
        break

    # Check that samples match between VCF and SEQ files
    idx_path = os.path.join(
      data_dir,
      genome,
      params["seq"] + ".fai")
    idx_df = pd.read_csv(
      idx_path,
      sep="\t",
      header=None)
    seq_chroms = idx_df.iloc[:, 0].astype(str).tolist()
    seq_sizes = idx_df.iloc[:, 1].tolist()

    if seq_chroms != chroms or seq_sizes != sizes:
        print(
          "Chromosomes do not match in size or names between "
          "vcf and seq files.")
        raise OSError

    progress(task_id, "Parsing phenotypes")
    # Save sample list, combine with phenotyped samples
    samples = list(vcf.header.samples)
    if not params["phe"]:
        params["phe"] = "phe.txt"
        phe_path = os.path.join(
          data_dir,
          genome,
          params["phe"])
        with open(phe_path, "w") as phe_out:
            phe_out.write("\n".join(["Acc"] + samples))
    else:
        phe_path = os.path.join(
          data_dir,
          genome,
          params["phe"])
        if os.path.islink(phe_path):
            link_path = os.readlink(phe_path)
            os.unlink(phe_path)
            shutil.copy2(link_path, phe_path)
        phe_df = pd.read_csv(
          phe_path,
          sep="\t",
          header=[0],
          index_col=[0])
        # Make trait names filename-friendly
        safe_cols = [
          re.sub("\W", "_", col)
          for col in phe_df.columns]
        phe_df = phe_df.set_axis(
          safe_cols,
          axis=1)
        for sample in phe_df.index:
            if sample not in samples:
                samples.append(sample)
        miss_df = pd.DataFrame(
          np.nan,
          index=[s for s in samples if s not in phe_df.index],
          columns=phe_df.columns)
        (pd
          .concat([phe_df, miss_df])
          .loc[samples]
          .rename_axis("Acc")
          .to_csv(
            phe_path,
            sep="\t"))

    # Rename samples
    if not params["samp"]:
        params["samp"] = "samp.txt"
        samp_path = os.path.join(
          data_dir,
          genome,
          params["samp"])
        (pd
          .DataFrame(
            samples,
            index=samples,
            columns=["Name"])
          .rename_axis("Acc")
          .to_csv(
            samp_path,
            sep="\t"))
    else:
        samp_path = os.path.join(
          data_dir,
          genome,
          params["samp"])
        samp_df = pd.read_csv(
          samp_path,
          sep="\t",
          header=[0],
          index_col=[0])

        miss_samples = [s for s in samples if s not in samp_df.index]
        miss_df = pd.DataFrame(
          miss_samples,
          index=miss_samples,
          columns=[samp_df.columns[0]])
        (pd
          .concat([phe_df, miss_df])
          .loc[samples]
          .rename_axis("Acc")
          .to_csv(
            samp_path,
            sep="\t"))

    # Initialise empty list for groups
    gdata["groups"] = dict()

    # Save genome data to JSON file
    with open(gdata_path, "w") as gdata_out:
        json.dump(gdata, gdata_out)

    # Save parameters as config.ini
    _ = [params.pop(v) for v in ("name", "icon", "snpeff", "sev", "csrf_token")]
    params["bam"] = json.dumps(params["bam"])
    params["enr"] = json.dumps(params["enr"])
    params["concat"] = json.dumps(params["concat"])
    params = {k: str(v) for k, v in params.items()}
    params["rm_markers"] = True

    config = configparser.ConfigParser()
    config["setup"] = params
    config_path = os.path.join(
      data_dir,
      genome,
      "config.ini")
    with open(config_path, "w") as config_out:
        config.write(config_out)

    # Run raw marker set QC
    stats_dir = os.path.join(
      data_dir,
      genome,
      "stats")
    if not os.path.exists(stats_dir):
        os.mkdir(stats_dir)
    run_qc(
      genome,
      task_id)
    assoc_dir = os.path.join(
      data_dir,
      genome,
      "assoc")
    if not os.path.exists(assoc_dir):
        os.mkdir(assoc_dir)

    init_path = os.path.join(
      data_dir,
      genome,
      "init")
    open(init_path, "w").close()

def annotate_db(
  genome,
  db,
  info,
  task_id):
    info_path = os.path.join(
      data_dir,
      genome,
      info)
    info_df = pd.read_csv(
      info_path,
      sep="\t",
      header=[0],
      index_col=[0])
    info_df = info_df.fillna("")
    info_df = info_df.astype(str)
    cols = info_df.columns

    new_entries = list()
    for gene in db.features_of_type("gene"):
        # GFFutils splits attributes by comma
        # Rejoin them for uniformity
        for attr, val in gene.attributes.items():
            gene.attributes[attr] = ",".join(val)
        # Save new attributes
        if gene.id in info_df.index:
            data = info_df.loc[gene.id]
            for col in cols:
                val = data[col]
                gene.attributes[col] = val# if not pd.isna(val) else ""
        else:
            for col in cols:
                gene.attributes[col] = ""
        new_entries.append(gene)
    db = db.update(
      new_entries,
      merge_strategy="replace",
      make_backup=False)

def calc_enr(
  genome,
  info,
  enr_cols,
  num_genes,
  task_id):
    info_path = os.path.join(
      data_dir,
      genome,
      info)
    info_df = pd.read_csv(
      info_path,
      sep="\t",
      header=[0],
      index_col=[0])

    enr = dict()
    enr["num_genes"] = num_genes
    for col in enr_cols:
        enr_vals = dict()
        attrs = [
          attr.strip() for row_attrs in info_df[col].dropna()
          for attr in row_attrs.split(",") if attr]
        enr[col] = Counter(attrs)

    enr_path = os.path.join(
      data_dir,
      genome,
      ".enr")
    with open(enr_path, "w") as enr_out:
        json.dump(enr, enr_out)
