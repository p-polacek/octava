#!/usr/bin/env python3

genome_dir = "data"
snpEff_url = "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"


import os
import re
import pysam
import subprocess

from app.dirs import data_dir
from app.task_worker import progress

from modules.assoc import popen

def snpeff_annotate(
  genome,
  vcf_file,
  gff_file,
  seq_file,
  cds_file,
  pep_file,
  task_id):

###############################################################################
### Download snpEff
###############################################################################

    if not os.path.exists("snpEff/"):
        progress(task_id, "Downloading snpEff")
        popen(
          task_id,
          ["wget",
           snpEff_url])
        popen(
          task_id,
          ["unzip",
           "-q",
           "snpEff_latest_core.zip"])
        popen(
          task_id,
          ["rm",
           "snpEff_latest_core.zip"])

###############################################################################
### Prepare input/output
###############################################################################

    in_root = os.path.join(
      data_dir,
      genome)
    # snpEff requires specific filenames for database building
    snpeff_root = os.path.join(
      "snpEff",
      "data",
      genome)
    seq_renamed = "sequences.fa"
    gff_renamed = "genes"
    cds_renamed = "cds.fa"
    pep_renamed = "protein.fa"

    if os.path.splitext(seq_file)[1] == ".gz":
        seq_renamed += ".gz"
    if cds_file:
        if os.path.splitext(cds_file)[1] == ".gz":
            cds_renamed += ".gz"
    if pep_file:
        if os.path.splitext(pep_file)[1] == ".gz":
            pep_renamed += ".gz"

    gff_root, gff_ext = os.path.splitext(re.sub("\.gz$", "", gff_file))
    gff_renamed += gff_ext
    if gff_ext == ".gff3":
        gff_renamed, gff_ext = gff_renamed[:-1], gff_ext[:-1]
    if gff_file.endswith(".gz"):
        gff_renamed += ".gz"

    if not os.path.exists(snpeff_root):
        os.makedirs(snpeff_root)

    with open(os.path.join("snpEff", "snpEff.config"), "r+") as database:
        for line in database.readlines():
            # Skip if database with given name is already present
            if line.startswith(genome + ".genome"):
                break
        else:
            progress(task_id, "Initializing snpEff database")
            database.write(f"\n# {genome}, created by user via octava\n")
            database.write(f"{genome}.genome : {genome}\n")

###############################################################################
### Copy input files
###############################################################################

    progress(task_id, "Initializing files")
    for source, target in (
      (gff_file, gff_renamed),
      (seq_file, seq_renamed),
      (cds_file, cds_renamed),
      (pep_file, pep_renamed)):
        if source:
            popen(
              task_id,
              ["cp",
               os.path.join(in_root, source), 
               os.path.join(snpeff_root, target)])

###############################################################################
### Create snpEff database
###############################################################################

    progress(task_id, "Building database")
    if gff_ext == ".gff":
        gff_option = "-gff3"
    elif gff_ext == ".gtf":
        gff_option = "-gtf22"
    else:
        raise OSError("Invalid gff_file extension.")
    cmd = [
      "snpEff",
      "build", gff_option,
      "-c", os.path.join("snpEff", "snpEff.config"),
      "-v", genome]
    if not cds_file:
        cmd.append("-noCheckCds")
    if not pep_file:
        cmd.append("-noCheckProtein")
    popen(
      task_id,
      cmd)

###############################################################################
### Annotate and index VCF file
###############################################################################

    progress(task_id, "Annotating VCF")
    pr = subprocess.Popen(
      ["snpEff",
       db_name,
       "-c", os.path.join("snpEff", "snpEff.config"),
       os.path.join(genome_dir, db_name, vcf_file)],
      stdout=subprocess.PIPE)
    ann, _ = pr.communicate()

    ann_path = os.path.join(
      in_root,
      f"ann.{vcf_file}")
    with open(ann_path, "wb") as ann_out:
        with io.BytesIO(ann) as ann_io:
            pysam.bgzip(ann_io, ann_out)
    progress(task_id, "Indexing VCF")
    pysam.tabix_index(
      ann_path,
      preset="vcf",
      force=True)

    # Remove copies of genome and annotation files
    # subprocess.run(f"rm snpEff/data/{db_name}/{gff_renamed}", shell=True)
    # subprocess.run(f"rm snpEff/data/{db_name}/{seq_renamed}", shell=True)
