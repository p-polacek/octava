#!/usr/bin/env python

snpEff_url = "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"

import io
import os
import re
import pysam
import shutil
import subprocess
from Bio.bgzf import BgzfWriter

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

    genome_safe = re.sub("\W", "_", genome)
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
      genome_safe)
    seq_renamed = seq_intmd = "sequences.fa"
    gff_renamed = gff_intmd = "genes"
    cds_renamed = cds_intmd = "cds.fa"
    pep_renamed = pep_intmd = "protein.fa"
    """
    if os.path.splitext(seq_file)[1] == ".gz":
        seq_intmd += ".gz"
    if cds_file and os.path.splitext(cds_file)[1] == ".gz":
        cds_intmd += ".gz"
    if pep_file and os.path.splitext(pep_file)[1] == ".gz":
        pep_intmd += ".gz"

    gff_root, gff_ext = os.path.splitext(re.sub("\.gz$", "", gff_file))
    gff_renamed += gff_ext
    gff_intmd += gff_ext
    if gff_ext == ".gff3":
        gff_renamed, gff_ext = gff_renamed[:-1], gff_ext[:-1]
        gff_intmd = gff_renamed
    if gff_file.endswith(".gz"):
        gff_intmd += ".gz"

    if not os.path.exists(snpeff_root):
        os.makedirs(snpeff_root)

    with open(os.path.join("snpEff", "snpEff.config"), "r+") as database:
        for line in database.readlines():
            # Skip if database with given name is already present
            if line.startswith(genome_safe + ".genome"):
                break
        else:
            progress(task_id, "Initializing snpEff database")
            database.write(f"\n# {genome_safe}, created by user via octava\n")
            database.write(f"{genome_safe}.genome : {genome_safe}\n")

###############################################################################
### Copy input files
###############################################################################

    progress(task_id, "Initializing files")
    for source, target in (
      (gff_file, gff_intmd),
      (seq_file, seq_intmd),
      (cds_file, cds_intmd),
      (pep_file, pep_intmd)):
        if source:
            popen(
              task_id,
              ["cp",
               os.path.join(in_root, source), 
               os.path.join(snpeff_root, target)])
            if target.endswith(".gz"):
                popen(
                  task_id,
                  ["gunzip",
                   "-f",
                   os.path.join(snpeff_root, target)])
    if os.path.exists(os.path.join(snpeff_root, vcf_file)):
        os.remove(os.path.join(snpeff_root, vcf_file))
    os.link(
      os.path.join(in_root, vcf_file),
      os.path.join(snpeff_root, vcf_file))

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
      "-v", genome_safe]
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
    """
    progress(task_id, "Annotating VCF file")
    pr = subprocess.Popen(
      ["java", "-jar",
       "snpEff/snpEff.jar",
       genome_safe,
       "-c", os.path.join("snpEff", "snpEff.config"),
       os.path.join(snpeff_root, vcf_file)],
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True)

    ann_path = os.path.join(
      in_root,
      f"ann.{vcf_file}")
    ann_out = BgzfWriter(ann_path, "w")
    for line in pr.stdout:
        ann_out.write(line)
    pr.wait()

    progress(task_id, "Indexing VCF file")
    pysam.tabix_index(
      ann_path,
      preset="vcf",
      force=True)

    # Remove copies of genome and annotation files
    #shutil.rmtree(snpeff_root)
