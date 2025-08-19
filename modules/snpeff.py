#!/usr/bin/env python

import io
import os
import re
import pysam
import shutil
import subprocess
from Bio.bgzf import BgzfWriter

from modules.mem import get_mem
from app.dirs import data_dir, tool_dir


def snpeff_annotate(genome, vcf_file, gff_file, seq_file, cds_file, pep_file, task_id):

###############################################################################
### Prepare input/output
###############################################################################

    genome_safe = re.sub("\W", "_", genome)
    in_root = os.path.join(data_dir, genome)
    # snpEff requires specific filenames for database building
    snpeff_root = os.path.join("snpEff", "data", genome_safe)
    seq_renamed = seq_intmd = "sequences.fa"
    gff_renamed = gff_intmd = "genes"
    cds_renamed = cds_intmd = "cds.fa"
    pep_renamed = pep_intmd = "protein.fa"

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

    print("<b>Initializing snpEff database</b>")
    with open(os.path.join("snpEff", "snpEff.config"), "w") as database:
        database.write(f"\n# {genome_safe}, created by user via octava\n")
        database.write(f"{genome_safe}.genome : {genome_safe}\n")

###############################################################################
### Copy input files
###############################################################################

    print("<b>Initializing files</b>")
    for source, target in (
        (gff_file, gff_intmd),
        (seq_file, seq_intmd),
        (cds_file, cds_intmd),
        (pep_file, pep_intmd),
    ):
        if source:
            shutil.copy(
                os.path.join(in_root, source), os.path.join(snpeff_root, target)
            )
            if target.endswith(".gz"):
                subprocess.run(
                    ["gunzip", "-f", os.path.join(snpeff_root, target)], check=True
                )
    if os.path.exists(os.path.join(snpeff_root, vcf_file)):
        os.remove(os.path.join(snpeff_root, vcf_file))

    # Link VCF or copy if cross-device
    try:
        os.link(
            os.path.join(in_root, vcf_file),
            os.path.join(snpeff_root, vcf_file)
        )
    except OSError as e:
        if e.errno == 18:
            shutil.copy2(
                os.path.join(in_root, vcf_file),
                os.path.join(snpeff_root, vcf_file)
            )
        else:
            raise

###############################################################################
### Create snpEff database
###############################################################################

    print("<b>Building database</b>")
    if gff_ext == ".gff":
        gff_option = "-gff3"
    elif gff_ext == ".gtf":
        gff_option = "-gtf22"
    else:
        raise OSError("Invalid gff_file extension.")
    cmd = [
        "java",
        "-jar",
        f"-Xms{get_mem()}m",
        f"-Xmx{get_mem()}m",
        os.path.join(tool_dir, "snpEff", "snpEff.jar"),
        "build",
        gff_option,
        "-c",
        os.path.join("snpEff", "snpEff.config"),
        "-v",
        genome_safe,
    ]
    if not cds_file:
        cmd.append("-noCheckCds")
    if not pep_file:
        cmd.append("-noCheckProtein")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

###############################################################################
### Annotate and index VCF file
###############################################################################

    print("<b>Annotating VCF file</b>")
    print(" ".join([
            "java",
            "-jar",
            f"-Xms{get_mem()}m",
            f"-Xmx{get_mem()}m",
            os.path.join(tool_dir, "snpEff", "snpEff.jar"),
            genome_safe,
            "-c",
            os.path.join("snpEff", "snpEff.config"),
            os.path.join(snpeff_root, vcf_file),
        ])
    pr = subprocess.Popen(
        [
            "java",
            "-jar",
            f"-Xms{get_mem()}m",
            f"-Xmx{get_mem()}m",
            os.path.join(tool_dir, "snpEff", "snpEff.jar"),
            genome_safe,
            "-c",
            os.path.join("snpEff", "snpEff.config"),
            os.path.join(snpeff_root, vcf_file),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    ann_path = os.path.join(in_root, f"ann.{vcf_file}")
    ann_out = BgzfWriter(ann_path, "w")
    for line in pr.stdout:
        ann_out.write(line)
    pr.wait()

    print("<b>Indexing VCF file</b>")
    pysam.tabix_index(ann_path, preset="vcf", force=True)

    # Remove copies of genome and annotation files
    shutil.rmtree(snpeff_root)
