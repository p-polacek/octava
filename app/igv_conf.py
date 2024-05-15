#!/usr/bin/env python

import os
import configparser

from flask import url_for
from app.serve import serve_file

from modules.load_input import load_config

from app.dirs import data_dir

def igv_am(genome):
    gconfig = load_config(genome)

    seq_fn  = gconfig["seq"]
    seq_idx = seq_fn + ".fai"
    gff_fn  = gconfig["ann"]
    gff_idx = gff_fn + ".tbi"

    if os.path.exists(
      os.path.join(
        data_dir,
        genome,
        "var_only.vcf.gz")):
        vcf_fn  = "var_only.vcf.gz"
        vcf_idx = "var_only.vcf.gz.tbi"
    else:
        vcf_fn  = gconfig["vcf"]
        vcf_idx = vcf_fn + ".tbi"

    igv_config = {
      "reference": {
        "id": genome,
        "fastaURL":  url_for(
          "serve_file",
          genome = genome,
          fn     = seq_fn),
        "indexURL":  url_for(
          "serve_file",
          genome = genome,
          fn     = seq_idx),
        "tracks": [
          {
            "name": "Genes",
            "type": "annotation",
            "format": "gff3",
            "url": url_for(
              "serve_file",
              genome = genome,
              fn     = gff_fn),
            "indexURL": url_for(
              "serve_file",
              genome = genome,
              fn     = gff_idx),
              "visibilityWindow": gconfig["vis_win"],
          },
          {
            "name": "Variants",
              "type": "variant",
              "format": "vcf",
              "url": url_for(
                "serve_file",
                genome = genome,
                fn     = vcf_fn),
              "indexURL": url_for(
                "serve_file",
                genome = genome,
                fn     = vcf_idx),
              "visibilityWindow": gconfig["vis_win"],
              "height": gconfig["height"],
          },
        ],
      },
    }

    if gconfig["bam"]:
        for bam_fn in gconfig["bam"].split(", "):
            bam_idx = f"{bam_fn}.bai"
            bam_track = {
            "name": bam_fn,
              "type": "alignment",
              "format": "bam",
              "url": url_for(
                "serve_file",
                genome = genome,
                fn     = bam_fn),
              "indexURL": url_for(
                "serve_file",
                genome = genome,
                fn     = bam_idx),
              #"visibilityWindow": gconfig["vis_win"],
              "height": gconfig["height"],
            }
            igv_config["reference"]["tracks"].append(bam_track)

    return igv_config

def igv_fm(
  genome,
  markers,
  model,
  trait,
  coord,
  pval,
  opt_hash):
    gconfig = load_config(genome)

    seq_fn  = gconfig["seq"]
    seq_idx = seq_fn + ".fai"
    gff_fn  = gconfig["ann"]
    gff_idx = gff_fn + ".tbi"
    vcf_fn  = gconfig["vcf"]
    vcf_idx = vcf_fn + ".tbi"

    gwas_fn = f"{model}.{trait}.tsv.gz"
    gwas_idx = gwas_fn + ".tbi"

    igv_config = {
      "reference": {
        "id": genome,
        "fastaURL":  url_for(
          "serve_file",
          genome = genome,
          fn     = seq_fn),
        "indexURL":  url_for(
          "serve_file",
          genome = genome,
          fn     = seq_idx),
        "tracks": [
          {
            "name": "Genes",
            "type": "annotation",
            "format": "gff3",
            "url": url_for(
              "serve_file",
              genome = genome,
              fn     = gff_fn),
            "indexURL": url_for(
              "serve_file",
              genome = genome,
              fn     = gff_idx),
            "visibilityWindow": gconfig["vis_win"],
            "order": 0,
          },
          {
            "name": "GWAS",
            "type": "gwas",
            "format": "gwas",
            "url": url_for(
              "serve_assoc",
              genome  = genome,
              markers = markers,
              fn      = gwas_fn),
            "indexURL": url_for(
              "serve_assoc",
              genome  = genome,
              markers = markers,
              fn      = gwas_idx),
            "visibilityWindow": gconfig["vis_win"],
            "height": gconfig["height"],
            "columns": {
              "chromosome": 1,
              "position": 2,
              "value": 3
            },
            "min": -0.5,
            "max": pval,
            "dotSize": 3,
            "cursor": "pointer",
            "color": "#3a8f1d",
          },
          {
            "name": "PIP",
            "type": "gwas",
            "format": "gwas",
            "posteriorProbability": True,
            "url": url_for(
              "serve_pip",
              opt_hash = opt_hash),
            "visibilityWindow": gconfig["vis_win"],
            "height": gconfig["height"],
            "columns": {
              "chromosome": 1,
              "position": 2,
              "value": 3
            },
            "min": -0.1,
            "max": 1.1,
            "dotSize": 3,
            "color": "#3a8f1d",
          },
        ],
     },
     "locus": coord
    }

    return igv_config