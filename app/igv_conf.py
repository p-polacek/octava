#!/usr/bin/env python

import os
import json
import configparser
import pandas as pd

from flask import url_for

from app.serve import (
  serve_file,
  serve_stats_file)

from modules.load_input import load_config

from app.dirs import data_dir, tmp_dir

def igv_am(genome):
    gconfig = load_config(genome)

    seq_fn  = gconfig["seq"]
    seq_idx = seq_fn + ".fai"
    seq_gzi = seq_fn + ".gzi"
    gff_fn  = gconfig["ann"]
    gff_idx = gff_fn + ".tbi"
    vcf_fn  = "var_only.vcf.gz"
    vcf_idx = "var_only.vcf.gz.tbi"

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
        "compressedIndexURL": url_for(
          "serve_file",
          genome = genome,
          fn     = seq_gzi),
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
              #"visibilityWindow": gconfig["igv_vis"],
          },
          {
            "name": "Variants",
              "type": "variant",
              "format": "vcf",
              "url": url_for(
                "serve_stats_file",
                genome = genome,
                fn     = vcf_fn),
              "indexURL": url_for(
                "serve_stats_file",
                genome = genome,
                fn     = vcf_idx),
              #"visibilityWindow": gconfig["igv_vis"],
              #"height": gconfig["igv_h"],
          },
        ],
      },
      "showCursorGuide": False,
      "minimumBases": 2,
    }

    bam_fns = json.loads(gconfig["bam"])
    if bam_fns:
        for bam_fn in bam_fns:
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
              #"visibilityWindow": gconfig["igv_vis"],
              #"height": gconfig["igv_h"],
            }
            igv_config["reference"]["tracks"].append(bam_track)

    return igv_config
