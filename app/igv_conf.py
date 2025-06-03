#!/usr/bin/env python

import os
import json
import configparser
import pandas as pd

from flask import url_for
from app.serve import serve_file

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
              "visibilityWindow": gconfig["igv_vis"],
          },
#          {
#            "name": "GWAS",
#            "type": "gwas",
#            "format": "gwas",
#            "url": url_for(
#              "serve_assoc",
#              genome  = genome,
#              markers = "Filtered",
#              fn      = "mlm.Capsaicin.tsv.gz"),
#            "indexURL": url_for(
#              "serve_assoc",
#              genome  = genome,
#              markers = "Filtered",
#              fn      = "mlm.Capsaicin.tsv.gz.tbi"),
#            "visibilityWindow": gconfig["igv_vis"],
#            "height": gconfig["igv_h"],
#            "columns": {
#              "chromosome": 1,
#              "position": 2,
#              "value": 3
#            },
#            "min": -0.5,
#            "max": 18,
#            "dotSize": 3,
#            "cursor": "pointer",
#            "color": "#3a8f1d",
#          },
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
              "visibilityWindow": gconfig["igv_vis"],
              "height": gconfig["igv_h"],
          },
        ],
      },
      "showCursorGuide": False,
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
              "height": gconfig["igv_h"],
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
  task_id):
    gconfig = load_config(genome)

    seq_fn  = gconfig["seq"]
    seq_idx = seq_fn + ".fai"
    seq_gzi = seq_fn + ".gzi"
    gff_fn  = gconfig["ann"]
    gff_idx = gff_fn + ".tbi"

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
            "visibilityWindow": gconfig["igv_vis"],
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
            "visibilityWindow": gconfig["igv_vis"],
            "height": gconfig["igv_h"],
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
            "name": "PIP-C",
            "type": "gwas",
            "format": "gwas",
            "posteriorProbability": True,
            "url": url_for(
              "serve_pip",
              task_id = task_id),
            "visibilityWindow": gconfig["igv_vis"],
            "height": gconfig["igv_h"],
            "columns": {
              "chromosome": 2,
              "position": 3,
              "value": 4
            },
            "min": -0.1,
            "max": 1.1,
            "dotSize": 3,
            "color": "#3a8f1d",
          },
          {
            "name": "PIP-RB",
            "type": "gwas",
            "format": "gwas",
            "posteriorProbability": True,
            "url": url_for(
              "serve_pip",
              task_id = task_id),
            "visibilityWindow": gconfig["igv_vis"],
            "height": gconfig["igv_h"],
            "columns": {
              "chromosome": 2,
              "position": 3,
              "value": 5
            },
            "min": -0.1,
            "max": 1.1,
            "dotSize": 3,
            "color": "#3a8f1d",
          },
        ],
      },
      "locus": coord,
      "showCursorGuide": False,
    }

    if model == "mlm":
        eff_max = (pd
          .read_csv(
            os.path.join(
              data_dir,
              genome,
              "assoc",
              markers,
              "res",
              f"{model}.{trait}.tsv.gz"),
            sep="\t",
            header=0,
            usecols=[3])
          .max().item())
        eff_track = {
          "name": "Effect",
          "type": "gwas",
          "format": "gwas",
          "posteriorProbability": True,
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
          "visibilityWindow": gconfig["igv_vis"],
          "height": gconfig["igv_h"],
          "columns": {
            "chromosome": 1,
            "position": 2,
            "value": 4
          },
          "min": -0.1,
          "max": eff_max + 1,
          "dotSize": 3,
          "cursor": "pointer",
          "color": "#3a8f1d",
        }
        igv_config["reference"]["tracks"].insert(2, eff_track)

    return igv_config


def igv_rf(
  genome,
  markers,
  model,
  trait,
  coord,
  pval,
  task_id):
    imp_max = (pd
      .read_csv(
        os.path.join(
          tmp_dir,
          f"{task_id}.fm"),
        sep="\t",
        header=0,
        usecols=[2])
      .max().item())
    gconfig = load_config(genome)

    seq_fn  = gconfig["seq"]
    seq_idx = seq_fn + ".fai"
    seq_gzi = seq_fn + ".gzi"
    gff_fn  = gconfig["ann"]
    gff_idx = gff_fn + ".tbi"

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
            "visibilityWindow": gconfig["igv_vis"],
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
            "visibilityWindow": gconfig["igv_vis"],
            "height": gconfig["igv_h"],
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
            "name": "Importance",
            "type": "gwas",
            "format": "gwas",
            "posteriorProbability": False,
            "url": url_for(
              "serve_pip",
              task_id = task_id),
            "visibilityWindow": gconfig["igv_vis"],
            "height": gconfig["igv_h"],
            "columns": {
              "chromosome": 1,
              "position": 2,
              "value": 3
            },
            "min": -0.1,
            "max": imp_max + 1,
            "dotSize": 3,
            "color": "#3a8f1d",
          },
        ],
      },
      "locus": coord,
      "showCursorGuide": False,
    }

    if model == "mlm":
        eff_max = (pd
          .read_csv(
            os.path.join(
              data_dir,
              genome,
              "assoc",
              markers,
              "res",
              f"{model}.{trait}.tsv.gz"),
            sep="\t",
            header=0,
            usecols=[3])
          .max().item())
        eff_track = {
          "name": "Effect",
          "type": "gwas",
          "format": "gwas",
          "posteriorProbability": True,
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
          "visibilityWindow": gconfig["igv_vis"],
          "height": gconfig["igv_h"],
          "columns": {
            "chromosome": 1,
            "position": 2,
            "value": 4
          },
          "min": -0.1,
          "max": eff_max + 1,
          "dotSize": 3,
          "cursor": "pointer",
          "color": "#3a8f1d",
        }
        igv_config["reference"]["tracks"].insert(2, eff_track)

    return igv_config

def igv_nofm(
  genome,
  markers,
  model,
  trait,
  coord,
  pval,
  task_id):
    gconfig = load_config(genome)

    seq_fn  = gconfig["seq"]
    seq_idx = seq_fn + ".fai"
    seq_gzi = seq_fn + ".gzi"
    gff_fn  = gconfig["ann"]
    gff_idx = gff_fn + ".tbi"

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
            "visibilityWindow": gconfig["igv_vis"],
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
            "visibilityWindow": gconfig["igv_vis"],
            "height": gconfig["igv_h"],
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
        ],
      },
      "locus": coord,
      "showCursorGuide": False,
    }

    if model == "mlm":
        eff_max = (pd
          .read_csv(
            os.path.join(
              data_dir,
              genome,
              "assoc",
              markers,
              "res",
              f"{model}.{trait}.tsv.gz"),
            sep="\t",
            header=0,
            usecols=[3])
          .max().item())
        eff_track = {
          "name": "Effect",
          "type": "gwas",
          "format": "gwas",
          "posteriorProbability": True,
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
          "visibilityWindow": gconfig["igv_vis"],
          "height": gconfig["igv_h"],
          "columns": {
            "chromosome": 1,
            "position": 2,
            "value": 4
          },
          "min": -0.1,
          "max": eff_max + 1,
          "dotSize": 3,
          "cursor": "pointer",
          "color": "#3a8f1d",
        }
        igv_config["reference"]["tracks"].insert(2, eff_track)

    return igv_config
