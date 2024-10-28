#!/usr/bin/env python

import os
from flask import send_from_directory

import app.path

from __main__ import app

from app.dirs import (
  tmp_dir,
  data_dir)

# Favicon
@app.route("/favicon.ico")
def favicon():
    return send_from_directory(
      "static",
      "favicon.ico",
      mimetype="image/vnd.microsoft.icon")

# Index page genome icon
@app.route("/icon/<fn>")
def serve_icon(fn):
    icon_dir = os.path.join(
      "static",
      "icons")

    if not os.path.exists(
      os.path.join(
        icon_dir, 
        f"{fn}.png")
    ):
        fn = "Icon_default.png"
    else:
        fn += ".png"

    return send_from_directory(
      icon_dir,
      fn)

# General image, e.g. GWAS static plots
@app.route("/image/<path:fn>")
def serve_image(fn):
    return send_from_directory(
      "",
      fn)

# General file, e.g. tabular data, IGV.js data
@app.route("/file/<genome>/<fn>")
def serve_file(genome, fn):
    fn_dir = os.path.join(
      data_dir,
      genome)

    if fn.endswith(".gz"):
        return send_from_directory(
          fn_dir,
          fn,
          mimetype="application/gzip")

    return send_from_directory(
      fn_dir,
      fn)

# GWAS output tabular data for locus fine-mapping
@app.route("/assoc/<genome>/<markers>/<fn>")
def serve_assoc(genome, markers, fn):
    fn_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers,
      "res")

    if fn.endswith(".gz"):
        return send_from_directory(
          fn_dir,
          fn,
          mimetype="application/gzip")

    return send_from_directory(
      fn_dir,
      fn)

# LD file for locus fine-mapping
@app.route("/ld_file/<task_id>")
def serve_ld(task_id):
    return send_from_directory(
      tmp_dir,
      f"{task_id}.r2")

# PIP file for locus fine-mapping
@app.route("/pip_file/<task_id>")
def serve_pip(task_id):
    return send_from_directory(
      tmp_dir,
      f"{task_id}.fm")

