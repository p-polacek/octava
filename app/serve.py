#!/usr/bin/env python

import os
from flask import send_from_directory

import app.path

from flask_app import app

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
@app.route("/icon/<genome>")
def serve_icon(genome):
    if genome == "add_new":
        icon_dir = os.path.join(
          "static",
          "icons")
        return send_from_directory(
          icon_dir,
          "add_new.png")

    icon_dir = os.path.join(
      data_dir,
      genome)
    for ext in (".jpg", ".jpeg", ".png"):
        if os.path.exists(os.path.join(icon_dir, "icon" + ext)):
            fn = "icon" + ext
            break
    else:
        icon_dir = os.path.join(
          "static",
          "icons")
        fn = "icon_default.png"

    return send_from_directory(
      icon_dir,
      fn)

# Index page genome icon
@app.route("/button/<button>")
def serve_button(button):
    icon_dir = os.path.join(
      "static",
      "icons")
    return send_from_directory(
      icon_dir,
      f"{button}.svg")

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

# QTLs table
@app.route("/qtl_table/<task_id>")
def serve_qtl(task_id):
    return send_from_directory(
      tmp_dir,
      "QTL_data.xlsx")

