#!/usr/bin/env python

import os
import gzip
import json
import math
import shutil
import hashlib
import numpy as np
import pandas as pd
from collections import Counter
from scipy.stats import hypergeom

from flask import (
  jsonify,
  request,
  url_for,
  redirect,
  render_template)

import app.path

from __main__ import app

from modules.forms import AssocResForm
from modules.vcf_reader import VcfReader

from modules.load_input import (
  load_gdata,
  load_config)


from app.dirs import (
  tmp_dir,
  data_dir)


# Gene enrichment data
@app.route("/enr", methods=["POST"])
def fetch_enr():
    data = request.form
    genome = data.get("genome")

    col = data.get("col")
    attrs = json.loads(data.get("attrs"))
    num_subset = len(attrs)
    attrs = [
      attr.strip() for row_attrs in attrs 
      for attr in row_attrs.split(",") if attr]
    num_tests = len(attrs)

    enr_path = os.path.join(
      data_dir,
      genome,
      ".enr")
    with open(enr_path, "r") as enr_in:
        enr = json.load(enr_in)
        enr, num_genes = enr[col], enr["num_genes"]
    enr_res = [
      [attr,
       np.log2((cnt / num_subset) / (enr[attr] / num_genes)),
       -np.log10(hypergeom(num_genes, enr[attr], num_subset).sf(cnt - 1)),
       enr[attr],
       cnt]
      for attr, cnt in Counter(attrs).items()]

    enr_res = sorted(enr_res, key=lambda x: x[2], reverse=True)
    return jsonify({
      "data": enr_res,
      "num_genes": num_genes,
      "num_subset": num_subset})

# Per-chromosome LD
@app.route("/chr_ld", methods=["POST"])
def fetch_chr_ld():
    form = request.form
    genome  = form["genome"]
    markers = form["markers"]
    chrom   = form["chrom"]

    if markers:
        stats_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "stats")
    else:
        stats_dir = os.path.join(
          data_dir,
          genome,
          "stats")
    r2_path = os.path.join(
      stats_dir,
      f"{chrom}.r2")
    with gzip.open(r2_path, "rt") as r2_in:
        chr_r2 = json.load(r2_in)
    return jsonify(chr_r2)

# Per-marker LD for locus fine-mapping
@app.route("/ld/<opt_hash>/<genome>/<markers>/<model>/<ch>/<pos>")
def fetch_ld(
  opt_hash,
  genome,
  markers,
  model,
  ch,
  pos):
    in_ld = os.path.join(
      tmp_dir,
      f"{opt_hash}.vcor2")
    in_var = os.path.join(
      tmp_dir,
      f"{opt_hash}.vcor2.vars")
    out_hash = hashlib.sha256(
                 str(
                   [markers,
                    model,
                    ch,
                    pos]).encode()
                  ).hexdigest()
    out_ld = os.path.join(
      tmp_dir,
      f"{out_hash}.vcor2")

    var = pd.read_csv(
      in_var,
      sep="\t",
      header=None)[0]
    var_idx = var[var == f"{ch}:{pos}"].index.values[0]

    ld = pd.read_csv(
      in_ld,
      sep="\t",
      header=None)

    ld = ld[[var_idx]].set_axis(
      ["R2"],
      axis=1)
    ld[["Chr", "Pos"]] = var.str.split(":", expand=True)
    ld = ld[["Chr", "Pos", "R2"]]
    ld.to_csv(
      out_ld,
      sep="\t",
      index=False)

    return jsonify(out_hash)

# Per-metric table data for GWAS QC filtering
@app.route("/filt", methods=["POST"])
def fetch_filt():
    form      = request.form
    genome    = form["genome"]
    metric    = form["metric"]
    start_row = int(form["start_row"])
    end_row   = int(form["end_row"])

    stats_root = os.path.join(
      data_dir,
      genome,
      "stats",
      "gt")

    fn = f"{stats_root}.{metric}.csv.gz"
    with gzip.open(fn, "rt") as metr_in:
        for i in range(1, start_row):
            next(metr_in)
        lines = [
          next(metr_in).strip().split(",")
          for l in range(start_row, end_row + 1)]
    return jsonify(lines)

# Filter histograms per marker set
@app.route("/hist", methods=["POST"])
def fetch_hist():
    form = request.form
    genome  = form["genome"]
    markers = form["markers"]

    stats_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers,
      "stats")
    hist_data = dict()
    for metric in (
      "smiss",
      "vmiss",
      "shet",
      "vhet",
      "maf",
      "hwe",
      "dens",
      "pca",
      "kin"):
        hist_path = os.path.join(stats_dir, metric)
        with open(hist_path, "r") as hist_in:
            hist_data[metric] = json.load(hist_in)

    with open(os.path.join(stats_dir, "cnts"), "r") as cnts_in:
        samp_num, var_num = json.load(cnts_in)

    return jsonify([
      hist_data,
      samp_num,
      var_num])

# Manhattan plot data for GWAS visualisation
@app.route("/mplot", methods=["POST"])
def fetch_mplot():
    form     = request.form
    genome   = form["genome"]
    markers  = form["markers"]
    model    = form["model"]
    trait    = form["trait"]

    assoc_fn   = f"{model}.{trait}.feather"
    assoc_path = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers,
      "res",
      assoc_fn)
    if os.path.exists(assoc_path):
        df = pd.read_feather(assoc_path)
        df = df.drop(["Effect"], axis=1).to_json(orient="values")
        return jsonify(df)
    return jsonify(success=True)

# Bubble Manhattan plot data for GWAS visualisation
@app.route("/bblplot", methods=["POST"])
def fetch_bblplot():
    form     = request.form
    genome   = form["genome"]
    window   = int(form["window"])
    thresh   = float(form["thresh"])
    groups   = json.loads(form["keys"])
    plot_avg = int(form["plot_avg"])
    min_size = math.ceil(float(form.get("min_size", 0)))

    gdata = load_gdata(genome)
    offsets = {
      ch: offs
      for ch, offs in zip(
        gdata["chr_names"],
        gdata["chr_offsets"])}

    win_data = list()
    for group in groups:
        group_data = list()
        win_group  = list()

        for key in group:
            markers, model, trait = key
            string_key = ", ".join(key)
            assoc_fn = f"{model}.{trait}.feather"
            assoc_path = os.path.join(
              data_dir,
              genome,
              "assoc",
              markers,
              "res",
              assoc_fn)

            if not os.path.exists(assoc_path):
                continue

            df = pd.read_feather(assoc_path)
            df = df.drop(["Effect"], axis=1)
            df = df.set_axis(["Chr", "Pos", string_key], axis=1)
            df = df[df[string_key] >= thresh]
            df = df.set_index(["Chr", "Pos"])
            group_data.append(df[string_key])

        if not group_data:
            win_data.append(win_group)
            continue
           
        group_df = pd.concat(group_data, axis=1)

        for ch, data, in group_df.groupby(level=["Chr"]):
            data.index = data.index.droplevel()
            data = data.sort_index()

            spans = data.index.map(lambda x: (x // window) * window).unique()
            spans = spans.union((spans + window)).unique()

            for i in spans:
                chunk = data[data.index.slice_indexer(i - window, i + window)]
                if not chunk.empty:
                    chunk = chunk.dropna(axis=1, how="all")
                    if chunk.shape[1] >= min_size:
                        coord = f"{ch}:{i-window:,}-{i+window:,}"
                        if plot_avg:
                            win_group.append(
                              [int(chunk.idxmax().mean()) + offsets[ch],
                               chunk.max().mean(),
                               chunk.shape[1],
                               coord,
                               chunk.columns.tolist()])
                        else:
                            agg_chunk = chunk.max(axis=1)
                            win_group.append(
                              [int(agg_chunk.idxmax()) + offsets[ch],
                               agg_chunk.max(),
                               chunk.shape[1],
                               coord,
                               chunk.columns.tolist()])
        win_data.append(win_group)
    return jsonify(win_data)

# Various types of data on mouse right-click
@app.route("/contextmenu", methods=["POST"])
def contextmenu():
    genome = request.form["genome"]
    action = request.form["action"]
    gconfig = load_config(genome)

    if action == "vars":
        genome = request.form["genome"]
        chrom  = request.form["chrom"]
        pos    = int(request.form["pos"])
        field  = request.form["field"]

        # Get the requested field for the given variant position
        vcf_path = os.path.join(
          data_dir,
          genome,
          gconfig["vcf"])
        vcf_reader = VcfReader(vcf_path)
        vcf_reader.set_region(chrom, pos - 1, pos)

        data = list()
        keys = list()
        for record in vcf_reader:
            if record.POS != pos:
                continue
            if field == "REF":
                keys += [["REF"]]
                data += [[record.REF]]
            elif field == "CSQ":
                csq   = record.CSQ
                keys  = list(csq[0].keys())
                data += [[ann[key] for key in keys] for ann in csq]
            elif field == "SEV":
                csq   = record.parsed_info_fields()["SEV"]
                desc  = vcf_reader.metadata["INFO"]["SEV"][0]["Description"]
                desc  = desc.split("'")[1].split("'")[0]
                keys  = [label.strip() for label in desc.split("|")]
                data += [ann.split("|") for ann in csq]
            else:
                gdata = load_gdata(genome)
                field = [
                  s_id for (s_id, s_name) in gdata["samples"]
                  if s_name == field][0]
                qual  = record.get_sample_call(field)
                keys  = list(qual.keys())
                data += [[qual[key] for key in keys]]

        return jsonify([data, keys])

    else:
        if action == "markers":
            params_json = os.path.join(
              data_dir,
              genome,
              "assoc",
              request.form["markers"],
              "stats",
              "params")

        elif action == "gwas_run":
            markers = request.form["markers"]
            model   = request.form["model"]
            trait   = request.form["trait"]
            fn = f"{model}.{trait}.params"
            params_json = os.path.join(
              data_dir,
              genome,
              "assoc",
              markers,
              "res",
              fn)
        try:
            with open(params_json, "r") as params_in:
                params = json.load(params_in)
            # Last element is phenotype variance kept for finemapping
            # Skipped here
            return jsonify(
            [[col[:-1] for col in params[1:]],
             params[0][:-1]])

        except:
            return jsonify(
            [[["Error fetching parameters."]],
             []])

# Append phenotype data to phenotype file
@app.route("/phe/<genome>", methods=["POST"])
def add_phe(genome):
    gconfig = load_config(genome)
    phe_path = os.path.join(
      data_dir,
      genome,
      gconfig["phe"])
    try:
        phe = pd.read_csv(
          phe_path,
          header    = [0],
          index_col = [0],
          sep       = "\t")
    except pd.errors.EmptyDataError:
        phe = pd.DataFrame()

    gdata = load_gdata(genome)
    errors = list()
    for file in request.files.getlist("new_phe"):
        fn = file.filename
        new_phe = pd.read_csv(
          file,
          header    = [0],
          index_col = [0],
          sep       = "," if fn.endswith(".csv") else "\t")

        col_dup = new_phe.columns.intersection(phe.columns)
        if not col_dup.empty:
            col_dup = ", ".join(col_dup.tolist())
            errors.append(f"<b>Duplicate traits in {fn}:</b><br>{col_dup}")
        else:
            phe = pd.concat([phe, new_phe], axis=1)

    if not errors:
        phe.to_csv(
          phe_path,
          sep="\t")
        errors = 0

    # Return back to GWAS options page
    form = AssocResForm(request.form)
    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")

    form.markers.choices = [
      d for d in os.listdir(markers_root)
      if os.path.isdir(os.path.join(markers_root, d))]

    form.traits.choices = phe.columns.tolist()
    form.chrom.choices = [""] + gdata["chr_names"]

    form.trait_groups.choices = [
      [json.dumps(group[0]), group[1]]
      for group in gdata["groups"].get("traits", [])]

    return render_template(
      "assoc/opts.html",
      genome  = genome,
      form    = form,
      phe     = phe,
      phe_upl = gconfig["phe_upload"],
      errors  = errors)

@app.route("/group", methods=["POST"])
def group():
    data = json.loads(request.get_data())
    genome = data["genome"]

    gdata_path = os.path.join(
      data_dir,
      genome,
      ".gdata")

    gdata = load_gdata(genome)
    gdata["groups"].setdefault(data["set"], [])

    # If making a group from a list of items
    if "name" in data.keys():
        group_entry = (
          data["items"],
          data["name"])

        # If a group with the given name exists, return error
        for group in gdata["groups"][data["set"]]:
            if group[1] == data["name"]:
                return "Group with this name alrady exists.", 400

        gdata["groups"][data["set"]].append(group_entry)

    # If deleting one or multiple groups
    else:
        filt = [
          group for group in gdata["groups"][data["set"]]
          if group[1] not in data["items"]]
        gdata["groups"][data["set"]] = filt

    # Re-save the edited genome data
    with open(gdata_path, "w") as data_out:
        json.dump(gdata, data_out)

    return "", 204

@app.route("/rm_markers", methods=["POST"])
def rm_markers():
    data = request.form
    genome = data["genome"]
    markers = data["markers"]

    if load_config(genome)["rm_markers"]:
        markers_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers)
        shutil.rmtree(markers_dir)

    return "", 204