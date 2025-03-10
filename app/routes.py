#!/usr/bin/env python

import app.path
from flask_app import (
  app,
  get_lock,
  task_queue,
  master_lock,
  task_status)

from app.dirs import (
  tmp_dir,
  data_dir,
  root_dir)

import os
import re
import glob
import json
import uuid
import pysam
import tabix
import shutil
import functools
import subprocess
import numpy as np
import pandas as pd

from time import time
import statsmodels.api as sm
from scipy.stats import boxcox
from glimix_core.lmm import LMM
from sklearn.decomposition import PCA
from numpy_sugar.linalg import economic_qs_linear

from flask import (
  flash,
  jsonify,
  request,
  session,
  url_for,
  redirect,
  Response,
  render_template)

from app.ajax import (
  group,
  fetch_filt,
  rm_markers,
  contextmenu,
  fetch_mplot,
  fetch_bblplot)

from app.igv_conf import (
  igv_am,
  igv_fm,
  igv_rf,
  igv_nofm)

from app.serve import (
  favicon,
  serve_ld,
  serve_pip,
  serve_file,
  serve_icon,
  serve_assoc,
  serve_image)

from modules.am import am_pipeline
from modules.tmp_clean import clean
from modules.sev import annotate_sev
from modules.setup import genome_setup

from modules.assoc import (
  mdbin,
  run_qc,
  filt_qc,
  run_assoc,
  assoc_finemap)

from modules.forms import (
  AddForm,
  EditForm,
  AssocQCform,
  AssocResForm,
  TraitViewForm,
  FuncAnalysisForm)
from modules.forms import models as models_dict

models_dict = dict(models_dict)

from modules.load_input import (
  load_gdata,
  load_config)

###############################################################################
### Long task progress
###############################################################################

@app.route("/status/<task_id>", methods=["GET"])
def get_status(task_id):
    # Remove tasks finished more than a week ago
    for _id in task_status.keys():
        if (
          _id != "specs" and
          "end" in task_status[_id].keys() and
          time() - task_status[_id]["end"] > 604800):
            _ = task_status.pop(_id)

    if task_id not in task_status.keys():
        if task_id not in [t[0] for t in list(task_queue.queue)]:
            return jsonify(""), 404
        else:
            return jsonify(""), 201
    else:
        with master_lock:
            lock = get_lock("__task_status__")
        with lock:
            progress = task_status[task_id]["progress"]
            task_status[task_id]["progress"] = ""
            code = task_status[task_id]["code"]
            if code == 200 or code == 500:
                end = time()
                task_status[task_id]["end"] = end
                duration = end - task_status[task_id]["start"]
                return jsonify([progress, round(duration)]), code
            return jsonify(progress), code

###############################################################################
### Index page
###############################################################################

@app.route("/")
def index():
    # Remove old files from tmp_dir
    with master_lock:
        lock = get_lock("__clean__")
    with lock:
        clean()
    genomes = sorted(
      [p for p in os.listdir(data_dir)
       if os.path.exists(
         os.path.join(
           data_dir,
           p,
           "init"))])
    links = {
      genome: {
        "gv": url_for("fa_opts", genome=genome),
        "ts": url_for("ts", genome=genome),
        "assoc": url_for("assoc", genome=genome),
        "edit": url_for("edit", genome=genome)}
      for genome in genomes}

    if os.path.exists(".noupdate"):
        no_update = True
        os.remove(".noupdate")
    else:
        no_update = False

    return render_template(
      "index.html",
      genomes = genomes,
      links   = links,
      genome  = "",
      no_update = no_update)

###############################################################################
### Error page
###############################################################################

@app.route("/error/<error>")
def error(error):
    return render_template(
      "error.html",
      error=error)

###############################################################################
### Add new/edit dataset
###############################################################################

@app.route(
  "/new_dataset",
  methods=["GET"])
def new():
    form = AddForm()
    genomes = sorted(
      [p for p in os.listdir(data_dir)
       if os.path.isdir(os.path.join(data_dir, p)) and
       not os.path.exists(
         os.path.join(
           data_dir,
           p,
           "init"))])
    genomes = [
      (url_for("setup_new", genome=genome),
       genome)
      for genome in genomes]
    form.name.choices = [(url_for("new"), "Select dataset")] + genomes
    return render_template(
      "new.html",
      form   = form,
      genome = "",
      edit   = False,
      action = "")

@app.route(
  "/new_dataset/<genome>",
  methods=["GET"])
def setup_new(genome):
    form = AddForm()
    genomes = sorted(
      [p for p in os.listdir(data_dir)
       if os.path.isdir(os.path.join(data_dir, p)) and
       not os.path.exists(
         os.path.join(
           data_dir,
           p,
           "init"))])
    genomes = [
      (url_for("setup_new", genome=genome),
       genome)
      for genome in genomes]
    form.name.choices = [(url_for("new"), "Select dataset")] + genomes
    form.name.data = url_for("setup_new", genome=genome)

    gdir = os.path.join(data_dir, genome)
    exts = {
      "vcf": (".vcf", ".vcf.gz", ".bcf", ".bcf.gz"),
      "ann": (".gff", ".gtf", ".gff3", ".gff.gz", ".gtf.gz", ".gff3.gz")}
    exts.update(
      dict.fromkeys(
        ["seq", "cds", "pep"],
        (".fa", ".fasta", ".fna", ".faa", ".fa.gz", ".fasta.gz", ".fna.gz", ".faa.gz")))
    exts.update(
      dict.fromkeys(
        ["phe", "exp", "run", "proj", "samp", "info", "sev"],
        (".txt", ".tsv", ".txt.gz", ".tsv.gz")))
    for ftype, ext in exts.items():
        form[ftype].choices = [""] + [
          f for f in os.listdir(gdir)
          if f.endswith(ext)]
    form.bam.choices = [
          f for f in os.listdir(gdir)
          if f.endswith(".bam")]
    prot_db_path = os.path.join(gdir, "prot_db")
    if os.path.exists(prot_db_path):
        form.db_prefix.choices = [""] + list(
          set(
            match.group(1)
            for fn in os.listdir(prot_db_path)
            if (match := re.match(
              r"(.+)(\.\d+)?\.phr$", 
              fn))))
        form.db_fasta.choices = [""] + list(
          set(
            match.group(1)
            for fn in os.listdir(prot_db_path)
            if (match := re.match(
              r"((.+)\.(fasta|fa)(.gz)?)$", 
              fn))))

    return render_template(
      "new.html",
      form   = form,
      genome = genome,
      edit   = False,
      action = url_for("setup", genome=genome))

@app.route(
  "/<genome>/setup",
  methods=["POST"])
def setup(genome):
    form = AddForm(request.form)
    task_id = str(uuid.uuid4())
    task_spec = ("setup", genome)
    task = (
      genome_setup,
      {"genome": genome,
      "params": form.data,
      "task_id": task_id})

    # Make sure the same task isn't running or queued
    if task_spec in task_status["specs"]:
        return render_template(
          "error.html",
          error="Dataset is already being processed.")
    else:
        task_queue.put((task_id, task_spec, task))
        return render_template(
          "loading.html",
          genome    = genome,
          task      = "setup",
          task_id   = task_id,
          task_spec = task_spec,
          form      = "")

@app.route(
  "/<genome>/Edit",
  methods=["GET", "POST"])
def edit(genome):
    errors = []
    config = load_config(genome)
    form = EditForm()
    allow = config["edit"]

    if request.method == "POST":
        if not config["edit"]:
            errors.append("Editing is not allowed for this dataset.")
        else:
            phe_path = os.path.join(
              data_dir,
              genome,
              config["phe"])
            with master_lock:
                lock = get_lock(phe_path)
            with lock:
                phe_df = pd.read_csv(
                  phe_path,
                  sep="\t",
                  header=[0],
                  index_col=[0])

                file = request.files["phe"]
                new_phe = pd.read_csv(
                  file,
                  header=[0],
                  index_col=[0],
                  sep="," if file.filename.endswith(".csv") else "\t")

                if len(new_phe.columns) == 0:
                    errors.append("<b>No columns found in file.</b>")
                else:
                    safe_cols = [
                      re.sub("[^\w\(\)_]", "_", col)
                      for col in new_phe.columns]
                    new_phe = new_phe.set_axis(
                      safe_cols,
                      axis=1)
                    col_dup = new_phe.columns.intersection(phe_df.columns)
                    if not col_dup.empty:
                        col_dup = ", ".join(col_dup.tolist())
                        errors.append(f"<b>Duplicate traits:</b><br>{col_dup}")
                    else:
                        new_phe = pd.concat(
                          [phe_df, new_phe],
                          axis=1)
                        new_phe.to_csv(
                          phe_path,
                          sep="\t")
                        # Update samples info to contain all samples
                        samples = new_phe.index
                        samp_path = os.path.join(
                          data_dir,
                          genome,
                          config["samp"])
                        with master_lock:
                            lock_2 = get_lock(samp_path)
                        with lock_2:
                            samp_df = pd.read_csv(
                              samp_path,
                              sep="\t",
                              header=[0],
                              index_col=[0])
                            miss_samples = [
                              s for s in samples
                              if s not in samp_df.index]
                            miss_df = pd.DataFrame(
                              miss_samples,
                              index=miss_samples,
                              columns=[samp_df.columns[0]])
                            (pd
                              .concat([samp_df, miss_df])
                              .rename_axis("Acc")
                              .to_csv(
                                samp_path,
                                sep="\t"))
                            errors = 0

    return render_template(
      "new.html",
      genome  = genome,
      form    = form,
      edit    = True,
      allow   = allow,
      errors  = errors,
      action  = url_for("edit", genome=genome))

###############################################################################
### Trait/sample view
###############################################################################

@app.route(
  "/<genome>/Trait_sample_view",
  methods=["GET", "POST"])
def ts(genome):
    gdata  = load_gdata(genome)
    config = load_config(genome)
    form   = TraitViewForm()
    errors = []
    phe_path = os.path.join(
      data_dir,
      genome,
      config["phe"])

    if request.method == "POST":
        form = TraitViewForm(request.form)
        params = form.data
        with master_lock:
            lock = get_lock(phe_path)
        with lock:
            phe_df = pd.read_csv(
              phe_path,
              sep="\t",
              header=[0],
              index_col=[0])
            if params["t_type"] == "tr":
                if params["transform"] == "log":
                    transform = np.log(phe_df[[params["trait_sgl"]]])
                elif params["transform"] == "pca":
                    traits = params["traits"]
                    for group in params["traits_group"]:
                        traits += json.loads(group)
                    pca = PCA(n_components=params["num_pc"])
                    phe = phe_df[traits].dropna()
                    columns = [
                      params["t_name"] + f"_PC{i}"
                      for i in range(1, params["num_pc"] + 1)]
                    transform = (pd
                      .DataFrame(
                        pca.fit_transform(phe),
                        index=phe.index,
                        columns=columns))
                elif params["transform"] == "bin":
                    transform = ((
                      phe_df[[params["trait_sgl"]]].dropna() >= float(params["binthr"]))
                      .astype(int))
                elif params["transform"] == "zs":
                    data = phe_df[[params["trait_sgl"]]].dropna()
                    transform = ((data - np.mean(data)) / np.std(data))
                elif params["transform"] == "bc":
                    transform = phe_df[[params["trait_sgl"]]][phe_df[params["trait_sgl"]] > 0]
                    transform[params["trait_sgl"]] = boxcox(transform[params["trait_sgl"]])[0]
                elif params["transform"] == "corr":




                    kq_root = os.path.join(
                      data_dir,
                      genome,
                      "assoc",
                      params["kq_from"],
                      "gt",
                      "gt")

                    q_matr = (pd
                      .read_csv(
                        f"{kq_root}.eigenvec",
                        sep="\t",
                        header=[0],
                        index_col=[0,1])
                      .droplevel(0))
                    q_matr = q_matr.loc[:, q_matr.columns.tolist()[:int(params["q_num"])]]

                    k_matr = pd.read_csv(
                      f"{kq_root}.rel",
                      sep="\t",
                      index_col=None,
                      header=None)
                    k_lab = pd.read_csv(
                      f"{kq_root}.rel.id",
                      sep="\t",
                      usecols=[1],
                      skiprows=[0],
                      header=None,
                      index_col=None)[1]
                    k_matr = (k_matr
                      .set_axis(k_lab, axis=0)
                      .set_axis(k_lab, axis=1))

                    phe = phe_df[params["trait_sgl"]]
                    phe_samples = phe[~pd.isna(phe)].index.intersection(q_matr.index)
                    phe = phe.loc[phe_samples]

                    phe_q = q_matr.loc[phe_samples]
                    phe_q = sm.add_constant(phe_q)
                    phe_k = k_matr.loc[phe_samples, phe_samples]
                    QS = economic_qs_linear(phe_k)
                    lmm = LMM(
                      phe,
                      phe_q,
                      QS=QS)
                    lmm.fit()
                    phe_corr = phe - lmm.mean()
                    transform = pd.DataFrame(phe_corr).set_index(phe_samples)

#                    phe_q = pd.read_csv("data/Pepper 423/covariates/Pun1_deletion.txt", sep="\t", header=0, index_col=0).loc[phe_samples]
#                    phe_q = sm.add_constant(phe_q)
#                    transform = pd.DataFrame(sm.OLS(phe, phe_q).fit().resid)


                transform = transform.set_axis([params["t_name"]], axis=1)

            elif params["t_type"] == "sel":
                samples = json.loads(params["samples_sel"])
                transform = phe_df.loc[samples, [params["trait_sgl"]]]
                transform = transform.set_axis([params["t_name"]], axis=1)
            safe_cols = [
              re.sub("[^\w\(\)]", "_", col)
              for col in transform.columns]
            transform = transform.set_axis(
              safe_cols,
              axis=1)
            col_dup = transform.columns.intersection(phe_df.columns)
            if not col_dup.empty:
                col_dup = ", ".join(col_dup.tolist())
                errors.append(f"<b>Duplicate traits:</b><br>{col_dup}")
            else:
                (pd
                  .concat(
                    [phe_df, transform],
                    axis=1)
                  .to_csv(
                    phe_path,
                    sep="\t"))
                errors = 0

        flash(errors)
        return redirect(url_for("ts", genome=genome))

    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        with open(phe_path, "r") as phe_in:
            phe = [
              line.strip("\n").split("\t")
              for line in phe_in if line.strip()]
    if phe:
        form.trait_sgl.choices = phe[0][1:]
        form.traits.choices = phe[0][1:]

    vcf_path = os.path.join(
      data_dir,
      genome,
      config["vcf"])
    vcf = pysam.VariantFile(vcf_path)
    gt_samples = list(vcf.header.samples)
    samp_path = os.path.join(
      data_dir,
      genome,
      config["samp"])

    with master_lock:
        lock = get_lock(samp_path)
    with lock:
        samp_df = (pd
          .read_csv(
            samp_path,
            sep="\t",
            header=0,
            index_col=[0])
          .fillna(""))
    samp_info = [["ID"] + samp_df.columns.tolist()] + samp_df.to_records().tolist()
    samples = (samp_df[[samp_df.columns[0]]]
      .to_records()
      .tolist())
    form.samples.choices = samples
    form.samples_group.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("samples", [])]
    form.traits_group.choices = [
      [json.dumps(group[0]), group[1]]
      for group in gdata["groups"].get("traits", [])]

    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")
    for m in os.listdir(markers_root):
        init_path = os.path.join(markers_root, m, "init")
        if os.path.exists(init_path):
            form.kq_from.choices.append(m)

    return render_template(
      "trt/trt.html",
      genome     = genome,
      form       = form,
      phe        = phe,
      gt_samples = gt_samples,
      samp_info  = samp_info,
      errors     = errors)

###############################################################################
### Functional options
###############################################################################

@app.route(
  "/<genome>/Gene_variant_view",
  methods=["GET", "POST"])
def fa_opts(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = FuncAnalysisForm()

    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")
    form.dataset.choices += [
      (d, d) for d in os.listdir(markers_root)
      if os.path.isdir(os.path.join(markers_root, d))
      and os.path.exists(os.path.join(markers_root, d, "init"))]

    form.chrom.choices = gdata["chr_names"]
    samp_path = os.path.join(
      data_dir,
      genome,
      config["samp"])
    with master_lock:
        lock = get_lock(samp_path)
    with lock:
        samples = (pd
          .read_csv(
            samp_path,
            sep="\t",
            header=0,
            index_col=[0],
            usecols=[0,1])
          .to_records()
          .tolist())
    form.samples.choices = samples

    form.samples_group.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("samples", [])]
    form.gene_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("genes", [])]
    form.coord_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("coords", [])]

    for key, value in request.args.items():
        form[key].data = value

    return render_template(
      "fa/opts.html",
      genome = genome,
      form   = form,
      id_fmt = gdata["id_fmt"])

###############################################################################
### Functional analysis results
###############################################################################

@app.route(
  "/<genome>/Gene_variant_view/Results",
  methods=["GET", "POST"])
def fa_res(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = FuncAnalysisForm(request.form)
    if form["fmt"].data is None:
        saved_form = json.loads(request.form["data"])
        form.process(data=saved_form)

    if request.method == "GET":
        form.process(data=session.get(genome + "_sev_form", {}))

    form.chrom.choices = gdata["chr_names"]
    samp_path = os.path.join(
      data_dir,
      genome,
      config["samp"])
    with master_lock:
        lock = get_lock(samp_path)
    with lock:
        samples = (pd
          .read_csv(
            samp_path,
            sep="\t",
            header=0,
            index_col=[0],
            usecols=[0,1])
          .to_records()
          .tolist())
    form.samples.choices = samples

    form.samples_group.choices = [
      [json.dumps(group[0]), group[1]]
       for group in gdata["groups"].get("samples", [])]
    form.gene_groups.choices = [
      [json.dumps(group[0]), group[1]]
       for group in gdata["groups"].get("genes", [])]
    form.coord_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("coords", [])]
    form.exprun_groups.choices = [
      [json.dumps(group[0]), group[1]]
       for group in gdata["groups"].get("expruns", [])]

#    form.traits_group.choices = [
#      [json.dumps(group[0]), group[1]]
#      for group in gdata["groups"].get("traits", [])]

    (gts,
     groups,
     genes,
     anns,
     errors) = am_pipeline(genome, form)

    if config["exp"]:
        exp_path = os.path.join(
          data_dir,
          genome,
          config["exp"])

        runs_path = os.path.join(
          data_dir,
          genome,
          config["run"])
        with open(runs_path, "r") as runs_in:
            runs_info = [
              line.strip().split("\t")
              for line in runs_in]
            form.exp_col.choices = [
              ("exp-by_" + str(idx), col)
              for idx, col in enumerate(runs_info[0])]

        if config["proj"]:
            proj_path = os.path.join(
              data_dir,
              genome,
              config["proj"])
            with open(proj_path, "r") as proj_in:
                proj_info = [
                  line.strip().split("\t")
                  for line in proj_in]
        else:
            proj_info = ""

        exp_data = list()
        tb = tabix.open(exp_path)
        if genes:
            # Skip header line
            for gene in genes[1:]:
                try:
                    for gene_entry in tb.query(*gene[1:4]):
                        if gene_entry[0] != gene[0]:
                            continue
                        exp_data.append([gene_entry[0]] + gene_entry[4:])
                except tabix.TabixError:
                    # Gene not found, ignore
                    pass
    else:
        exp_data  = ""
        runs_info = ""
        proj_info = ""

    try:
        with open("test_colors.json", "r") as in_colors:
            colors = json.load(in_colors)
    except:
        colors = ""

    phe_path = os.path.join(
      data_dir,
      genome,
      config["phe"])
    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        with open(phe_path, "r") as phe_in:
            phe = [
              line.strip("\n").split("\t")
              for line in phe_in if line.strip()]
    if config["enr"]:
        form.enr_col.choices = json.loads(config["enr"])
    igv_conf = igv_am(genome)

    return render_template(
      "fa/res.html",
      genome    = genome,
      gts       = gts,
      groups    = groups,
      genes     = genes,
      anns      = anns,
      errors    = errors,
      form      = form,
      id_fmt    = gdata["id_fmt"],
      phe       = phe,
      runs_info = runs_info,
      proj_info = proj_info,
      exp_data  = exp_data,
      colors    = colors,
      igv_conf  = igv_conf)

###############################################################################
### Annotate severity
###############################################################################

@app.route(
  "/<genome>/Annotate_severity",
  methods=["POST"])
def ann_sev(genome):
    #prot_db = os.path.join(
    #  "res",
    #  "prot_db",
    #  "prot")
    prot_dir = os.path.join(
      data_dir,
      genome,
      "prot_db")
    if not os.path.exists(prot_dir):
        return render_template(
          "error.html",
          error="Missing required input files.")


    form = FuncAnalysisForm(request.form)
    markers = json.loads(form.data["ann_sev"])

    marker_data = list()
    for marker in markers:
        ch, pos = marker.split(":")
        pos = int(pos.split(".")[0])
        marker_data.append((ch, pos))

    task_id = str(uuid.uuid4())
    task_spec = ("sev", genome, marker_data)
    task = (
      annotate_sev,
      {"genome": genome,
       "markers": marker_data,
       "task_id": task_id})

    task_queue.put((task_id, task_spec, task))
    form.ann_sev.data = ""
    return render_template(
      "loading.html",
      genome    = genome,
      task      = "sev",
      task_id   = task_id,
      task_spec = task_spec,
      warning   = "",
      form      = form.data)

###############################################################################
### Association analysis redirect
###############################################################################

@app.route(
  "/<genome>/Association_analysis",
  methods=["GET"])
def assoc(genome):
    assoc_dir = os.path.join(
      data_dir,
      genome,
      "assoc")

    if not os.path.exists(assoc_dir):
        os.mkdir(assoc_dir)
    if not os.listdir(assoc_dir):
        return redirect(
          url_for(
            "assoc_qc",
            genome=genome))
    elif not any(
      os.path.exists(
        os.path.join(
          assoc_dir,
          markers,
          "res"))
      for markers in os.listdir(assoc_dir)):
        return redirect(
          url_for(
            "assoc_opts",
            genome=genome))
    else:
        return redirect(
          url_for(
            "assoc_res",
            genome=genome))

###############################################################################
### Choose filters for defining a new marker set
###############################################################################

@app.route(
  "/<genome>/Association_analysis/QC",
  methods=["GET"])
def assoc_qc(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = AssocQCform()
    stats_dir = os.path.join(
      data_dir,
      genome,
      "stats")

    hist_data = dict()
    for metric in (
      "smiss",
      "vmiss",
      "shet",
      "vhet",
      "maf",
      "hwe",
      "dens"):
        hist_path = os.path.join(
          stats_dir,
          metric)
        with open(hist_path, "r") as hist_in:
            hist_data[metric] = json.load(hist_in)

    with open(os.path.join(stats_dir, "cnts"), "r") as cnts_in:
        samp_num, var_num = json.load(cnts_in)

    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")
    saved_sets = list()
    for d in os.listdir(markers_root):
        init_path = os.path.join(
          markers_root,
          d,
          "init")
        if os.path.exists(init_path):
            saved_sets.append(d.casefold())
            params_path = os.path.join(
              markers_root,
              d,
              "stats",
              "params")
            with open(params_path) as params_in:
                params = json.load(params_in)
            if params[1][10] == "No":
                form.filt_from.choices.append(d)

    form.chr.choices = gdata["chr_names"]

    form.samples_group.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("samples", [])]

    return render_template(
      "assoc/qc.html",
      genome       = genome,
      form         = form,
      hist_data    = hist_data,
      samp_num     = samp_num,
      var_num      = var_num,
      mdbin        = mdbin,
      chr_names    = gdata["chr_names"],
      chr_sizes    = gdata["chr_sizes"],
      chr_offsets  = gdata["chr_offsets"],
      concat       = json.loads(config["concat"]),
      saved_sets   = saved_sets)

###############################################################################
### Apply filters to create new marker set, redirect
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Filter",
  methods=["POST"])
def assoc_filt(genome):
    form = AssocQCform(request.form)

    stats_root = os.path.join(
      data_dir,
      genome,
      "stats",
      "gt")

    markers = form.data["name"]
    task_id = str(uuid.uuid4())
    task_spec = ("filt_qc", genome, form.data["name"])
    task = (
      filt_qc,
      {"genome": genome,
      "data": form.data,
      "task_id": task_id})

    # Make sure the same task isn't running or queued
    if task_spec in task_status["specs"]:
        return render_template(
          "error.html",
          error="Set name unavailable. Choose another name.")
    else:
        task_queue.put((task_id, task_spec, task))
        return render_template(
          "loading.html",
          genome    = genome,
          task      = "filt",
          task_id   = task_id,
          task_spec = task_spec,
          warning   = "",
          form      = "")

###############################################################################
### Choose parameters for association analysis
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Options",
  methods=["GET"])
def assoc_opts(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = AssocResForm()

    phe_path = os.path.join(
      data_dir,
      genome,
      config["phe"])
    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        with open(phe_path, "r") as phe_in:
            phe_header = next(phe_in)
    form.traits.choices = phe_header.strip().split("\t")[1:]

    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")
    haps = dict()
    for m in os.listdir(markers_root):
        init_path = os.path.join(markers_root, m, "init")
        if os.path.exists(init_path):
            params_path = os.path.join(
              markers_root,
              m,
              "stats",
              "params")
            with open(params_path, "r") as params_in:
                hap = json.load(params_in)[1][10] == "Yes"
            haps[m] = hap

            form.markers_sgl.choices.append(m)
            form.kq_from.choices.append(m)

    form.chrom.choices = [""] + gdata["chr_names"]
    form.traits_group.choices = [
      [json.dumps(group[0]), group[1]]
      for group in gdata["groups"].get("traits", [])]
    form.samples_group.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("samples", [])]
    form.plot.choices += [
      ("kin",
       "Kinship matrix"),
      ("pca",
       "Principal components")]
    concat = json.loads(config["concat"])
    form.chr.choices = [
      c for c in gdata["chr_names"]
      if c not in concat]

    return render_template(
      "assoc/opts.html",
      genome      = genome,
      form        = form,
      phe_upl     = True,
      rm_markers  = True,
      haps        = haps,
      mdbin       = mdbin,
      chr_names   = gdata["chr_names"],
      chr_sizes   = gdata["chr_sizes"],
      chr_offsets = gdata["chr_offsets"],
      concat      = concat,
      errors      = None)

###############################################################################
### Optimize random forest parameters
###############################################################################

@app.route(
  "/Optimize_RF/<task_id>",
  methods=["GET", "POST"])
def rf_optim_view(task_id):
    optim_path = os.path.join(
      tmp_dir,
      f"{task_id}.rfopt")

    if not os.path.exists(optim_path):
        return render_template(
          "error.html",
          error="Results not found.")
    with open(optim_path) as optim_in:
        rf_data = json.load(optim_in)
    return render_template(
      "assoc/rf_optim.html",
      rf_data=rf_data)

###############################################################################
### Run association analysis
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Run",
  methods=["POST"])
def assoc_run(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    errors = list()
    form = AssocResForm(request.form)
    data = form.data

    res_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      data["markers_sgl"],
      "res")
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    error = []

    if data["model_sgl"] != "mtmlm":
        for group in data["traits_group"]:
            group, group_name = json.loads(group)
            for trait in group:
                if trait not in data["traits"]:
                    data["traits"].append(trait)

        existing_traits = list()
        for trait in data["traits"][:]:
            if not (data["overwrite"] or data["optimize"]):
                # Skip existing runs
                trait_out = os.path.join(
                  data_dir,
                  genome,
                  "assoc",
                  data["markers_sgl"],
                  "res",
                  f"{data['model_sgl']}.{trait}.feather")
                if os.path.exists(trait_out):
                    data["traits"].remove(trait)
                    existing_traits.append(trait)
        if existing_traits:
            errors.append("Results already exist for:")
            errors.append(", ".join(existing_traits))
        # Skip queued runs
        queued_traits = list()
        for task_spec in task_status["specs"]:
            if (
              task_spec[0] == "run_assoc" and
              task_spec[1] == genome and
              task_spec[2] == data["markers_sgl"] and
              task_spec[3] == data["model_sgl"]):
                for trait in task_spec[4]:
                    if trait in data["traits"]:
                        data["traits"].remove(trait)
                        queued_traits.append(trait)
        if queued_traits:
            errors.append("Results are already being generated for:")
            errors.append(", ".join(queued_traits))

        if data["traits"]:
            data["traits"] = sorted(data["traits"])

            task_id = str(uuid.uuid4())
            task_spec = (
              "run_assoc",
              genome,
              data["markers_sgl"],
              data["model_sgl"],
              data["traits"])
            task = (
              run_assoc,
              {"genome": genome,
              "params": data,
              "task_id": task_id})

            task_queue.put((task_id, task_spec, task))
            return render_template(
              "loading.html",
              genome    = genome,
              task      = "assoc" if not data["optimize"] else "rf_optim",
              task_id   = task_id,
              task_spec = task_spec,
              errors    = errors,
              form      = data)

        else:
            errors = ["Results already exist for all selected traits."]
            flash(errors)
            return redirect(url_for("assoc_res", genome=genome))
    else:
        if data["traits"]:
            pass
            #errors = ["To use multi-trait MLM, traits must be grouped."]
        else:
            data["traits_group"] = [
              json.loads(g)
              for g in data["traits_group"]]
            existing_groups = list()
            for group, group_name in data["traits_group"][:]:
                # Skip existing runs
                trait_out = os.path.join(
                  data_dir,
                  genome,
                  "assoc",
                  data["markers_sgl"],
                  "res",
                  f"{data['model_sgl']}.{group_name}.feather")
                if not data["overwrite"]:
                    if os.path.exists(trait_out):
                        data["traits_group"].remove([group, group_name])
                        existing_groups.append(group_name)
            if existing_groups:
                errors.append("Results already exist for:")
                errors.append(", ".join(existing_groups))
            # Skip queued runs
            queued_groups = list()
            for task_spec in task_status["specs"]:
                if (
                  task_spec[0] == "run_assoc" and
                  task_spec[1] == genome and
                  task_spec[2] == data["markers_sgl"] and
                  task_spec[3] == data["model_sgl"]):
                    for group in task_spec[4]:
                        if group in data["traits_group"]:
                            data["traits_group"].remove(group)
                            queued_groups.append(group[1])
            if queued_groups:
                errors.append("Results are already being generated for:")
                errors.append(", ".join(queued_groups))

            if data["traits_group"]:
                task_id = str(uuid.uuid4())
                task_spec = (
                  "run_assoc",
                  genome,
                  data["markers_sgl"],
                  data["model_sgl"],
                  data["traits_group"])
                task = (
                  run_assoc,
                  {"genome": genome,
                  "params": data,
                  "task_id": task_id})

                task_queue.put((task_id, task_spec, task))
                return render_template(
                  "loading.html",
                  genome    = genome,
                  task      = "assoc",
                  task_id   = task_id,
                  task_spec = task_spec,
                  errors    = errors,
                  form      = data)
            else:
                errors = ["Results already exist for all selected groups."]
                flash(errors)
                return redirect(url_for("assoc_res", genome=genome))

###############################################################################
### Association analysis results
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Results",
  methods=["GET", "POST"])
def assoc_res(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    errors = list()
    form = AssocResForm()
    saved_form = None
    if request.method == "POST":
        saved_form = json.loads(request.form["data"])
        form.process(data=saved_form)
    if saved_form and not form["markers"].data:
        form["markers"].data = [form["markers_sgl"].data]
        form["models"].data = [form["model_sgl"].data]
        form["traits"].data = form["traits"].data
        form["traits_group"].data = [
          json.dumps(json.loads(group.replace("'", '"'))[0])
          for group in form["traits_group"].data]

        
    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")

    form.markers.choices = [
      d for d in os.listdir(markers_root)
      if os.path.isdir(os.path.join(markers_root, d))
      and os.path.exists(os.path.join(markers_root, d, "init"))]

    analysed_traits = list()
    analysed_models = list()
    for marker_set in form.markers.choices:
        assoc_root = os.path.join(
          data_dir,
          genome,
          "assoc",
          marker_set,
          "res")
        if not os.path.exists(assoc_root):
            continue

        for fn in os.listdir(assoc_root):
            if not fn.endswith(".feather"):
                continue

            fn_data = fn.rsplit(".", 1)[0]
            model, trait = fn_data.split(".", 1)
            if trait not in analysed_traits:
                analysed_traits.append(trait)
            if model not in analysed_models:
                analysed_models.append(model)

    form.traits.choices = analysed_traits
    form.models.choices = [(m, models_dict[m]) for m in analysed_models]

    form.traits_group.choices = [
      [json.dumps(group[0]), group[1]]
      for group in gdata["groups"].get("traits", [])]

    form.assoc_loci.choices = [
      [json.dumps(group[0]), group[1]]
      for group in gdata["groups"].get("assoc-locus", [])]

    form.gene_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("genes", [])]

    chroms  = gdata["chr_names"]
    sizes   = gdata["chr_sizes"]
    offsets = gdata["chr_offsets"]
    concat  = json.loads(config["concat"])

    fm_data = list()
    form["fm_ids"].choices = form["fm_ids"].data
    if form["fm_ids"].data is not None:
        for task_id in form["fm_ids"].data:
            fm_path = os.path.join(
              tmp_dir,
              f"{task_id}.fm")
            params_path = os.path.join(
              tmp_dir,
              f"{task_id}.params")
            with open(params_path) as params_in:
                params = json.load(params_in)
            input_params = params[1][:3]
            fm_data.append([
              input_params,
              task_id,
              "",
              pd
                .read_feather(
                  fm_path)
                .values.tolist()])

    saved_runs = list()
    saved_runs_dir = os.path.join(
      data_dir,
      genome,
      "fm")
    for run_fn in os.listdir(saved_runs_dir):
        if run_fn.endswith(".nm"):
            run_path = os.path.join(
              saved_runs_dir,
              run_fn)
            with open(os.path.join(saved_runs_dir, run_fn)) as name_in:
                run_name = name_in.read().strip()
            saved_runs.append((
              os.path.splitext(run_fn)[0],
              run_name))
    form.fm_runs.choices = saved_runs
    is_hap = dict()
    for marker_set in form.markers.choices:
        with open(os.path.join(markers_root, marker_set, "stats", "params")) as params_in:
            is_hap[marker_set] = json.load(params_in)[1][10] == "Yes"

    return render_template(
      "assoc/res.html",
      genome  = genome,
      form    = form,
      errors  = errors,
      chroms  = chroms,
      sizes   = sizes,
      offsets = offsets,
      models  = models_dict,
      concat  = concat,
      is_hap  = is_hap,
      fm_data = fm_data)

###############################################################################
### Association analysis locus fine-mapping
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Locus",
  methods=["POST"])
def assoc_locus(genome):
    form = AssocResForm(request.form)

    task_id = str(uuid.uuid4())
    task_spec = ("finemap", genome)
    task = (
      assoc_finemap,
      {"genome": genome,
      "params": form.data,
      "task_id": task_id})

    task_queue.put((task_id, task_spec, task))
    if not form.data["optimize"]:
        form["fm_ids"].data.append(task_id)
    return render_template(
      "loading.html",
      genome    = genome,
      task      = "finemap" if not form.data["optimize"] else "rf_optim",
      task_id   = task_id,
      task_spec = task_spec,
      form      = form.data)
