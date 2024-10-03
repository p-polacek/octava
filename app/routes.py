#!/usr/bin/env python

import app.path
from __main__ import (
  app,
  lock,
  task_queue,
  task_status)

from app.dirs import (
  tmp_dir,
  data_dir,
  root_dir)

import os
import re
import json
import time
import uuid
import magic
import tabix
import shutil
import functools
import subprocess
import pandas as pd

from flask import (
  jsonify,
  request,
  url_for,
  redirect,
  Response,
  render_template)

from app.ajax import (
  group,
  fetch_ld,
  fetch_filt,
  rm_markers,
  contextmenu,
  fetch_mplot,
  fetch_bblplot)

from app.igv_conf import (
  igv_am,
  igv_fm)

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
from modules.setup import genome_setup

from modules.assoc import (
  mdbin,
  run_qc,
  filt_qc,
  run_assoc)

from modules.forms import (
  AddForm,
  EditForm,
  AssocQCform,
  AssocResForm,
  TraitViewForm,
  AlleleMinerForm)

from modules.load_input import (
  load_gdata,
  load_config)

###############################################################################
### Long task progress
###############################################################################

@app.route("/status/<task_id>", methods=["GET"])
def get_status(task_id):
    # Remove tasks finished more than a week ago
    for task_id in task_status.keys():
        if (
          task_id != "specs" and
          "date" in task_status[task_id].keys() and
          time.time() - task_status[task_id]["date"] > 604800):
            _ = task_status.pop(task_id)

    if task_id not in task_status.keys():
        if task_id not in [t[0] for t in list(task_queue.queue)]:
            return jsonify(""), 404
        else:
            return jsonify(""), 201
    else:
        with lock:
            progress = task_status[task_id]["progress"]
            task_status[task_id]["progress"] = ""
            code = task_status[task_id]["code"]
            if code == 200 or code == 500:
                task_status[task_id]["date"] = time.time()
            return jsonify(progress), code

###############################################################################
### Index page
###############################################################################

@app.route("/")
def index():
    # Remove old files from tmp_dir
    clean()
    genomes = list(
      enumerate(
        sorted(
          [p for p in os.listdir(data_dir)
           if os.path.exists(
             os.path.join(
               data_dir,
               p,
               "init"))])))
    links = {
      genome[1]: {
        "Allele_miner": url_for("am_opts", genome=genome[1]),
        "Trait_view": url_for("trt", genome=genome[1]),
        "Association_analysis": url_for("assoc", genome=genome[1]),
        "Edit_dataset": url_for("edit", genome=genome[1])}
      for genome in genomes}
    return render_template(
      "index.html",
      genomes = genomes,
      links   = json.dumps(links),
      genome  = "")

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
      "ann": (".gff", ".gtf", ".gff3", ".gff.gz", ".gtf.gz", ".gff3.gz"),
      "seq": (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")}
    exts.update(
      dict.fromkeys(
        ["phe", "exp", "run", "proj", "samp", "info"],
        (".txt", ".tsv", ".txt.gz", ".tsv.gz")))
    for ftype, ext in exts.items():
        form[ftype].choices = [""] + [
          f for f in os.listdir(gdir)
          if f.endswith(ext)]
    form.bam.choices = [
          f for f in os.listdir(gdir)
          if f.endswith(".bam")]

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
    icon = request.files["icon"]#
    warning = ""
    error = ""
    if icon.filename:
        icon.seek(0, os.SEEK_END)
        file_size = icon.tell()
        icon.seek(0)
        if file_size > 102400:
            error = "Submitted icon is too large."
        else:
            mime_types = ["image/jpg", "image/jpeg", "image/png"]
            ext = os.path.splitext(icon.filename)[-1]
            if (
              magic.from_buffer(icon.read(2048), mime=True) not in mime_types
              or ext not in (".png", ".jpg", ".jpeg")):
                error = "Icon is the wrong file type."
            else:
                icon.seek(0) 
                icon.save(
                  os.path.join(
                    app.config["UPLOAD_FOLDER"],
                    genome + ext))
    if error:
        return render_template(
          "error.html",
          error=error)

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
          warning   = warning)

@app.route(
  "/<genome>/Edit_dataset",
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

            safe_cols = [
              re.sub("\W", "_", col)
              for col in new_phe.columns]
            new_phe = new_phe.set_axis(
              safe_cols,
              axis=1)
            col_dup = new_phe.columns.intersection(phe_df.columns)
            if not col_dup.empty:
                col_dup = ", ".join(col_dup.tolist())
                errors.append(f"<b>Duplicate traits:</b><br>{col_dup}")
            else:
                (pd
                  .concat(
                    [phe_df, new_phe],
                    axis=1)
                  .to_csv(
                    phe_path,
                    sep="\t"))
                errors = 0

    return render_template(
      "new.html",
      genome  = genome,
      form    = form,
      edit    = True,
      allow   = allow,
      errors   = errors,
      action  = url_for("edit", genome=genome))

###############################################################################
### Trait/sample view
###############################################################################

@app.route(
  "/<genome>/Trait_view",
  methods=["GET", "POST"])
def trt(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = TraitViewForm()

    phe_path = os.path.join(
      data_dir,
      genome,
      config["phe"])
    with open(phe_path, "r") as phe_in:
        phe = [
          line.strip("\n").split("\t")
          for line in phe_in if line.strip()]
        if phe:
            form.trait_sgl.choices = phe[0][1:]
            form.traits.choices = phe[0][1:]

    form.trait_groups.choices = [
      [json.dumps(group[0]), group[1]]
      for group in gdata["groups"].get("traits", [])]
    form.samples.choices = [line[0] for line in phe[1:]]
    form.sample_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("samples", [])]

    return render_template(
      "trt/trt.html",
      genome  = genome,
      form    = form,
      phe     = phe,
      samples = gdata["samples"])

###############################################################################
### Allele miner options
###############################################################################

@app.route(
  "/<genome>/Allele_miner",
  methods=["GET", "POST"])
def am_opts(genome):
    gdata = load_gdata(genome)
    form = AlleleMinerForm()

    form.chrom.choices = gdata["chr_names"]
    form.samples.choices = gdata["samples"]

    form.sample_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("samples", [])]
    form.gene_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("genes", [])]

    for key, value in request.args.items():
        form[key].data = value

    return render_template(
      "am/opts.html",
      genome = genome,
      form   = form,
      id_fmt = gdata["id_fmt"])

###############################################################################
### Allele miner results
###############################################################################

@app.route(
  "/<genome>/Allele_miner/Results",
  methods=["POST"])
def am_res(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = AlleleMinerForm(request.form)

    form.chrom.choices = gdata["chr_names"]
    form.samples.choices = gdata["samples"]

    form.sample_groups.choices = [
      [json.dumps(group[0]), group[1]]
       for group in gdata["groups"].get("samples", [])]
    form.gene_groups.choices = [
      [json.dumps(group[0]), group[1]]
       for group in gdata["groups"].get("genes", [])]
    form.exprun_groups.choices = [
      [json.dumps(group[0]), group[1]]
       for group in gdata["groups"].get("expruns", [])]

#    form.trait_groups.choices = [
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
        # Skip header line
        if genes:
            for gene in genes[1:]:
                try:
                    for gene_entry in tb.query(
                      gene[1],
                      gene[2],
                      gene[3]):
                        if gene_entry[0] != gene[0]:
                            continue
                        exp_data.append([gene_entry[0]] + gene_entry[4:])
                except tabix.TabixError:
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
    with open(phe_path, "r") as phe_in:
        phe = [
          line.strip("\n").split("\t")
          for line in phe_in if line.strip()]
    if config["enr"]:
        form.enr_col.choices = json.loads(config["enr"])
    igv_conf = igv_am(genome)

    return render_template(
      "am/res.html",
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
    saved_sets = [
      d for d in os.listdir(markers_root)
      if os.path.isdir(os.path.join(markers_root, d)) and
      os.path.exists(
        os.path.join(
          markers_root,
          d,
          "init"))]
    form.ld_chr.choices = gdata["chr_names"]

    return render_template(
      "assoc/qc.html",
      genome      = genome,
      form        = form,
      hist_data   = hist_data,
      samp_num    = samp_num,
      var_num     = var_num,
      mdbin       = mdbin,
      chr_names   = gdata["chr_names"],
      chr_sizes   = gdata["chr_sizes"],
      chr_offsets = gdata["chr_offsets"],
      concat      = json.loads(config["concat"]),
      saved_sets  = saved_sets)

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
          warning   = "")

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
    with open(phe_path, "r") as phe_in:
        phe_header = next(phe_in)
        form.traits.choices = phe_header.strip().split("\t")[1:]

    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")
    form.markers_sgl.choices = [
      m for m in os.listdir(markers_root)
      if os.path.isdir(os.path.join(markers_root, m)) and
      os.path.exists(
        os.path.join(
          markers_root,
          m,
          "init"))]
    form.kq_from.choices = form.markers_sgl.choices

    form.chrom.choices = [""] + gdata["chr_names"]
    form.trait_groups.choices = [
      [json.dumps(group), group[1]]
      for group in gdata["groups"].get("traits", [])]

    form.plot.choices += [
      ("kin",
       "Kinship matrix"),
      ("pca",
       "Principal components")]
    form.ld_chr.choices = gdata["chr_names"]

    return render_template(
      "assoc/opts.html",
      genome      = genome,
      form        = form,
      phe_upl     = True,
      rm_markers  = True,
      mdbin       = mdbin,
      chr_names   = gdata["chr_names"],
      chr_sizes   = gdata["chr_sizes"],
      chr_offsets = gdata["chr_offsets"],
      concat      = json.loads(config["concat"]),
      errors      = None)

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
    if request.method == "GET":
        form = AssocResForm()
    else:
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
            for group in data["trait_groups"]:
                group = json.loads(group)[0]
                for trait in group:
                    if trait not in data["traits"]:
                        data["traits"].append(trait)

            existing_traits = list()
            for trait in data["traits"][:]:
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
                  task      = "assoc",
                  task_id   = task_id,
                  task_spec = task_spec,
                  errors    = errors)

            else:
                errors = ["Results already exist for all selected traits."]
        else:
            if data["traits"]:
                errors = ["To use multi-trait MLM, traits must be grouped."]
            else:
                data["trait_groups"] = [
                  json.loads(g)
                  for g in data["trait_groups"]]
                existing_groups = list()
                for group in data["trait_groups"][:]:
                    group_name = group[1]
                    # Skip existing runs
                    trait_out = os.path.join(
                      data_dir,
                      genome,
                      "assoc",
                      data["markers_sgl"],
                      "res",
                      f"{data['model_sgl']}.{group_name}.feather")
                    if os.path.exists(trait_out):
                        data["trait_groups"].remove(group)
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
                            if group[1] in data["trait_groups"]:
                                data["trait_groups"].remove(group[1])
                                queued_groups.append(group[1])
                if queued_groups:
                    errors.append("Results are already being generated for:")
                    errors.append(", ".join(queued_groups))

                if data["trait_groups"]:
                    task_id = str(uuid.uuid4())
                    task_spec = (
                      "run_assoc",
                      genome,
                      data["markers_sgl"],
                      data["model_sgl"],
                      data["trait_groups"])
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
                      errors    = errors)
                else:
                    errors = ["Results already exist for all selected groups."]

    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")

    form.markers.choices = [
      d for d in os.listdir(markers_root)
      if os.path.isdir(os.path.join(markers_root, d))]
    form.markers.data = [form.markers_sgl.data]

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
    form.models.choices = analysed_models

    form.trait_groups.choices = [
      [json.dumps(group[0]), group[1]]
      for group in gdata["groups"].get("traits", [])]

    chroms  = gdata["chr_names"]
    sizes   = gdata["chr_sizes"]
    offsets = gdata["chr_offsets"]
    concat  = json.loads(config["concat"])

    return render_template(
      "assoc/res.html",
      genome  = genome,
      form    = form,
      errors  = errors,
      chroms  = chroms,
      sizes   = sizes,
      offsets = offsets,
      concat  = concat)

###############################################################################
### Association analysis locus fine-mapping
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Locus",
  methods=["POST"])
def assoc_locus(genome):
    config  = load_config(genome)
    form    = AssocResForm(request.form)
    markers = form.data["mrk"]
    model   = form.data["mod"]
    trait   = form.data["trt"]
    pval    = form.data["pval"]
    coord   = form.data["coord"]
    gwas_fn = f"{model}.{trait}.tsv.gz"

    r2  = pd.DataFrame()
    pip = pd.DataFrame()
    for locus in coord.split(", "):
        chrom, span = locus.split(":")
        start, end  = span.replace(",", "").split("-")

        # Extract summary statistics for the region
        

        # Calculate LD matrix within the regions
        plink_path = os.path.join(
          root_dir,
          "modules",
          "tools",
          "plink2")
        gt_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "gt")
        gt_root = os.path.join(
          gt_dir,
          "gt")
        opt_hash = hashlib.sha256(
                     str(
                       [markers,
                        model,
                        coord]).encode()
                     ).hexdigest()
        out_root = os.path.join(
          tmp_dir,
          opt_hash)

        subprocess.run(
          [plink_path,
           "--bfile", gt_root,
           "--allow-extra-chr",
           "--r2-phased", "square",
           "--chr", chrom,
           "--from-bp", str(start),
           "--to-bp", str(end),
           "--out", out_root],
          check=True)

        # Run FM using SuSiE
        assoc_root = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "res")

        susie_path = os.path.join(
          root_dir,
          "modules",
          "R",
          "SuSiE.R")
        with open(
          os.path.join(
            assoc_root,
            f"{model}.{trait}.params"), "r") as params_in:
            params = json.load(params_in)
            var = params[1][4]
            num_phe = params[1][1]

        subprocess.run(
          ["Rscript", susie_path,
           "-a", os.path.join(assoc_root, gwas_fn),
           "-r", f"{out_root}.phased.vcor2",
           "-v", str(var),
           "-n", str(num_phe),
           "-c", chrom,
           "-s", str(start),
           "-e", str(end),
           "-o", assoc_root],
          check=True)

        r2_locus = pd.read_csv(
          f"{out_root}.phased.vcor2",
          sep       = "\t",
          header    = None,
          index_col = None)
        r2_locus_vars = pd.read_csv(
          f"{out_root}.phased.vcor2.vars",
          header = None)[0]
        r2_locus = r2_locus.set_axis(
          r2_locus_vars,
          axis=0)
        r2_locus = r2_locus.set_axis(
          r2_locus_vars,
          axis=1)
        r2 = pd.merge(
          r2,
          r2_locus,
          left_index  = True,
          right_index = True,
          how         = "outer")
        r2_vars = r2.index

        pip = pd.concat(
          [pip, 
           pd.read_csv(
             os.path.join(assoc_root, "pip.txt"),
             sep="\t")])

    pip.to_csv(
      f"{out_root}.pip",
      sep   = "\t",
      index = False)
    r2.to_csv(
      f"{out_root}.vcor2",
      sep    = "\t",
      index  = False,
      header = False)
    r2_vars.to_series().to_csv(
      f"{out_root}.vcor2.vars",
      index  = False,
      header = False)

    igv_conf = igv_fm(
      genome,
      markers,
      model,
      trait,
      coord,
      pval,
      opt_hash)

    return render_template(
      "assoc/locus.html",
      genome   = genome,
      markers  = markers,
      model    = model,
      igv_conf = igv_conf,
      vis_win  = config["igv_win"],
      height   = config["igv_h"],
      opt_hash = opt_hash)
