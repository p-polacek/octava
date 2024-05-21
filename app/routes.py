#!/usr/bin/env python

import app.path
from __main__ import (
  app,
  executor)

from app.dirs import (
  tmp_dir,
  data_dir,
  root_dir)

import os
import json
import tabix
import shutil
import hashlib
import functools
import subprocess
import pandas as pd

from flask import (
  request,
  url_for,
  redirect,
  Response,
  render_template)

from app.ajax import (
  group,
  add_phe,
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

from modules.assoc import (
  run_qc,
  filt_qc,
  run_assoc)

from modules.forms import (
  AssocResForm,
  AssocQCform,
  AlleleMinerForm)

from modules.load_input import (
  load_gdata,
  load_config)


###############################################################################
### Handling long tasks
###############################################################################

def task_stream(future):
    try:
        if future.result() == "Done.":
            return "data: Done.\n\n"
    except Exception as e:
        print(e)
        e_args = "\n".join([str(arg) for arg in e.args]).strip().split("\n")
        return f"data: {', '.join(e_args)}\n\n"

@app.route("/stream/<task_id>")
def stream(task_id):
    #future = executor.futures._futures[task_id]
    return Response(
      task_stream(future),
      content_type="text/event-stream")

def filt_done(
  genome,
  markers,
  task_id,
  future):
    if future.exception():
        markers_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers)
        if os.path.exists(markers_dir):
            shutil.rmtree(markers_dir)

    if task_id in executor.futures._futures.keys():
        _ = executor.futures.pop(task_id)

def assoc_done(
  genome,
  markers,
  model,
  traits,
  task_id,
  future):
    if task_id in executor.futures._futures.keys():
        _ = executor.futures.pop(task_id)

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
               ".gdata"))])))
    return render_template(
      "index.html",
      genomes = genomes,
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
    gconfig = load_config(genome)
    form = AlleleMinerForm(request.form)

    form.chrom.choices = gdata["chr_names"]
    form.samples.choices = gdata["samples"]

    form.sample_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("samples", [])]
    form.gene_groups.choices = [
      (json.dumps(group[0]), group[1])
       for group in gdata["groups"].get("genes", [])]

    (gts,
     groups,
     genes,
     anns,
     errors) = am_pipeline(genome, form)

    if gconfig["exp"]:
        exp_path = os.path.join(
          data_dir,
          genome,
          gconfig["exp"])

        runs_path = os.path.join(
          data_dir,
          genome,
          gconfig["runs"])
        with open(runs_path, "r") as runs_in:
            runs_info = [
              line.strip().split("\t")
              for line in runs_in]
            form.exp_col.choices = [
              ("exp-by_" + str(idx), col)
              for idx, col in enumerate(runs_info[0])]

        if gconfig["proj"]:
            proj_path = os.path.join(
              data_dir,
              genome,
              gconfig["proj"])
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
      gconfig["phe"])
    with open(phe_path, "r") as phe_in:
        phe = [
          line.strip().split("\t")
          for line in phe_in if line.strip()]
    if gconfig["enr_cols"]:
        form.enr_col.choices = gconfig["enr_cols"].split(", ")

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
### Choose filters for defining new marker set
###############################################################################

@app.route(
  "/<genome>/Association_analysis/QC",
  methods=["GET"])
def assoc_qc(genome):
    gdata = load_gdata(genome)
    gconfig = load_config(genome)
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
        hist_path = os.path.join(stats_dir, metric)
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
          "gt",
          "gt.raw.gz"))]
    form.ld_chr.choices = gdata["chr_names"]

    return render_template(
      "assoc/qc.html",
      genome      = genome,
      form        = form,
      hist_data   = hist_data,
      samp_num    = samp_num,
      var_num     = var_num,
      chr_names   = gdata["chr_names"],
      chr_sizes   = gdata["chr_sizes"],
      chr_offsets = gdata["chr_offsets"],
      vis_names   = gconfig["chr_vis"],
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

    # Prepare heterozygosity filters
    if form.data["shet"]:
        shet_exc = pd.read_csv(
          f"{stats_root}.shet.csv.gz",
          compression="gzip",
          header=None)
        shet_exc = shet_exc.set_axis(
          ["ID", "Het"],
          axis=1)
        shet_exc = shet_exc[
          shet_exc["Het"] > form.data["shet"]]
        shet_exc = shet_exc["ID"].tolist()

    if form.data["vhet"]:
        vhet_exc = pd.read_csv(
          f"{stats_root}.vhet.csv.gz",
          compression="gzip",
          header=None)
        vhet_exc = vhet_exc.set_axis(
          ["ID", "Het"],
          axis=1)
        vhet_exc = vhet_exc[
          vhet_exc["Het"] > form.data["vhet"]]
        vhet_exc = vhet_exc["ID"].tolist()

    markers = form.data["name"]
    task = "filt"
    hash_key = hashlib.sha256(
      str(
        [markers,
         task]).encode()
    ).hexdigest()
    if hash_key not in executor.futures._futures.keys():
        future = executor.submit_stored(
          hash_key,
          filt_qc,
          genome,
          form.data,
          shet_exc,
          vhet_exc)
        future.add_done_callback(
          functools.partial(
            filt_done,
            genome,
            markers,
            hash_key))

        return render_template(
          "loading.html",
          genome  = genome,
          task_id = hash_key,
          task    = task,
          warning = "")
    else:
        return render_template(
          "error.html",
          error = (
            "Marker set with the same name "
            "is being created by another user. "
            "Choose another name."))

###############################################################################
### Choose parameters for association analysis
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Options",
  methods=["GET"])
def assoc_opts(genome):
    gdata = load_gdata(genome)
    gconfig = load_config(genome)
    form = AssocResForm()

    markers_root = os.path.join(
      data_dir,
      genome,
      "assoc")
    form.markers_sgl.choices = [
      d for d in os.listdir(markers_root)
      if os.path.isdir(os.path.join(markers_root, d)) and
      os.path.exists(
        os.path.join(
          markers_root,
          d,
          "gt",
          "gt.raw.gz"))]
    form.kq_from.choices = form.markers_sgl.choices

    phe_path = os.path.join(
      data_dir,
      genome,
      gconfig["phe"])
    with open(phe_path, "r") as phe_in:
        with open(phe_path, "r") as phe_in:
            try:
                phe_header = next(phe_in).strip()
            except StopIteration:
                pass
            else:
                form.traits.choices = phe_header.split("\t")[1:]

    form.chrom.choices = [""] + gdata["chr_names"]
    form.trait_groups.choices = [
      [json.dumps(group[0]), group[1]]
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
      phe_upl     = gconfig["phe_upload"],
      rm_markers  = gconfig["rm_markers"],
      chr_names   = gdata["chr_names"],
      chr_sizes   = gdata["chr_sizes"],
      chr_offsets = gdata["chr_offsets"],
      vis_names   = gconfig["chr_vis"],
      errors      = None)

###############################################################################
### Association analysis results
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Results",
  methods=["GET", "POST"])
def assoc_res(genome):
    gdata = load_gdata(genome)
    gconfig = load_config(genome)
    errors = list()
    if request.method == "GET":
        form = AssocResForm()
    else:
        form = AssocResForm(request.form)
        res_dir = os.path.join(
          data_dir,
          genome,
          "assoc",
          form.data["markers_sgl"],
          "res")
        if not os.path.exists(res_dir):
            os.makedirs(res_dir)
        error = []
        # Skip existing runs
        for trait_group in form.data["trait_groups"]:
            for trait in trait_group.split(","):
                if trait not in form.data["traits"]:
                    form.data["traits"].append(trait)

        for trait in form.data["traits"][:]:
            trait_out = os.path.join(
              data_dir,
              genome,
              "assoc",
              form.data["markers_sgl"],
              "res",
              f"{form.data['model_sgl']}.{trait}.feather")
            if os.path.exists(trait_out):
                form.data["traits"].remove(trait)
                errors.append(trait)
            else:
                _ = open(trait_out, "w").close()

        if form.data["traits"]:
            form.data["traits"] = sorted(form.data["traits"])
            hash_key = hashlib.sha256(
              str(
                [form.data["markers_sgl"],
                 form.data["model_sgl"],
                 form.data["traits"]]).encode()
            ).hexdigest()

            future = executor.submit_stored(
              hash_key,
              run_assoc,
              genome,
              form.data)
            future.add_done_callback(
              functools.partial(
                assoc_done,
                genome,
                form.data["markers_sgl"],
                form.data["model_sgl"],
                form.data["traits"],
                hash_key))

            if errors:
                errors = ["Results alread exist for traits:"] + errors
            return render_template(
              "loading.html",
              genome  = genome,
              task_id = hash_key,
              task    = "assoc",
              errors  = errors)

        else:
            errors = ["Results already exist for all selected traits."]

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
    chr_vis = (
      [ch.strip() for ch in gconfig["chr_vis"].split(",")]
      if gconfig["chr_vis"] else [])         

    return render_template(
      "assoc/res.html",
      genome  = genome,
      form    = form,
      errors  = errors,
      chroms  = chroms,
      sizes   = sizes,
      offsets = offsets,
      chr_vis = chr_vis)

###############################################################################
### Association analysis locus fine-mapping
###############################################################################

@app.route(
  "/<genome>/Association_analysis/Locus",
  methods=["POST"])
def assoc_locus(genome):
    gconfig = load_config(genome)
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
      vis_win  = gconfig["vis_win"],
      height   = gconfig["height"],
      opt_hash = opt_hash)
