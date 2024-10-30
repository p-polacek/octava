#!/usr/bin/env python

import os
import re
import gzip
import json
import math
import uuid
import pysam
import shutil
import hashlib
import bed_reader
import subprocess
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

from flask_app import (
  app,
  get_lock,
  master_lock)

from modules.forms import AssocResForm
from modules.vcf_reader import VcfReader

from modules.load_input import (
  load_gdata,
  load_config)


from app.dirs import (
  tmp_dir,
  data_dir,
  tool_dir,
  plink2_path)


# Fetch columns for possible enrichment and chromosomes to show/hide
# for selected VCF/INFO files in new dataset setup
@app.route("/enr_cols/<genome>", methods=["POST"])
def enr_cols(genome):
    info = request.form.get("info")
    cols = pd.read_csv(
        os.path.join(
            data_dir,
            genome,
            info),
        sep="\t",
        compression="infer",
        header=[0],
        index_col=[0],
        nrows=0)
    return jsonify(cols.columns.tolist())

@app.route("/concat/<genome>", methods=["POST"])
def concat(genome):
    vcf = request.form.get("vcf")
    vcf = pysam.VariantFile(
      os.path.join(
        data_dir,
        genome,
        vcf))
    return jsonify(list(vcf.header.contigs))

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

# Per-locus LD
@app.route("/locus_ld", methods=["POST"])
def fetch_locus_ld():
    form = request.form
    genome = form["genome"]
    coord  = form["coord"]
    task_id = str(uuid.uuid4())

    gt_root = os.path.join(
      data_dir,
      genome,
      "stats",
      "gt")
    out_root = os.path.join(
      tmp_dir,
      task_id)

    # Extract genotype data at QTLs
    with open(f"{out_root}.bed1", "w") as coord_out:
        for locus in coord.split(", "):
            chrom, start = locus.split(":")
            start = start.replace(",", "")
            if "-" in start:
                start, end = start.split("-")
            else:
                end = start
            _ = coord_out.write("\t".join((chrom, start, end)) + "\n")
    pr = subprocess.run(
      [plink2_path,
       "--bfile", gt_root,
       "--extract", "bed1", f"{out_root}.bed1",
       "--thin-count", "1000",
       "--r2-unphased", "square0",
       "--out", out_root],
      check=True,
      text=True,
      stdout=subprocess.PIPE)
    output = pr.stdout
    try:
        num_loaded = int(re.search("(\d+) variants loaded", output).group(1))
        num_removed = int(re.search("\-\-extract bed1\: (\d+)", output).group(1))
        num_markers = num_loaded - num_removed
    except:
        num_markers = None

    r2_vals = pd.read_csv(
      f"{out_root}.unphased.vcor2",
      sep="\t",
      header=None,
      index_col=None).round(3).fillna("")
    r2_idx = pd.read_csv(
      f"{out_root}.unphased.vcor2.vars",
      sep="\t",
      header=None,
      index_col=None)

    return jsonify([r2_vals.values.tolist(), r2_idx[0].tolist(), num_markers])

# Per-chromosome block distribution
@app.route("/chr_bdist", methods=["POST"])
def fetch_chr_bdist():
    form = request.form
    genome  = form["genome"]
    markers = form["markers"]
    chrom   = form["chrom"]
    stats_dir = os.path.join(
      data_dir,
      genome,
      "assoc",
      markers,
      "stats")
    bdist_path = os.path.join(
      stats_dir,
      f"{chrom}.bdist.gz")
    chr_bdist = pd.read_csv(
      bdist_path,
      sep="\t",
      header=None,
      index_col=None).values.tolist()
    return jsonify(chr_bdist)

#Per-marker LD for locus fine-mapping
@app.route("/ld", methods=["POST"])
def fetch_ld():
    form = request.form
    marker_id = form["ch"] + ":" + form["pos"]
    task_id = str(uuid.uuid4())

    in_root = os.path.join(
      data_dir,
      form["genome"],
      "assoc",
      form["markers"],
      "gt",
      "gt")
    out_root = os.path.join(
      tmp_dir,
      task_id)
    subprocess.run(
      [plink2_path,
       "--bfile", in_root,
       "--r2-unphased",
       "--ld-snp", marker_id,
       "--ld-window-kb", "1000",
       "--ld-window-r2", "0",
       "--out", out_root],
      check=True)
    (pd
      .read_csv(
        f"{out_root}.vcor",
        sep="\t",
        header=0,
        index_col=None,
        usecols=[3,4,6])
      .to_csv(
        f"{out_root}.r2",
        sep="\t",
        header=True,
        index=False))
    return jsonify(task_id)

# Per-metric table data for GWAS QC filtering
@app.route("/filt", methods=["POST"])
def fetch_filt():
    form      = request.form
    genome    = form["genome"]
    metric    = form["metric"]
    start_row = int(form["start_row"])
    end_row   = int(form["end_row"])
    if "markers" in form.keys() and form["markers"]:
        markers = form["markers"]
        stats_root = os.path.join(
          data_dir,
          genome,
          "assoc",
          markers,
          "gt",
          "gt")
    else:
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
    if form["hap"] == "true":
        metrics = (
          "smiss",
          "vmiss",
          "shet",
          "vhet",
          "bsp",
          "bmc",
          "bhc",
          "bfhc",
          "baf",
          "pca",
          "kin")
    else:
        metrics = (
          "smiss",
          "vmiss",
          "shet",
          "vhet",
          "maf",
          "hwe",
          "dens",
          "pca",
          "kin")
    for metric in metrics:
        hist_path = os.path.join(stats_dir, metric)
        with open(hist_path, "r") as hist_in:
            hist_data[metric] = json.load(hist_in)

#    with open(os.path.join(stats_dir, "cnts"), "r") as cnts_in:
#        samp_num, var_num = json.load(cnts_in)

    return jsonify(hist_data)

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
        df = (
          pd
            .read_feather(assoc_path)
#            .drop(["Effect"], axis=1)
            .to_json(orient="values"))
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
            df = df.loc[:, df.columns[:3]]
            #df = df.drop(["Effect", "std Error"], axis=1)
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
    config = load_config(genome)

    if action == "vars":
        genome = request.form["genome"]
        chrom  = request.form["chrom"]
        pos    = int(request.form["pos"].replace(",", ""))
        field  = request.form["field"]

        # Get the requested field for the given variant position
        if field == "SEV":
            sev_path = os.path.join(
              data_dir,
              genome,
              f"sev.{config['vcf']}")
            with master_lock:
                lock = get_lock(sev_path)
            with lock:
                reader = VcfReader(sev_path)
                reader.set_region(chrom, pos - 1, pos)
                records = list(reader)
        else:
            vcf_path = os.path.join(
              data_dir,
              genome,
              config["vcf"])
            reader = VcfReader(vcf_path)
            reader.set_region(chrom, pos - 1, pos)
            records = list(reader)

        data = list()
        keys = list()
        for record in records:
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
                desc  = reader.metadata["INFO"]["SEV"][0]["Description"]
                desc  = desc.split("'")[1].split("'")[0]
                keys  = [label.strip() for label in desc.split("|")]
                data += [ann.split("|") for ann in csq]
            else:
                samp_path = os.path.join(
                  data_dir,
                  genome,
                  config["samp"])
                with master_lock:
                    lock = get_lock(samp_path)
                with lock:
                    samp_df = pd.read_csv(
                      samp_path,
                      sep="\t",
                      header=0,
                      index_col=None,
                      usecols=[0,1])
                    samples = dict(samp_df.to_records(index=False))
                field = samples.get(field, field)
                qual  = record.get_sample_call(field)
                keys  = list(qual.keys())
                data += [[qual[key] for key in keys]]

        return jsonify([data, keys])

    elif action == "pve":
        data = request.form
        if data["model"] != "mlm":
            return jsonify(
              [[["PVE calculation is only supported for MLM."]],
               []])
        chrom, pos = data["pos"].split(":")
        pos = pos.replace(",", "")
        pos_id = ":".join((chrom, pos))
        pos = int(pos)

        gt_root = os.path.join(
          data_dir,
          data["genome"],
          "assoc",
          data["markers"],
          "gt",
          "gt")
        bed = bed_reader.open_bed(f"{gt_root}.bed")

        phe_path = os.path.join(
          data_dir,
          data["genome"],
          config["phe"])
        with master_lock:
            lock = get_lock(phe_path)
        with lock:
            phe_df = pd.read_csv(
              phe_path,
              sep="\t",
              header=0,
              index_col=0,
              usecols=["Acc", data["trait"]]).loc[bed.iid].dropna()

        gts = bed.read(dtype="int8")[
            [i for i, n in enumerate(bed.iid) if n in phe_df.index],
            np.where(bed.sid == pos_id)[0][0]]
        gt_var = np.var(gts)
        phe_var = np.var(phe_df[data["trait"]])

        res_fn = os.path.join(
          data_dir,
          data["genome"],
          "assoc",
          data["markers"],
          "res",
          f"{data['model']}.{data['trait']}.tsv.gz")
        res = pysam.TabixFile(res_fn)
        records = list(res.fetch(chrom, pos - 1, pos))
        if len(records) != 1:
            return jsonify(
              [[["Error calculating PVE."]],
               []])
        beta = float(records[0].split("\t")[3])
        pve = (beta**2*gt_var) / phe_var

        if 1:
            af = sum(gts) / (gts.size * 2)
            se = float(records[0].split("\t")[4])
            pve2 = (2*(beta**2) * af * (1 - af)) / (2*(beta**2) * af * (1 - af) + (se**2) * 2 * gts.size * af * (1 - af))
            #from scipy.stats import linregress
            #slope, intercept, r_value, p_value, std_err = linregress(gts, phe_df[data["trait"]])
            #print(r_value**2)

        return jsonify(
          [[[round(pve * 100, 2)], [round(pve2 * 100, 2)]],
          ["PVE"]])

        
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
            return jsonify([params[1:], params[0]])

        except:
            return jsonify(
            [[["Error fetching parameters."]],
             []])

@app.route("/group", methods=["POST"])
def group():
    data = json.loads(request.get_data())
    genome = data["genome"]

    if data["set"] == "trait_sgl":
        config = load_config(genome)
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
            (phe_df
              .drop(
                [data["items"]],
                axis=1)
              .to_csv(
                phe_path,
                sep="\t"))
        return "", 204

    gdata_path = os.path.join(
      data_dir,
      genome,
      ".gdata")

    with master_lock:
        lock = get_lock(gdata_path)
    with lock:
        with open(gdata_path, "r") as gdata_in:
            gdata = json.load(gdata_in)
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
