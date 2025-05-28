#!/usr/bin/env python

import io
import os
import re
import gzip
import json
import math
import uuid
import pysam
import shutil
import hashlib
import gffutils
import bed_reader
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO
import statsmodels.api as sm
from collections import Counter
from glimix_core.lmm import LMM
from scipy.stats import hypergeom
from numpy_sugar.linalg import economic_qs_linear

from flask import (
  jsonify,
  request,
  url_for,
  redirect,
  send_file,
  render_template)

import app.path

from flask_app import (
  app,
  get_lock,
  master_lock)

from modules.forms import AssocResForm

from modules.forms import models as models_dict
models_dict = dict(models_dict)

from modules.vcf_reader import VcfReader

from modules.load_input import (
  load_gdata,
  load_config)


from app.dirs import (
  tmp_dir,
  data_dir,
  tool_dir,
  plink2_path)

from modules.assoc.kq import kq_process

# Fetch columns for possible enrichment from INFO file
# and chromosomes to show/hide from VCF file in new dataset setup
@app.route("/enr_cols/<genome>", methods=["POST"])
def enr_cols(genome):
    info_file = request.form.get("info")
    cols = pd.read_csv(
        os.path.join(
            data_dir,
            genome,
            info_file),
        sep="\t",
        header=0,
        index_col=0,
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
    raw_attrs = json.loads(data.get("attrs"))
    num_subset = len(raw_attrs)

    attrs = []
    for gene in raw_attrs:
        attrs += list(
        set(
        [attr.strip() for attr
         in gene.strip().split(";") if attr]))

    enr_path = os.path.join(
      data_dir,
      genome,
      ".enr")
    with open(enr_path, "r") as enr_in:
        enr, num_genes = json.load(enr_in)

    enr = enr[col]
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
    markers = form["markers"]
    coord  = form["coord"]
    task_id = str(uuid.uuid4())

    gt_root = os.path.join(
      data_dir,
      genome,
      #"stats",
      "assoc",
      markers,
      "gt",
      "gt.full")
    out_root = os.path.join(
      tmp_dir,
      task_id)

    # Extract genotype data at QTLs
    span = 0
    with open(f"{out_root}.bed1", "w") as coord_out:
        for locus in coord.split(", "):
            chrom, start = locus.split(":")
            start = start.replace(",", "")
            if "-" in start:
                start, end = start.split("-")
                span += int(end) - int(start) + 1
            else:
                end = start
            _ = coord_out.write("\t".join((chrom, start, end)) + "\n")
    pr = subprocess.run(
      [plink2_path,
       "--bfile", gt_root,
       "--extract", "bed1", f"{out_root}.bed1",
       #"--thin-count", "1000",
       "--bp-space", str(max(1, span // 2000)),
       "--r2-unphased", "square0",
       #"--silent",
       "--out", out_root],
      check=True,
      text=True,
      stdout=subprocess.PIPE)
    output = pr.stdout
    try:
        rm, rt = map(int, re.search(
          "\-\-bp\-space\: (\d+) variants removed \((\d+) remaining\)\.",
          output).groups())
        num_markers = f"LD for {rt} of {rm + rt} markers in region"
    except:
        num_markers = "(Could not determine marker number)"

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
    markers = form["markers"]
    marker_id = form["marker_id"].replace(",", "")
    window = int(form["window"])
    r2 = float(form["r2"])
    task_id = str(uuid.uuid4())
    is_hap = False
    if "-" in marker_id and "_" in marker_id:
        is_hap = True
        in_root = os.path.join(
          data_dir,
          form["genome"],
          "assoc",
          markers,
          "gt",
          "gt.orig")
    else:
#        in_root = os.path.join(
#          data_dir,
#          form["genome"],
#          "stats",
#          "gt")
        in_root = os.path.join(
          data_dir,
          form["genome"],
          "assoc",
          markers,
          "gt",
          "gt.orig")
    out_root = os.path.join(
      tmp_dir,
      task_id)
    if is_hap:
        ch, pos = marker_id.split("_")[0].split(":")
        start, end = map(int, pos.split("-"))
    else:
        ch, pos = marker_id.split(":")
        start = end = int(pos)
    subprocess.run(
      [plink2_path,
       #"--bfile", "chr6",
       "--bfile", in_root,
       "--chr", ch,
       "--from-bp", str(start - window*1000),
       "--to-bp", str(end + window*1000),
       "--r2-unphased", "cols=id" if is_hap else "cols=pos", #"cols=chrom,pos,id,dprimeabs",
       "--ld-snp", marker_id,
       "--ld-window-kb", str(window),
       "--ld-window-r2", "0",
       #"--silent",
       "--out", out_root],
      check=True)
    ld = (pd
      .read_csv(
        f"{out_root}.vcor",
        sep="\t",
        header=0,
        index_col=None,
        usecols=[1,2]))

    if r2 > 0:
        # Remove all entries below r2
        ld.loc[ld["UNPHASED_R2"] < r2, "UNPHASED_R2"] = 0
        # Keep the entries for plotting, reduce each 0-block to first and last
        zero_rows = ld["UNPHASED_R2"] == 0
        block = (zero_rows != zero_rows.shift()).cumsum()
        mask = ~zero_rows
        zero_groups = ld[zero_rows].groupby(block)
        keep_idx = zero_groups.head(1).index.union(zero_groups.tail(1).index)
        ld = ld[mask | ld.index.isin(keep_idx)].reset_index(drop=True)

    if is_hap:
        ld.set_axis(
          ["id", "val"],
          axis=1)
        ld[["chr", "start", "end", "idx"]] = ld["id"].str.split("[\:\-\_]", expand=True)
        ld["start"] = ld["start"].astype(int)
        ld["end"] = ld["end"].astype(int)

        ch, start, end, idx = re.split("[\:\-\_]", marker_id)
        ld.loc[len(ld)] = [marker_id, 1, ch, int(start), int(end), idx]
        ld = ld.sort_values(by="start").groupby("start", as_index=False).apply(lambda x: x.loc[x["val"].idxmax()])
        ld = (ld
          [["chr", "start", "val", "end", "idx"]]
          .values.tolist())
        ld = sorted(ld, key=lambda x: x[1])
        
    else:
        ch, pos = marker_id.split(":")
        ld = (ld
          .assign(ch=ch)
          [["ch","POS_B", "UNPHASED_R2"]]
          .values.tolist())
        ld.append([ch, int(pos), 1])
        ld = sorted(ld, key=lambda x: x[1])
    return jsonify(ld)

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

# Correlation values for >=10 traits
@app.route("/correl", methods=["POST"])
def fetch_correl():
    form   = request.form
    genome = form["genome"]
    traits = json.loads(form["traits"])
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
          header=0,
          index_col=0)

    return jsonify(
      (phe_df
        [traits]
        .corr()
        .mask(np.triu(np.ones((len(traits), len(traits)), dtype=bool), k=1))
        .replace(np.nan, "")
        .values
        .tolist()))

# Filter histograms per marker set
@app.route("/hist", methods=["POST"])
def fetch_hist():
    form = request.form
    genome  = form["genome"]
    markers = form["markers"]

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
        metrics = [
          "smiss",
          "vmiss",
          "shet",
          "vhet",
          "maf",
          "hwe",
          "dens"]
        if markers:
            metrics += [
              "pca",
              "kin"]
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

            params_fn = os.path.join(
              data_dir,
              genome,
              "assoc",
              markers,
              "stats",
              "params")
            with open(params_fn) as params_in:
                is_hap = json.load(params_in)[1][10] == "Yes"
#            if is_hap:
#                continue

            df = pd.read_feather(assoc_path)
            df = df.loc[:, df.columns[:3]]
            #df = df.drop(["Effect", "std Error"], axis=1)
            df = df.set_axis([string_key, "Chr", "Pos"], axis=1)
            df = df[df[string_key] >= thresh]
            if is_hap:
                df = df.groupby(["Chr", "Pos"]).max()
            else:
                df = df.set_index(["Chr", "Pos"])
            group_data.append(df[string_key])
            #print(df[string_key])

        if not group_data:
            win_data.append(win_group)
            continue
           
        group_df = pd.concat(group_data, axis=1)

        for ch, data, in group_df.groupby(level="Chr"):
            data.index = data.index.droplevel()
            data = data.sort_index()

            spans = data.index.map(lambda x: (x // window) * window).unique()
            spans = spans.union((spans + window)).unique()

            for i in spans:
                chunk = data[data.index.slice_indexer(i, i + window)]
                if not chunk.empty:
                    chunk = chunk.dropna(axis=1, how="all")
                    if chunk.shape[1] >= min_size:
                        coord = f"{ch}:{i:,}-{i+window:,}"
                        agg_chunk = chunk.max(axis=1)
                        win_group.append(
                          [i + offsets[ch],
                           agg_chunk.max(),
                           coord,
                           chunk.columns.tolist()])
        win_data.append(win_group)
    return jsonify(win_data)

@app.route("/genes", methods=["POST"])
def fetch_genes():
    form     = request.form
    genome   = form["genome"]
    ranges   = json.loads(form["ranges"])

    config = load_config(genome)
    ann_path = os.path.join(
      data_dir,
      genome,
      config["ann"])
    ann = pysam.TabixFile(ann_path)

    gene_data = list()
    for r in ranges:
        for gene in ann.fetch(r[0], int(r[1]), int(r[2])):
            entry = gene.split("\t")
            if entry[2] == "gene":
                gene_data.append(
                  [entry[8].split("ID=")[1].split(";")[0],
                   entry[0],
                   int(entry[3]),
                   int(entry[4])])
    return jsonify(gene_data)

@app.route("/gene_pos", methods=["POST"])
def fetch_gene_pos():
    form     = request.form
    genome   = form["genome"]
    genes   = json.loads(form["genes"])

    db_path = os.path.join(
      data_dir,
      genome,
      "gff.db")
    db = gffutils.FeatureDB(db_path)

    recs = list()
    for gene in genes:
        try:
            rec = db[gene]
            if rec not in recs:
                recs.append(rec)
        except:
            continue
    gene_data = [
      [rec.id, rec.chrom, rec.start, rec.end]
      for rec in recs]

    return jsonify(gene_data)

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

        chrom, pos = data["pos"].split(":")
        pos = pos.replace(",", "")
        pos_id = ":".join((chrom, pos))
        #pos = int(pos)

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
            phe = (pd
              .read_csv(
                phe_path,
                sep="\t",
                header=0,
                index_col=0)
              .loc[bed.iid, data["trait"]]
              .dropna())

        marker_idx = [
          idx for idx, marker in enumerate(open(f"{gt_root}.bim"))
          if marker and marker.split("\t", 2)[1] == pos_id]
        if len(marker_idx) != 1:
            return jsonify(
              [[["Error calculating PVE."]],
               []])
        #marker_idx = marker_idx[0]

        samples = set(phe.index)
        sample_idxs = [
          idx for idx, sample in enumerate(open(f"{gt_root}.fam"))
          if sample and sample.split("\t", 2)[1] in samples]
        if len(sample_idxs) != len(samples):
            return jsonify(
              [[["Error calculating PVE."]],
               []])

        gts = bed.read(dtype="int8", index=(sample_idxs, marker_idx)).flatten()

        q_matr, k_matr = bkq_process(
          data["genome"],
          data["markers"],
          "SNPs",
          phe.index,
          True,
          3,
          None)

        phe_q = sm.add_constant(q_matr.loc[phe.index])
        phe_k = phe_k = k_matr.loc[phe.index, phe.index]
        QS = economic_qs_linear(phe_k)
        lmm = LMM(
          phe,
          phe_q,
          QS=QS)
        lmm.fit()
        phe_corr = phe - lmm.mean()

        #pve = sm.OLS(phe, sm.add_constant(gts)).fit().rsquared
        pve = sm.OLS(phe_corr, sm.add_constant(gts)).fit().rsquared

        if 0:
            return jsonify(
              [[["Error calculating PVE."]],
               []])

        return jsonify(
          [[[round(pve * 100, 2)]],
          ["PVE"]])

        """
        if data["model"] != "MLM":
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
        """

    elif action == "gene_seq":
        gene = request.form["gene"]
        gff_db = gffutils.FeatureDB(os.path.join(data_dir, genome, "gff.db"))
        transcripts = [
          t.id for t in gff_db.children(gene)
          if t.featuretype == "mRNA"]

        ids = list()
        seqs = list()
        seq_types = list()
        for transcript in transcripts:
            cds_path = os.path.join(
              data_dir,
              genome,
              "prot_db",
              "cds",
              f"{transcript}.fasta")
            if os.path.exists(cds_path):
                seq = next(SeqIO.parse(
                  cds_path,
                  "fasta"))
                #seq = re.sub(r"(\w{1,60})", r"\1<br>", seq)
                seqs.append(str(seq.seq))
                seq_types.append("CDS")
                ids.append(seq.id)

            prot_path = os.path.join(
              data_dir,
              genome,
              "prot_db",
              "prot",
              f"{transcript}.fasta")
            if os.path.exists(prot_path):
                seq = next(SeqIO.parse(
                  prot_path,
                  "fasta"))
                #seq = re.sub(r"(\w{1,60})", r"\1<br>", seq)
                seqs.append(str(seq.seq))
                seq_types.append("Protein")
                ids.append(seq.id)

        return jsonify(
          [[seq_types, seqs],
           ids])

    elif action == "cov":
        cov_cols = (pd
          .read_csv(
            os.path.join(data_dir, genome, "covariates", request.form["cov"]),
            sep="\t",
            header=0,
            index_col=0,
            nrows=0)
          .columns
          .tolist())
        return jsonify(
          [[[c] for c in cov_cols],
           ["Columns"]])

    elif action == "markers":
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

    # If deleting a single trait
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
              header=0,
              index_col=0)
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
        gdata["group"].setdefault(data["set"], [])

        # If making a group from a list of items
        if "name" in data.keys():
            group_entry = (
              data["items"],
              data["name"])
            # If a group with the given name exists, return error
            for group in gdata["group"][data["set"]]:
                if group[1] == data["name"]:
                    return "Group with this name alrady exists.", 400
            gdata["group"][data["set"]].append(group_entry)
        # If deleting one or multiple groups

        else:
            filt = [
              group for group in gdata["group"][data["set"]]
              if group[1] not in data["items"]]
            gdata["group"][data["set"]] = filt
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

@app.route("/save_fm", methods=["POST"])
def save_fm():
    data = json.loads(request.get_data())
    genome = data["genome"]
    task_id = data["run"]
    name = data["name"]
    fm_dir = os.path.join(data_dir, genome, "fm")

    try:
        for ext in ("fm", "params"):
            shutil.copy(
              os.path.join(tmp_dir, f"{task_id}.{ext}"),
              fm_dir)
        with open(os.path.join(fm_dir, f"{task_id}.nm"), "w") as name_out:
            _ = name_out.write(name)
        return "", 204
    except:
        return "Could not save analysis.", 500