#!/usr/bin/env python

import io
import os
import sys
import json
import math
import tabix
#import redis
#import shutil
import secrets
import gffutils
import subprocess
import numpy as np
import pandas as pd
from time import time

from flask import (Flask,
                   jsonify,
                   request,
                   session,
                   url_for,
                   redirect,
                   Response,
                   render_template,
                   send_from_directory)

from flask_compress import Compress

from modules.Forms import (GWASForm,
                           EmailForm,
                           GWASQCform,
                           HaploSmasherForm)

from modules.parseVcf import VcfReader
from modules.makeDeleteGroup import makeDeleteGroup
from modules.GWAS import run_QC, filter_QC, run_GWAS
from modules.HaploSmasher import HaploSmasherPipeline
from modules.pl_manhattan import manhattanplot, qqplot
from modules.getQueryOptions import getFile, getDbDir, getGenomeData


#allowed_extensions = ['.jpg', '.jpeg', '.png', '.gif']
##### MONGODB ANITIZATION       https://stackoverflow.com/questions/43925397/what-is-the-best-way-to-sanitize-inputs-with-flask-and-when-using-mongodb
#####                           https://www.google.com/search?client=firefox-b-d&q=flask+input+sanitization
##### FLASK OFFLOAD             https://www.google.com/search?client=firefox-b-d&q=flask+offload
##### IS THERE ANYTHING IN FLASK TO DO ABOUT HTTPS?

#redis_client = redis.StrictRedis(host="localhost", port=6379, db=0)
##### REDIS HSET HGET (hash, key, value) and JSON.DUMPS/LOADS

##### IMPLEMENT LOADING SCREEN: https://stackoverflow.com/questions/14525029/display-a-loading-message-while-a-time-consuming-function-is-executed-in-flask

##### --rm-dup remove duplicates from genotype data usink plink2


##### CREATE GFF>GZ FILES
##### remove blank lines, zip and tabix
##### grep "\S" ann.gff3 | bgzip -c ann.gff3.gz
##### tabix ann.gff3.gz -f


# Set up Flask app
app = Flask(__name__)
#app.config["UPLOAD_FOLDER"] = "static/icons"
app.secret_key = "275b250b2816223102a047d297676ce0"
app.config['MAX_CONTENT_LENGTH'] = 5 * 1024 * 1024
Compress(app)

# Define working directories
root_dir = os.getcwd()
data_dir = os.path.join(root_dir, "data")


@app.route("/favicon.ico")
def favicon():
    return send_from_directory("static",
                               "favicon.ico",
                               mimetype="image/vnd.microsoft.icon")


@app.route("/icon/<filename>")
def serve_icon(filename):
    icon_dir = os.path.join("static",
                            "icons")
    if not os.path.exists(os.path.join(icon_dir, f"{filename}.png")):
        filename = "Icon_default.png"
    else:
        filename += ".png"

    return send_from_directory(icon_dir,
                               filename)


@app.route("/image/<path:filename>")
def serve_image(filename):
    image_dir = os.getcwd()
    image_dir = ""
    return send_from_directory(image_dir,
                               filename)


@app.route("/data/<genome>/<fn>")
def serve_file(genome, fn):
    ### Correct this later
#    if fn.startswith("SL"):
#        return send_from_directory(os.path.join(data_dir, genome, "GWAS", "snp", "imp_filt"), "SL4.0ch11.53252504.vcor")
#    fn_dir = os.path.join(data_dir,
#                          genome,
#                          "GWAS",
#                          "snp",
#                          "assoc")
#    else:
    fn_dir = os.path.join(data_dir,
                          genome)
    if fn.endswith(".gz"):
        return send_from_directory(fn_dir,
                                   fn,
                                   mimetype="application/gzip")
    return send_from_directory(fn_dir,
                               fn)


@app.route("/markerdata/<genome>/<markers>/<fn>")
def serve_markerdata(genome, markers, fn):
    ### Correct this later
#    if fn.startswith("GAPIT"):
#    model = fn.strip("GAPIT.processed.").split(".", 1)[0]
    fn_dir = os.path.join(data_dir,
                          genome,
                          "GWAS",
                          markers,
                          "assoc")
#    else:
#        fn_dir = os.path.join(data_dir,
#                              genome)
    if fn.endswith(".gz"):
        return send_from_directory(fn_dir,
                                   fn,
                                   mimetype="application/gzip")
    return send_from_directory(fn_dir,
                               fn)


@app.route("/ld_file/<genome>/<markers>")
def serve_ld_file(genome, markers):
    return send_from_directory("tmp", "snp.ld.single.vcor2")


@app.route("/pip_file/<genome>/<markers>")
def serve_pip_file(genome, markers):
    return send_from_directory("tmp", "pip_combined.txt")


@app.route("/ld/<genome>/<markers>/<ch>/<pos>")
def get_ld(genome, markers, ch, pos):
    """
    plink_path = os.path.join(root_dir,
                              "modules",
                              "tools",
                              "plink2")
    gt_dir = os.path.join(data_dir,
                          genome,
                          "GWAS",
                          markers,
                          "imp_filt")
    gt_root = os.path.join(gt_dir,
                           "imp.filt")
    out_root = os.path.join(gt_dir,
                            "snp.ld")
    out_file = f"{out_root}.vcor"

    subprocess.run([plink_path,
                    "--bfile", gt_root,
                    "--allow-extra-chr",
                    "--r2-phased",
                    "--ld-snp", ":".join((ch, pos)),
                    "--ld-window-r2", "0",
                    "--out", out_root],
                   check=True)

    pd.read_csv(out_file, sep="\t").iloc[:, [3, 4, 6]].set_axis(["Chr", "Pos", "R2"], axis=1).to_csv(out_file, sep="\t", index=False)
    """

#    gt_dir = os.path.join(data_dir,
#                          genome,
#                          "GWAS",
#                          markers,
#                          "imp_filt")
    gt_dir = "tmp"
    in_file = os.path.join(gt_dir,
                           #"snp.ld.phased.vcor2")
                           "r2_combined.txt")
    in_vars = os.path.join(gt_dir,
                           #"snp.ld.phased.vcor2.vars")
                           "r2_combined.vars")
    out_file = os.path.join(gt_dir,
                            "snp.ld.single.vcor2")

    variants = pd.read_csv(in_vars, sep="\t", header=None)[0]
    var_index = variants[variants == f"{ch}:{pos}"].index.values[0]
    r2 = pd.read_csv(in_file, sep="\t", header=None)[[var_index]].set_axis(["R2"], axis=1)
    r2[["Chr", "Pos"]] = variants.str.split(":", expand=True)
    r2[["Chr", "Pos", "R2"]].to_csv(out_file, sep="\t", index=False)

    return jsonify(success=True)


@app.route("/fetch_ajax", methods=["POST"])
def fetch_ajax():
    genome = request.form["genome"]
    if "plot_choice" in request.form:
        metric = request.form["plot_choice"]
        start_row = int(request.form["start_row"])
        end_row = int(request.form["end_row"])

        stats_root = os.path.join(data_dir,
                                  genome,
                                  "stats",
                                  "plink_out")

        with open(f"{stats_root}.{metric}.csv", "r") as metr_in:
            for i in range(1, start_row):
                next(metr_in)
            lines = [next(metr_in).strip().split(",")
                     for l in range(start_row, end_row + 1)]

    elif "metric" in request.form:
        model = request.form["model"]
        trait = request.form["trait"]
        metric = request.form["metric"]

        start_row = int(request.form["start_row"])
        end_row = int(request.form["end_row"])

        gwas_dir = os.path.join(data_dir,
                                genome,
                                "GWAS",
                                "snp",
                                "assoc")

        gwas_fn = os.path.join(gwas_dir, f"{model}.{trait}.tsv")

        with open(gwas_fn, "r") as gwas_in:
            next(gwas_in)
            for i in range(1, start_row):
                next(gwas_in)
            lines = [next(gwas_in).strip().split("\t")
                     for l in range(start_row, end_row + 1)]

    elif "trait" in request.form:
        trait = request.form["trait"]
        phe_fn = getFile(genome, "phe")

        phe_df = pd.read_csv(phe_fn, sep="\t")
        phe_df[trait] = phe_df[trait].round(5)
        lines = phe_df[[phe_df.columns[0], trait]].values.tolist()

    return jsonify(lines)


@app.route("/fetch_mplot", methods=["POST"])
def fetch_mplot():
    genome   = request.form["genome"]
    markers  = request.form["markers"]
    model    = request.form["model"]
    trait    = request.form["trait"]

    assoc_fn = f"{model}.{trait}.feather"
    assoc_path = os.path.join(data_dir,
                              genome,
                              "GWAS",
                              markers,
                              "assoc",
                              assoc_fn)

    if os.path.exists(assoc_path):
        df = pd.read_feather(assoc_path)
        return jsonify(df.drop(["Effect"], axis=1).to_json(orient="values"))

    return jsonify(success=True)


@app.route("/fetch_bubbleplot", methods=["POST"])
def fetch_bubbleplot():
    #t = time()
    form = request.form
    genome   = form["genome"]
    window   = int(form["window"])
    thresh   = float(form["thresh"])
    keys     = json.loads(form["keys"])
    min_size = math.ceil(float(form["min_size"])) if "min_size" in form.keys() else 0
    plot_avg = int(form["plot_avg"])

    genome_data = getGenomeData(genome)
    offsets = {k: v["offset"] for k, v in genome_data["Chromosome_data"].items()}

    windowed_data = list()
    for group in keys:
        group_data = list()
        windowed_group = list()

        for key in group:
            string_key = ", ".join(key)
            assoc_fn = f"{key[1]}.{key[2]}.feather"
            assoc_path = os.path.join(data_dir,
                                      genome,
                                      "GWAS",
                                      key[0],
                                      "assoc",
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
            windowed_data.append(windowed_group)
            continue
            
        group_df = pd.concat(group_data, axis=1)#.max(axis=1)

        for chrom, data, in group_df.groupby(level=["Chr"]):
            data.index = data.index.droplevel()
            data = data.sort_index()

            spans = data.index.map(lambda x: (x // window) * window).unique()
            spans = spans.union((spans + window)).unique()

            #for i in range(data.index[0] - data.index[0] % window, data.index[-1] + window, window):
            for i in spans:
                chunk = data[data.index.slice_indexer(i - window, i + window)]
                if not chunk.empty:
                    chunk = chunk.dropna(axis=1, how="all")
                    if chunk.shape[1] >= min_size:
                        if plot_avg:
                            windowed_group.append([int(chunk.idxmax().mean()) + offsets[chrom],
                                                   chunk.max().mean(),
                                                   chunk.shape[1],
                                                   f"{chrom}:{i-window:,}-{i+window:,}",
                                                   chunk.columns.tolist()])
                        else:
                            agg_chunk = chunk.max(axis=1)
                            windowed_group.append([int(agg_chunk.idxmax()) + offsets[chrom],
                                                   agg_chunk.max(),
                                                   chunk.shape[1],
                                                   f"{chrom}:{i-window:,}-{i+window:,}",
                                                   chunk.columns.tolist()])

        windowed_data.append(windowed_group)

    #print(time() - t)
    return jsonify(windowed_data)


@app.route("/contextmenu_data", methods=["POST"])
def contextmenu_data():
    genome = request.form["genome"]
    action = request.form["action"]

    if action == "markers":
        params_json = os.path.join(data_dir,
                                   genome,
                                   "GWAS",
                                   request.form["markers"],
                                   "imp_filt",
                                   "imp.filt.params.json")
        try:
            with open(params_json, "r") as params_in:
                params = json.load(params_in)
            return jsonify([params[1:], params[0]])

        except:
            return jsonify(["Error fetching parameters.", None])

    elif action == "gwas_run":
        fn = f"{request.form['model']}.{request.form['trait']}.params.json"
        params_json = os.path.join(data_dir,
                                   genome,
                                   "GWAS",
                                   request.form["markers"],
                                   "assoc",
                                   fn)
        try:
            with open(params_json, "r") as params_in:
                params = json.load(params_in)
            return jsonify([params[1:], params[0]])

        except:
            return jsonify(["Error fetching parameters.", None])

    elif action == "vars":
        genome = request.form["genome"]
        chrom  = request.form["chrom"]
        pos    = int(request.form["pos"])
        field  = request.form["field"]

        # Get the requested field for the given variant position
        vcf_file = getFile(genome, "vcf")
        vcf_reader = VcfReader(vcf_file)
        vcf_reader.set_region(chrom, pos - 1, pos)
        for record in vcf_reader:
            if record.POS != pos:
                continue
            if field == "REF":
                data = record.REF
                keys = None
            elif field == "CSQ":
                csq  = record.CSQ
                keys = list(csq[0].keys())
                data = [[ann[key] for key in keys] for ann in csq]
            elif field == "PROV":
                csq  = record.parsed_info_fields()["PROV"]
                desc = vcf_reader.metadata["INFO"]["PROV"][0]["Description"]
                desc = desc.split("'")[1].split("'")[0]
                keys = [label.strip() for label in desc.split("|")]
                data = [ann.split("|") for ann in csq]
            else:
                qual = record.get_sample_call(field)
                keys = list(qual.keys())
                data = [[qual[key] for key in keys]]

        return jsonify([data, keys])


@app.route("/save_session/<user_id>/<genome>/<tool>")
def save_session(user_id, genome, tool):
    session[user_id] = {"genome": genome}
    return redirect(url_for(tool,
                            user_id=user_id))


@app.route("/addPhenotype/<user_id>/<genome>", methods=["POST"])
def add_phenotype(user_id, genome):
    """Add a new phenotype file"""
    phe_file = getFile(genome, "phe")
    phe = pd.read_csv(phe_file, header=[0], index_col=[0], sep="\t")

    genome_data = getGenomeData(genome)
    samples = [sample[0] for sample in genome_data["Samples"]]
    errors = list()

    for file in request.files.getlist("new_phe"):
        fn = file.filename
        new_phe = pd.read_csv(file, header=[0], index_col=[0],
                              sep="," if fn.endswith(".csv") else "\t")

#        i_diff = new_phe.index.difference(phe.index)
#        if not i_diff.empty:
#            errors.append(f"<b>Ungenotyped samples in {fn}:</b><br>{', '.join(i_diff.tolist())}")

        c_same = new_phe.columns.intersection(phe.columns)
        if not c_same.empty:
            errors.append(f"<b>Duplicate traits in {fn}:</b><br>{', '.join(c_same.tolist())}")

        if not errors:
            phe = pd.concat([phe, new_phe], axis=1)

    if not errors:
        phe.to_csv(phe_file, sep="\t")
        errors = 1

    # Return back to GWAS options page
    form = GWASForm(request.form)
    markers_root = os.path.join(data_dir,
                                genome,
                                "GWAS")
    form.markers.choices = [d for d in os.listdir(markers_root)
                            if os.path.isdir(os.path.join(markers_root, d))]
    form.traits.choices = phe.columns
    form.chrom.choices = [""] + genome_data["Chromosomes"]
    if "traits" in genome_data["Groups"].keys():
        form.trait_groups.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]["traits"]]

    return render_template("GWAS/GWAS_options.html",
                           user_id=user_id,
                           genome=genome,
                           form=form,
                           phe=phe,
                           errors=errors)

    # file.save(fn)


@app.route("/")
def index():
    user_id = secrets.token_urlsafe(16)
    genomes = list(enumerate(sorted([p for p in os.listdir(data_dir)
                                     if os.path.isdir(os.path.join(data_dir, p))])))

    return render_template("index.html",
                           user_id=user_id,
                           genomes=genomes,
                           genome="")

"""
@app.route("/login", methods=["GET", "POST"])
def login():
    form = EmailForm()

#    if form.validate_on_submit():
    if request.method == "POST":
        email = form.email.data

    return render_template("login.html",
                           form=form)
"""

@app.route("/<user_id>/Allele_browser", methods=["GET", "POST"])
def Allele_browser(user_id):
    form = HaploSmasherForm()
    genome = session.get(user_id)["genome"]

    try:
        genome_data = getGenomeData(genome)
    except OSError as error:
        return render_template("index.html",
                               errors = error)

    form.chrom.choices = genome_data["Chromosomes"]
    if "samples" in genome_data["Groups"].keys():
        form.sample_groups.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]["samples"]]
    if "genes" in genome_data["Groups"].keys():
        form.gene_groups.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]["genes"]]
    form.samples.choices = genome_data["Samples"]

    for key, value in request.args.items():
        form[key].data = value

    return render_template("Allele_browser/allele_browser.html",
                           user_id        = user_id,
                           genome         = genome,
                           form           = form,
                           sample_id      = genome_data["Sample_ID"])


@app.route("/makeDeleteGroup", methods=["POST"])
def Group():
    """Create or delete groups of samples"""
    data = json.loads(request.get_data())
    genome = session.get(data["user_id"])["genome"]
    result = makeDeleteGroup(genome, data)
    if result:
        return result, 400
    return "", 204


@app.route("/<user_id>/Allele_viewer", methods=["POST"])
def Allele_viewer(user_id):
    """Results page showing output and further options"""
    #if form.validate_on_submit():
    genome = session.get(user_id)["genome"]
    form = HaploSmasherForm(request.form)

    try:
        genome_data = getGenomeData(genome)
    except OSError as error:
        return render_template("index.html",
                               errors = error)

    form.chrom.choices = genome_data["Chromosomes"]
    if "samples" in genome_data["Groups"].keys():
        form.sample_groups.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]["samples"]]
    if "genes" in genome_data["Groups"].keys():
        form.gene_groups.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]["genes"]]
    if "expruns" in genome_data["Groups"].keys():
        form.exprun_groups.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]["expruns"]]
    form.samples.choices = genome_data["Samples"]

    (genotypes,
     groups,
     genes,
     annotations,
     errors) = HaploSmasherPipeline(genome, form)

    # Convert genes and annotations from dict to list
    genes = list(genes.values())
    genes = [{k: v for k, v in gene.items() if k != "Region"} for gene in genes]


    if genome == "Pepper v1.6 G2PSOL":
        genome_dir = os.path.join(data_dir, genome)
        exp_file = "all_quant.combined.quantiles.txt.gz"
        run_file = "run_info.txt"
        project_file = "bp_info.txt"

        with open(os.path.join(genome_dir, run_file), "r") as in_file:
            run_info = in_file.read()

        with open(os.path.join(genome_dir, project_file), "r") as in_file:
            project_info = in_file.read()

        exp_data = ""
        tb = tabix.open(os.path.join(genome_dir, exp_file))
        for gene in genes:
            for exp_entry in tb.query(gene["Chromosome"],
                                      gene["Start"],
                                      gene["End"]):
                if exp_entry[0] != gene["ID"]:
                    continue
                exp_data += "\t".join([exp_entry[0]] + exp_entry[4:]) + "\n"
    else:
        exp_data = ""
        run_info = ""
        project_info = ""

    try:
        with open("test_colors.json", "r") as in_colors:
            colors = json.load(in_colors)
    except:
        colors = ""

    seq_file = getFile(genome, "seq").split("/")[-1]
    seq_index = seq_file + ".fai"
    gff_file = getFile(genome, "gff").split("/")[-1]
    vcf_file = getFile(genome, "vcf").split("/")[-1]
    vcf_index = vcf_file + ".tbi"
    gff_index = gff_file + ".tbi"

    igv_config = {
        "reference": {
            "id": genome,
            "fastaURL":  url_for("serve_file", genome=genome, fn=seq_file),
            "indexURL":  url_for("serve_file", genome=genome, fn=seq_index),
            "tracks": [
                {
                    "name": "Genes",
                    "type": "annotation",
                    "format": "gff3",
                    "url": url_for("serve_file", genome=genome, fn=gff_file),
                    "indexURL": url_for("serve_file", genome=genome, fn=gff_index),
                    "visibilityWindow": 10000000,
                    "height": 150,
                },
                {
                    "name": "Variants",
                    "type": "variant",
                    "format": "vcf",
                    "url": url_for("serve_file", genome=genome, fn=vcf_file),
                    "indexURL": url_for("serve_file", genome=genome, fn=vcf_index),
                    "height": 150,
                    "visibilityWindow": 10000000,
                },
            ],
       },
    }

    phe_fn = getFile(genome, "phe")
    with open(phe_fn, "r") as phe_in:
        phe = [line.strip().split("\t") for line in phe_in]

    return render_template("Allele_browser/results.html",
                           user_id        = user_id,
                           genome         = genome,
                           genotypes      = genotypes,
                           groups         = groups,
                           genes          = genes,
                           annotations    = annotations,
                           errors         = errors,
                           form           = form,
                           sample_id      = genome_data["Sample_ID"],

                           phe            = phe,

                           run_info       = run_info,
                           project_info   = project_info,
                           exp_data       = exp_data,
                           colors         = colors,

                           igv_config     = igv_config,
                           )


@app.route("/<user_id>/Genome_browser", defaults={"gene_id": None, "position": None})
@app.route("/<user_id>/Genome_browser/<gene_id>/<position>", methods=["GET", "POST"])
def Genome_browser(user_id, gene_id, position):
    return

@app.route("/<user_id>/locus/<locus>")
def update_locus(user_id, locus):
    genome = session.get(user_id)["genome"]
    try:
        genome_data = getGenomeData(genome)
    except OSError as error:
        return render_template("index.html",
                               errors = error)

    form = HaploSmasherForm()
    form.chrom.choices = genome_data["Chromosomes"]
    form.groups.choices  = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]]
    form.samples.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Samples"]]

    form.fmt.data = "region"
    chrom, span = locus.split(":")
    start, end = span.split("-")
    form.chrom.data = chrom
    form.start.data, form.end.data = math.floor(float(start)), math.ceil(float(end))

    return render_template("Allele_browser/allele_browser.html",
                           user_id        = user_id,
                           genome         = genome,
                           form           = form,
                           genome_data    = genome_data)


@app.route("/<user_id>/Association_analysis", methods=["GET"])
def Association_analysis(user_id):
    genome = session.get(user_id)["genome"]

    # Define GWAS input files
    GWAS_root = os.path.join(data_dir,
                             genome,
                             "GWAS")

    if not os.listdir(GWAS_root):

        return redirect(url_for("GWAS_QC",
                                user_id=user_id))

    else:
        return redirect(url_for("GWAS_options",
                                user_id=user_id))


@app.route("/<user_id>/Association_analysis/QC", methods=["GET"])
def GWAS_QC(user_id):
    genome = session.get(user_id)["genome"]

    form = GWASQCform()

    stats_root = os.path.join(data_dir,
                              genome,
                              "stats",
                              "plink_out")

    if any((not os.path.exists(f"{stats_root}.{metric}.hist.json")
           for metric in ("smiss",
                          "vmiss",
                          "shet",
                          "vhet",
                          "maf",
                          "hwe"))):
        run_QC(genome)

    hist_data = dict()
    for metric in ("smiss", "vmiss", "shet", "vhet", "maf", "hwe"):
        with open(f"{stats_root}.{metric}.hist.json", "r") as hist_in:
            hist_data[metric] = json.load(hist_in)

    return render_template("GWAS/GWAS_QC.html",
                           user_id=user_id,
                           genome=genome,
                           form=form,
                           hist_data=hist_data)
                           #qc_data=qc_data,
                           #s_ids=s_ids,
                           #v_ids=v_ids)


@app.route("/<user_id>/Association_analysis/Filter", methods=["POST"])
def GWAS_filter(user_id):
    genome = session.get(user_id)["genome"]
    form = GWASQCform(request.form)

    stats_root = os.path.join(data_dir,
                              genome,
                              "stats",
                              "plink_out")

    # Prepare heterozygosity filters
    if form.data["shet"]:
        shet_exclude = pd.read_csv(f"{stats_root}.shet.csv", header=None)
        shet_exclude = shet_exclude.set_axis(["ID", "Het"], axis=1)
        shet_exclude = shet_exclude[shet_exclude["Het"] > form.data["shet"]]
        shet_exclude = shet_exclude["ID"].tolist()

    if form.data["vhet"]:
        vhet_exclude = pd.read_csv(f"{stats_root}.vhet.csv", header=None)
        vhet_exclude = vhet_exclude.set_axis(["ID", "Het"], axis=1)
        vhet_exclude = vhet_exclude[vhet_exclude["Het"] > form.data["vhet"]]
        vhet_exclude = vhet_exclude["ID"].tolist()

    filter_QC(genome,
              form.data,
              shet_exclude,
              vhet_exclude)

    return redirect(url_for("GWAS_options", user_id=user_id))


@app.route("/<user_id>/Association_analysis/Options", methods=["GET"])
def GWAS_options(user_id):
    genome = session.get(user_id)["genome"]

    form = GWASForm()
    markers_root = os.path.join(data_dir,
                                genome,
                                "GWAS")
    form.markers.choices = [d for d in os.listdir(markers_root)
                            if os.path.isdir(os.path.join(markers_root, d))]

    phe_fn = getFile(genome, "phe")
    with open(phe_fn, "r") as phe_in:
        phe = [line.split("\t") for line in phe_in.read().strip().split("\n")]
    form.traits.choices = phe[0][1:]

    try:
        genome_data = getGenomeData(genome)
    except OSError as error:
        return render_template("index.html",
                               errors = error)
    form.chrom.choices = [""] + genome_data["Chromosomes"]
    if "traits" in genome_data["Groups"].keys():
        form.trait_groups.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]["traits"]]

    return render_template("GWAS/GWAS_options.html",
                           user_id=user_id,
                           genome=genome,
                           form=form,
                           phe=phe,
                           errors=None)


@app.route("/<user_id>/Association_analysis/Results", methods=["POST"])
def GWAS_results(user_id):
    genome = session.get(user_id)["genome"]
    form = GWASForm(request.form)
    run_GWAS(genome, form.data)

    markers_root = os.path.join(data_dir,
                                genome,
                                "GWAS")
    form.markers.choices = [d for d in os.listdir(markers_root)
                            if os.path.isdir(os.path.join(markers_root, d))]

    analysed_traits = list()
    analysed_models = list()
    for marker_set in form.markers.choices:
        assoc_root = os.path.join(data_dir,
                                  genome,
                                  "GWAS",
                                  marker_set,
                                  "assoc")
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

    try:
        genome_data = getGenomeData(genome)
    except OSError as error:
        return render_template("index.html",
                               errors = error)

    if "traits" in genome_data["Groups"].keys():
        form.trait_groups.choices = [[json.dumps(g[0]), g[1]]  for g in genome_data["Groups"]["traits"]]

    chromosomes = genome_data["Chromosomes"]
    chromosome_data = genome_data["Chromosome_data"]

    return render_template("GWAS/GWAS_results.html",
                           user_id=user_id,
                           genome=genome,
                           form=form,
                           chromosomes=chromosomes,
                           chromosome_data=chromosome_data)


@app.route("/<user_id>/Association_analysis/Locus", methods=["POST"])
def GWAS_locus(user_id):
    genome = session.get(user_id)["genome"]

    form = GWASForm(request.form)
    markers = form.data["mrk"]
    model   = form.data["mod"]
    trait   = form.data["trt"]
    pval    = form.data["pval"]
    coord   = form.data["coord"]

    gwas_file = f"{model}.{trait}.tsv.gz"
    gwas_index = gwas_file + ".tbi"

    r2  = pd.DataFrame()
    pip = pd.DataFrame()

    for locus in coord.split(", "):
        chrom, span = locus.split(":")
        start, end  = span.replace(",", "").split("-")

        # Calculate LD matrix within the regions
        plink_path = os.path.join(root_dir,
                                  "modules",
                                  "tools",
                                  "plink2")
        gt_dir = os.path.join(data_dir,
                              genome,
                              "GWAS",
                              markers,
                              "imp_filt")
        gt_root = os.path.join(gt_dir,
                               "imp.filt")
        out_root = os.path.join(gt_dir,
                                "snp.ld")

        subprocess.run([plink_path,
                        "--bfile", gt_root,
                        "--allow-extra-chr",
                        "--r2-phased", "square",
                        "--chr", chrom,
                        "--from-bp", str(start),
                        "--to-bp", str(end),
                        "--out", out_root],
                       check=True)

        # Run FM using SuSiE
        phe_file = os.path.join(root_dir, getFile(genome, "phe"))

        assoc_root = os.path.join(data_dir,
                                  genome,
                                  "GWAS",
                                  markers,
                                  "assoc")

        R_dir = os.path.join(root_dir,
                             "modules",
                             "R")
        susie_path = os.path.join(R_dir,
                                  "SuSiE.R")

        command = ["Rscript", susie_path,
                   "-g", os.path.join(assoc_root, gwas_file),
                   "-r", f"{out_root}.phased.vcor2",
                   "-p", phe_file,
                   "-t", trait,
                   "-c", chrom,
                   "-s", str(start),
                   "-e", str(end),
                   "-o", assoc_root]
        subprocess.run(command, check=True)

        r2_locus = pd.read_csv(f"{out_root}.phased.vcor2", sep="\t", header=None, index_col=None)
        r2_locus_vars  = pd.read_csv(f"{out_root}.phased.vcor2.vars", header=None)[0]
        r2_locus = r2_locus.set_axis(r2_locus_vars, axis=0)
        r2_locus = r2_locus.set_axis(r2_locus_vars, axis=1)
        r2 = pd.merge(r2, r2_locus, left_index=True, right_index=True, how="outer")
        r2_vars = r2.index

        pip = pd.concat([pip, pd.read_csv(os.path.join(assoc_root, "pip.txt"), sep="\t")])

    pip.to_csv("tmp/pip_combined.txt", sep="\t", index=False)
    r2.to_csv("tmp/r2_combined.txt", sep="\t", index=False, header=False)
    r2_vars.to_series().to_csv("tmp/r2_combined.vars", index=False, header=False)
    

    seq_file = getFile(genome, "seq").split("/")[-1]
    seq_index = seq_file + ".fai"
    gff_file = getFile(genome, "gff").split("/")[-1]
    gff_index = gff_file + ".tbi"

    igv_config = {
        "reference": {
            "id": genome,
            "fastaURL":  url_for("serve_file", genome=genome, fn=seq_file),
            "indexURL":  url_for("serve_file", genome=genome, fn=seq_index),
            "tracks": [
                {
                    "name": "Genes",
                    "type": "annotation",
                    "format": "gff3",
                    "url": url_for("serve_file", genome=genome, fn=gff_file),
                    "indexURL": url_for("serve_file", genome=genome, fn=gff_index),
                    "visibilityWindow": 1000000,
                    "order": 0,
                },
                {
                    "name": "GWAS",
                    "type": "gwas",
                    "format": "gwas",
                    "url": url_for("serve_markerdata", genome=genome, markers=markers, fn=gwas_file),
                    "indexURL": url_for("serve_markerdata", genome=genome, markers=markers, fn=gwas_index),
                    "visibilityWindow": 1000000,
                    "columns": {
                        "chromosome": 1,
                        "position": 2,
                        "value": 3
                    },
                    "min": -0.5,
                    "max": pval,
                    "height": 200,
                    "dotSize": 3,
                    "cursor": "pointer",
                    "color": "#3a8f1d",
                },
                {
                    "name": "PIP",
                    "type": "gwas",
                    "format": "gwas",
                    "posteriorProbability": True,
                    "url": url_for("serve_pip_file", genome=genome, markers=markers),
                    "visibilityWindow": 1000000,
                    "columns": {
                        "chromosome": 1,
                        "position": 2,
                        "value": 3
                    },
                    "min": -0.1,
                    "max": 1.1,
                    "height": 200,
                    "dotSize": 3,
                    "color": "#3a8f1d",
                },
#                {
#                    "name": "LD",
#                    "type": "gwas",
#                    "format": "gwas",
#                    "posteriorProbability": True,
#                    "url": url_for("serve_ld_file"),
#                    "visibilityWindow": 1000000,
#                    "columns": {
#                        "chromosome": 1,
#                        "position": 2,
#                        "value": 3
#                    },
#                    "min": -0.1,
#                    "max": 1.1,
#                    "height": 200,
#                    "dotSize": 3,
#                    "color": "#3a8f1d",
#                },
#                {
#                    "name": "Variants",
#                    "type": "variant",
#                    "format": "vcf",
#                    "url": url_for("serve_file", genome=genome, fn=vcf_file),
#                    "indexURL": url_for("serve_file", genome=genome, fn=vcf_index),
#                    "url": url_for("in_memory_vcf"),
#                    "height": 50,
#                    "visibilityWindow": 10000000,
#                },
            ],
       },
       "locus": coord
    }

    return render_template("GWAS/GWAS_Locus.html",
                           user_id    = user_id,
                           genome     = genome,
                           markers    = markers,
                           igv_config = igv_config)


    """
    lz_data = {
          "url": url_for("serve_file", genome=genome, fn="GAPIT.processed.BLINK.Inflorescence.tsv.gz"),
          "type": "tabix",
        }

    lz_config = {
        "state": {
            "chr": "1",
            "start": 1e6,
            "end": 2e6,
            "width": 800,
        },
        "panels": [
            {
                "id": "association",
                "data_layers": [
                    {
                        "id": "association",
                        "type": "scatter",
                        "fields": ["position", "pvalue"],
                        "data_layer": {
                            "url": url_for("serve_file", genome=genome, fn="GAPIT.processed.BLINK.Inflorescence.tsv.gz"),
                            "type":"tabix",
                        },
                    },
                ],
            },
        ],
    }
            {
                "id": "annotation",
                "data_layers": [
                    {
                        "id": "genes",
                        "type": "annotation",
                        "glyph": "triangle",
                        "label": "Genes",
                        "fill": "#3182bd",
                        "border": "#08519c",
                        "height": 8,
                        "data_layer": {
                            "url": url_for("serve_file", genome=genome, fn=getFile(genome, "gff")),
                            "type": "gff",
                        },
                    },
                ],
            },
    """


@app.route("/<user_id>/Phylogeny_analysis", methods=["GET", "POST"])
def Phylogeny_analysis(user_id):
    return

if __name__ == "__main__":
    app.run(debug=True,
            #host="0.0.0.0",
            port=8888)


"""

    igv_config = {
        "reference": {
            "id": "Capsicum_annuum",
            "fastaURL":  url_for('serve_file', genome=genome, filename="noChr.Annuum.v1.6.Total.seq.fa"),
            "indexURL":  url_for('serve_file', genome=genome, filename="/noChr.Annuum.v1.6.Total.seq.fa.fai"),
            "tracks": [
                {
                    "name": "Refseq Genes",
                    "type": "annotation",
                    "format": "gff3",
                    "url": url_for('serve_file', genome=genome, filename="nameAdded.noChr.Annuum.v.2.0.chromosome.gff3"),
                },
                {
                    "name": "VCF Track",
                    "type": "variant",
                    "format": "vcf",
                    #"url": url_for('serve_file', genome=genome, filename="var_only.ann.G2P-SOL_Pepper_CC_Geno_data_2023.vcf.gz"),
                    #"indexURL": url_for('serve_file', genome=genome, filename="var_only.ann.G2P-SOL_Pepper_CC_Geno_data_2023.vcf.gz.tbi"),
                    "url": url_for('in_memory_vcf'),
                    "height": 50,
                    "visibilityWindow": 10000000,
                },
            ],
       },
        "roi": [
            {
                "name": "ROI set 2",
                "color": "rgba(3,52,249,0.25)",
                "features": [
                    {
                        "chr": "1",
                        "start": 60936999,
                        "end": 60936999,
                    },
                ],
            },
        ],
        "locus": "1:60490000-60520000",
    }
"""