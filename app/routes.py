#!/usr/bin/env python

import app.path
from flask_app import (
    app,
    get_lock,
    running_pr,
    task_queue,
    master_lock,
    task_status,
    cancelled_tasks,
)

from app.dirs import data_dir

import os
import re
import json
import uuid
import pysam
import tabix
import shutil
import signal
import numpy as np
import pandas as pd

from time import time
from sklearn.decomposition import PCA
from scipy.stats import boxcox, zscore

from flask import (
    flash,
    jsonify,
    request,
    session,
    url_for,
    redirect,
    render_template,
)

from app.ajax import (
    group,
    rm_markers,
    contextmenu,
    fetch_mplot,
    fetch_bblplot,
)

from app.igv_conf import igv_am

from app.serve import (
    favicon,
    serve_ld,
    serve_icon,
    serve_assoc,
    serve_image,
)

from modules.am import am_pipeline
from modules.tmp_clean import clean
from modules.safe_col import safe_col

from modules.assoc.assoc_qc import mdbin

from modules.assoc.kq import kq_correct

from modules.forms import (
    AddForm,
    EditForm,
    AssocQCform,
    AssocResForm,
    TraitViewForm,
    FuncAnalysisForm,
)

from modules.forms import models as models_dict

from modules.load_input import load_gdata, load_config

models_dict = dict(models_dict)

###############################################################################
### Long task progress
###############################################################################


@app.route("/status/<task_id>", methods=["GET"])
def get_status(task_id):
    # Remove tasks finished more than a week ago
    with master_lock:
        status_lock = get_lock("__task_status__")
    with status_lock:
        for _id in task_status.keys():
            if (
                _id != "specs"
                and "end" in task_status[_id].keys()
                and time() - task_status[_id]["end"] > 604800
            ):
                _ = task_status.pop(_id)

        if task_id not in task_status.keys():
            with master_lock:
                queue_lock = get_lock("__task_queue__")
            with queue_lock:
                if task_id not in [t[0] for t in list(task_queue.queue)]:
                    return jsonify(""), 404
                else:
                    return jsonify(""), 201
        else:
            progress = task_status[task_id]["progress"]
            task_status[task_id]["progress"] = ""
            code = task_status[task_id]["code"]
            if code == 200 or code == 500:
                end = time()
                task_status[task_id]["end"] = end
                duration = end - task_status[task_id]["start"]
                return jsonify([progress, round(duration)]), code
            return jsonify(progress), code


@app.route("/remove_task/<task_id>", methods=["POST"])
def remove_task(task_id):
    with master_lock:
        status_lock = get_lock("__task_status__")

    (pr, task_spec) = running_pr.get(task_id, (None, None))
    if pr and pr.poll() is None:
        with status_lock:
            try:
                os.killpg(os.getpgid(pr.pid), signal.SIGTERM)
            except ProcessLookupError:
                pass
            task_status[task_id]["code"] = 410
            if task_spec in task_status["specs"]:
                task_status["specs"].remove(task_spec)
    else:
        cancelled_tasks.append(task_id)
        with status_lock:
            task_status[task_id] = {"progress": "", "code": 410}

    return jsonify("")


###############################################################################
### Index page
###############################################################################


@app.route("/")
def index():
    # Remove files older than 7 days from tmp_dir and output/
    clean()

    genomes = sorted(
        [
            genome
            for genome in os.listdir(data_dir)
            if os.path.exists(os.path.join(data_dir, genome, "stats", "init"))
        ]
    )
    links = {
        genome: {
            "gv": url_for("fa_opts", genome=genome),
            "ts": url_for("ts", genome=genome),
            "assoc": url_for("assoc", genome=genome),
            "edit": url_for("edit", genome=genome),
        }
        for genome in genomes
    }

    update_failed = os.path.exists("updatefailed")
    if update_failed:
        os.remove("updatefailed")

    return render_template(
        "index.html",
        genomes=genomes,
        links=links,
        genome="",
        update_failed=update_failed,
    )


###############################################################################
### Error page
###############################################################################


@app.route("/error/<error>")
def error(error):
    return render_template("error.html", error=error)


###############################################################################
### Add new/edit dataset
###############################################################################


@app.route("/new_dataset", methods=["GET"])
def new():
    form = AddForm()
    genomes = sorted(
        [
            genome
            for genome in os.listdir(data_dir)
            if os.path.isdir(os.path.join(data_dir, genome))
            and not os.path.exists(os.path.join(data_dir, genome, "stats", "init"))
        ]
    )
    genomes = [(url_for("setup_new", genome=genome), genome) for genome in genomes]
    form.name.choices = [(url_for("new"), "Select dataset")] + genomes
    return render_template("new.html", form=form, genome="", edit=False, action="")


@app.route("/new_dataset/<genome>", methods=["GET"])
def setup_new(genome):
    form = AddForm()
    avail = sorted(
        [
            genome
            for genome in os.listdir(data_dir)
            if os.path.isdir(os.path.join(data_dir, genome))
            and not os.path.exists(os.path.join(data_dir, genome, "stats", "init"))
        ]
    )
    avail = [(url_for("setup_new", genome=genome), genome) for genome in avail]
    form.name.choices = [(url_for("new"), "Select dataset")] + avail
    form.name.data = url_for("setup_new", genome=genome)

    form.setup_from.choices = [""] + sorted(
        [
            genome
            for genome in os.listdir(data_dir)
            if os.path.isdir(os.path.join(data_dir, genome))
            and os.path.exists(os.path.join(data_dir, genome, "stats", "init"))
        ]
    )

    genome_dir = os.path.join(data_dir, genome)
    exts = {
        "vcf": (".vcf", ".vcf.gz", ".bcf", ".bcf.gz"),
        "ann": (".gff", ".gtf", ".gff3", ".gff.gz", ".gtf.gz", ".gff3.gz"),
    }
    exts.update(
        dict.fromkeys(
            ["seq", "cds", "pep"],
            (
                ".fa",
                ".fasta",
                ".fna",
                ".faa",
                ".fa.gz",
                ".fasta.gz",
                ".fna.gz",
                ".faa.gz",
            ),
        )
    )
    exts.update(
        dict.fromkeys(
            ["phe", "exp", "run", "samp", "info", "sev"],  # "proj",
            (".txt", ".tsv", ".txt.gz", ".tsv.gz"),
        )
    )
    for ftype, ext in exts.items():
        form[ftype].choices = [""] + [
            f for f in os.listdir(genome_dir) if f.endswith(ext)
        ]
    form.bam.choices = [f for f in os.listdir(genome_dir) if f.endswith(".bam")]
    prot_db_path = os.path.join(genome_dir, "prot_db")
    if os.path.exists(prot_db_path):
        form.db_prefix.choices = [""] + list(
            set(
                match.group(1)
                for fn in os.listdir(prot_db_path)
                if (match := re.match(r"(.+)(\.\d+)?\.phr$", fn))
            )
        )
        # form.db_fasta.choices = [""] + list(
        #  set(
        #    match.group(1)
        #    for fn in os.listdir(prot_db_path)
        #    if (match := re.match(
        #      r"((.+)\.(fasta|fa)(.gz)?)$",
        #      fn))))

    return render_template(
        "new.html",
        form=form,
        genome=genome,
        edit=False,
        action=url_for("setup", genome=genome),
    )


@app.route("/<genome>/setup", methods=["POST"])
def setup(genome):
    form = AddForm(request.form)
    task_id = str(uuid.uuid4())
    task_spec = ("setup", genome)

    task_params = [("Task:", "Setting up a new genome"), ("Genome:", genome)]
    param_keys = [
        "setup_from",
        "vcf",
        "ann",
        "seq",
        "phe",
        "exp",
        "run",
        "samp",
        "info",
        "snpeff",
        "cds",
        "pep",
        "sev",
        "db_prefix",
        "concat",
    ]
    for key in param_keys:
        task_params.append((form._fields[key].label.text, str(form.data[key])))
        task_params[0] = ("From dataset:", form.setup_from.data)

        if form.setup_from.data:
            task_params = [p for i, p in enumerate(task_params) if i in (0,1,4,7,9)]

    task = [
        "python3",
        os.path.join("modules", "setup.py"),
        genome,
        json.dumps(form.data),
        task_id,
    ]

    log_dir = os.path.join(data_dir, genome, "logs")
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    # Make sure the same task isn't running or queued
    with master_lock:
        status_lock = get_lock("__task_status__")
    with status_lock:
        if task_spec in task_status["specs"] or task_spec in [
            t[1] for t in list(task_queue.queue)
        ]:
            return render_template(
                "error.html", error="Dataset is already being processed."
            )
        else:
            #            queue_lock = get_lock("__task_queue__")
            #            with queue_lock:
            task_queue.put((task_id, task_spec, task_params, task))
            return render_template(
                "loading.html",
                genome=genome,
                task="setup",
                task_id=task_id,
                task_spec=task_spec,
                form="",
            )


@app.route("/<genome>/Edit", methods=["GET", "POST"])
def edit(genome):
    errors = []
    config = load_config(genome)
    form = EditForm()

    if request.method == "POST":
        phe_file = request.files["phe"]
        if phe_file:
            phe_path = os.path.join(data_dir, genome, config["phe"])
            with master_lock:
                lock = get_lock(phe_path)
            with lock:
                phe_df = pd.read_csv(phe_path, sep="\t", header=0, index_col=0)

                new_phe = pd.read_csv(
                    phe_file,
                    header=0,
                    index_col=0,
                    sep="," if phe_file.filename.endswith(".csv") else "\t",
                )

                if len(new_phe.columns) == 0:
                    errors.append("<b>No columns found in file.</b>")
                else:
                    safe_cols = [safe_col(col) for col in new_phe.columns]
                    new_phe = new_phe.set_axis(safe_cols, axis=1)
                    col_dup = new_phe.columns.intersection(phe_df.columns)
                    if not col_dup.empty:
                        col_dup = ", ".join(col_dup.tolist())
                        errors.append(f"<b>Duplicate traits:</b><br>{col_dup}")
                    else:
                        new_phe = pd.concat([phe_df, new_phe], axis=1)
                        new_phe.to_csv(phe_path, sep="\t")
                        # Update samples info to contain all samples
                        samples = new_phe.index
                        samp_path = os.path.join(data_dir, genome, config["samp"])
                        with master_lock:
                            lock_2 = get_lock(samp_path)
                        with lock_2:
                            samp_df = pd.read_csv(
                                samp_path, sep="\t", header=0, index_col=0
                            )
                            miss_samples = [
                                s for s in samples if s not in samp_df.index
                            ]
                            miss_df = pd.DataFrame(
                                miss_samples,
                                index=miss_samples,
                                columns=[samp_df.columns[0]],
                            )
                            (
                                pd.concat([samp_df, miss_df])
                                .rename_axis("Acc")
                                .to_csv(samp_path, sep="\t")
                            )
                            errors = 0

        samp_file = request.files["samp"]
        if samp_file:
            samp_path = os.path.join(data_dir, genome, config["samp"])
            with master_lock:
                lock = get_lock(samp_path)
            with lock:
                samp_df = pd.read_csv(samp_path, sep="\t", header=0, index_col=0)

                new_samp = pd.read_csv(
                    samp_file,
                    header=0,
                    index_col=0,
                    sep="," if samp_file.filename.endswith(".csv") else "\t",
                )

                if len(new_samp.columns) == 0:
                    errors.append("<b>No columns found in file.</b>")
                else:
                    new_samp = pd.concat([samp_df, new_samp], axis=1)
                    new_samp.to_csv(samp_path, sep="\t")
                    errors = 0

        cof_file = request.files["cof"]
        if cof_file:
            cof_path = os.path.join(data_dir, genome, "covariates", cof_file.filename)
            if os.path.exists(cof_path):
                errors.append("Covariate file with the same name already exists.")
            else:
                (pd
                    .read_csv(
                        cof_file,
                        header=0,
                        index_col=0,
                        sep="," if cof_file.filename.endswith(".csv") else "\t",
                    )
                    .to_csv(
                        cof_path,
                        sep="\t",
                    )
                )

    return render_template(
        "new.html",
        genome=genome,
        form=form,
        edit=True,
        errors=errors,
        action=url_for("edit", genome=genome),
    )


###############################################################################
### Trait/sample view
###############################################################################


@app.route("/<genome>/Trait_sample_view", methods=["GET", "POST"])
def ts(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = TraitViewForm()
    errors = []
    phe_path = os.path.join(data_dir, genome, config["phe"])

    if request.method == "POST":
        form = TraitViewForm(request.form)
        params = form.data
        traits = params["trait"]
        for group in params["trait_group"]:
            traits += [t for t in json.loads(group) if t not in traits]

        if params["t_name"]:
            with master_lock:
                lock = get_lock(phe_path)
            with lock:
                phe_df = pd.read_csv(phe_path, sep="\t", header=0, index_col=0)

                samples = json.loads(params["sample_sel"]) or phe_df.index.tolist()
                if params["transform"] == "pca":
                    pca = PCA(n_components=params["num_pc"])
                    phe = phe_df[traits].dropna()
                    columns = [
                        params["t_name"] + f"_PC{i}"
                        for i in range(1, params["num_pc"] + 1)
                    ]
                    transforms = pd.DataFrame(
                        pca.fit_transform(phe), index=phe.index, columns=columns
                    )
                else:
                    transforms = list()
                    for trait in traits:
                        raw_phe = phe_df.loc[samples, [trait]].dropna()

                        if params["transform"] == "log":
                            transform = np.log(raw_phe)
                        elif params["transform"] == "bin":
                            transform = (raw_phe >= float(params["binthr"])).astype(int)
                        elif params["transform"] == "zs":
                            transform = pd.DataFrame(zscore(raw_phe))
                        elif params["transform"] == "outl":
                            transform = raw_phe.mask(
                                np.abs(zscore(raw_phe)) > float(params["num_sd"])
                            )
                        elif params["transform"] == "thr":
                            transform = raw_phe.mask(
                                (
                                    raw_phe
                                    < (
                                        params["thr_min"]
                                        if params["thr_min"] is not None
                                        else -np.inf
                                    )
                                )
                                | (
                                    raw_phe
                                    > (
                                        params["thr_max"]
                                        if params["thr_max"] is not None
                                        else np.inf
                                    )
                                )
                            )
                        elif params["transform"] == "bc":
                            transform[params["trait_sgl"]] = boxcox(
                                transform[params["trait_sgl"]]
                            )[0]
                        elif params["transform"] == "corr":
                            phe_corr = kq_correct(
                                genome,
                                raw_phe[params["trait_sgl"]],
                                None,
                                params["k"],
                                int(params["q_num"]) if params["q"] else 0,
                                params["kq_from"],
                                params["cov"],
                            )
                            transform = pd.DataFrame(phe_corr)

                        else:
                            transform = pd.DataFrame(raw_phe)
                        transform = transform.set_axis(
                            [trait + "_" + params["t_name"]],
                            axis=1)
                        transforms.append(transform)
                    transforms = pd.concat(transforms, axis=1)

                safe_cols = [safe_col(col) for col in transforms.columns]
                transforms = transforms.set_axis(safe_cols, axis=1)
                col_dup = transforms.columns.intersection(phe_df.columns)
                if not col_dup.empty:
                    col_dup = ", ".join(col_dup.tolist())
                    errors.append(f"<b>Duplicate traits after processing:</b><br>{col_dup}")
                else:
                    (pd.concat([phe_df, transforms], axis=1).to_csv(phe_path, sep="\t"))
                    errors = 0

            flash(errors)
            return redirect(url_for("ts", genome=genome))

        else:
            form.trait_sgl.data = form.trait.data[0]

    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        with open(phe_path, "r") as phe_in:
            phe = [line.strip("\n").split("\t") for line in phe_in if line.strip()]
    if phe:
        form.trait_sgl.choices = phe[0][1:]
        form.trait.choices = phe[0][1:]

    vcf_path = os.path.join(data_dir, genome, config["vcf"])
    vcf = pysam.VariantFile(vcf_path)
    gt_samples = list(vcf.header.samples)
    samp_path = os.path.join(data_dir, genome, config["samp"])

    with master_lock:
        lock = get_lock(samp_path)
    with lock:
        samp_df = pd.read_csv(samp_path, sep="\t", header=0, index_col=0).fillna("")
    samp_info = [["ID"] + samp_df.columns.tolist()] + samp_df.to_records().tolist()
    samples = samp_df[[samp_df.columns[0]]].to_records().tolist()
    form.sample.choices = samples
    form.sample_group.choices = [
        (json.dumps(group[0]), group[1]) for group in gdata["group"].get("sample", [])
    ]
    form.trait_group.choices = [
        [json.dumps(group[0]), group[1]] for group in gdata["group"].get("trait", [])
    ]

    markers_root = os.path.join(data_dir, genome, "assoc")
    for m in os.listdir(markers_root):
        init_path = os.path.join(markers_root, m, "stats", "init")
        if os.path.exists(init_path):
            form.kq_from.choices.append(m)
    form.cov.choices = [
        fn for fn in os.listdir(os.path.join(data_dir, genome, "covariates"))
    ]

    return render_template(
        "trt/trt.html",
        genome=genome,
        form=form,
        phe=phe,
        gt_samples=gt_samples,
        samp_info=samp_info,
        errors=errors,
    )


###############################################################################
### Functional options
###############################################################################


@app.route("/<genome>/Gene_variant_view", methods=["GET", "POST"])
def fa_opts(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = FuncAnalysisForm()

    markers_root = os.path.join(data_dir, genome, "assoc")
    form.dataset.choices += [
        (d, d)
        for d in os.listdir(markers_root)
        if os.path.isdir(os.path.join(markers_root, d))
        and os.path.exists(os.path.join(markers_root, d, "stats", "init"))
    ]

    form.chrom.choices = gdata["chr_names"]
    samp_path = os.path.join(data_dir, genome, config["samp"])
    with master_lock:
        lock = get_lock(samp_path)
    with lock:
        samples = (
            pd.read_csv(samp_path, sep="\t", header=0, index_col=0, usecols=[0, 1])
            .to_records()
            .tolist()
        )
    form.sample.choices = samples

    form.sample_group.choices = [
        (json.dumps(group[0]), group[1]) for group in gdata["group"].get("sample", [])
    ]
    form.gene_group.choices = [
        (json.dumps(group[0]), group[1]) for group in gdata["group"].get("gene", [])
    ]
    form.coord_group.choices = [
        (json.dumps(group[0]), group[1]) for group in gdata["group"].get("coord", [])
    ]

    for key, value in request.args.items():
        form[key].data = value

    return render_template(
        "fa/opts.html",
        genome=genome,
        form=form,
        anns=None,
        genes=None,
        is_res=False,
        id_fmt=gdata["id_fmt"],
    )


###############################################################################
### Functional analysis results
###############################################################################


@app.route("/<genome>/Gene_variant_view/Results", methods=["GET", "POST"])
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
    samp_path = os.path.join(data_dir, genome, config["samp"])
    with master_lock:
        lock = get_lock(samp_path)
    with lock:
        samples = (
            pd.read_csv(samp_path, sep="\t", header=0, index_col=0, usecols=[0, 1])
            .to_records()
            .tolist()
        )
    form.sample.choices = samples

    form.sample_group.choices = [
        [json.dumps(group[0]), group[1]] for group in gdata["group"].get("sample", [])
    ]
    form.gene_group.choices = [
        [json.dumps(group[0]), group[1]] for group in gdata["group"].get("gene", [])
    ]
    form.coord_group.choices = [
        (json.dumps(group[0]), group[1]) for group in gdata["group"].get("coord", [])
    ]
    form.exprun_group.choices = [
        [json.dumps(group[0]), group[1]] for group in gdata["group"].get("exprun", [])
    ]
    #    form.trait_group.choices = [
    #      [json.dumps(group[0]), group[1]]
    #      for group in gdata["group"].get("trait", [])]

    (gts, groups, genes, anns, errors) = am_pipeline(genome, form)

    # Get attribute names of additional information for genes
    if config["info"]:
        info_path = os.path.join(data_dir, genome, config["info"])
        with open(info_path, "r") as info_in:
            gene_attrs = next(info_in).strip().split("\t")[1:]
    else:
        gene_attrs = []

    genes.sort(key=lambda g: (gdata["chr_names"].index(g[2]), g[3], g[4]))
    gene_attrs = [
        "Gene ID",
        "Transcript ID",
        "Chromosome",
        "Start",
        "End",
        "Length",
        "Strand",
        "Variants",
    ] + gene_attrs
    genes = [gene_attrs] + genes

    if config["exp"]:
        exp_path = os.path.join(data_dir, genome, config["exp"])

        runs_path = os.path.join(data_dir, genome, config["run"])
        with open(runs_path, "r") as runs_in:
            runs_info = [line.strip().split("\t") for line in runs_in]
            form.exp_col.choices = [
                ("exp-by_" + str(idx), col) for idx, col in enumerate(runs_info[0])
            ]

        # if config["proj"]:
        #    proj_path = os.path.join(
        #      data_dir,
        #      genome,
        #      config["proj"])
        #    with open(proj_path, "r") as proj_in:
        #        proj_info = [
        #          line.strip().split("\t")
        #          for line in proj_in]
        # else:
        #    proj_info = ""

        exp_data = list()
        tb = tabix.open(exp_path)
        if genes:
            # Skip header line
            for gene in genes[1:]:
                try:
                    for t_entry in tb.query(gene[1], 0, 1):
                        #if t_entry[0] != gene[1]:
                        #    continue
                        exp_data.append([t_entry[0]] + t_entry[2:])
                except tabix.TabixError:
                    # Gene not found, ignore
                    pass
    else:
        exp_data = None
        runs_info = ""
        # proj_info = ""

    # try:
    #    with open("test_colors.json", "r") as in_colors:
    #        colors = json.load(in_colors)
    # except:
    #    colors = ""

    phe_path = os.path.join(data_dir, genome, config["phe"])
    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        with open(phe_path, "r") as phe_in:
            phe = [line.strip("\n").split("\t") for line in phe_in if line.strip()]
    # if config["enr"]:
    #    form.enr_col.choices = json.loads(config["enr"])
    igv_conf = igv_am(genome)

    return render_template(
        "fa/res.html",
        genome=genome,
        gts=gts,
        groups=groups,
        genes=genes,
        anns=anns,
        errors=errors,
        form=form,
        is_res=True,
        id_fmt=gdata["id_fmt"],
        phe=phe,
        runs_info=runs_info,
        # proj_info = proj_info,
        exp_data=exp_data,
        # colors    = colors,
        igv_conf=igv_conf,
    )


###############################################################################
### Annotate severity
###############################################################################


@app.route("/<genome>/Annotate_severity", methods=["POST"])
def ann_sev(genome):
    # prot_db = os.path.join(
    #  "res",
    #  "prot_db",
    #  "prot")
    prot_dir = os.path.join(data_dir, genome, "prot_db")
    if not os.path.exists(prot_dir):
        return render_template("error.html", error="Missing required input files.")

    form = FuncAnalysisForm(request.form)
    markers = json.loads(form.data["ann_sev"])

    marker_data = list()
    for marker in markers:
        ch, pos = marker.replace(",", "").rsplit(":", 1)
        pos = int(pos.split(".")[0])
        marker_data.append((ch, pos))

    task_id = str(uuid.uuid4())
    task_spec = ("sev", genome, marker_data)

    task_params = [
        ("Task:", "Annotating marker severity"),
        ("Genome:", genome),
        ("Markers:", json.dumps(markers)),
    ]

    task = [
        "python3",
        os.path.join("modules", "sev.py"),
        genome,
        json.dumps(marker_data),
        task_id,
    ]

    task_queue.put((task_id, task_spec, task_params, task))
    form.ann_sev.data = ""
    return render_template(
        "loading.html",
        genome=genome,
        task="sev",
        task_id=task_id,
        task_spec=task_spec,
        warning="",
        form=form.data,
    )


###############################################################################
### Association analysis redirect
###############################################################################


@app.route("/<genome>/Association_analysis", methods=["GET"])
def assoc(genome):
    assoc_dir = os.path.join(data_dir, genome, "assoc")
    if not os.path.exists(assoc_dir):
        os.mkdir(assoc_dir)
    if not os.listdir(assoc_dir):
        return redirect(url_for("assoc_qc", genome=genome))
    elif not any(
        os.path.exists(os.path.join(assoc_dir, markers, "res"))
        for markers in os.listdir(assoc_dir)
    ):
        return redirect(url_for("assoc_opts", genome=genome))
    else:
        return redirect(url_for("assoc_res", genome=genome))


###############################################################################
### Choose filters for defining a new marker set
###############################################################################


@app.route("/<genome>/Association_analysis/QC", methods=["GET"])
def assoc_qc(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = AssocQCform()
    stats_dir = os.path.join(data_dir, genome, "stats")

    hist_data = dict()
    for metric in ("smiss", "vmiss", "shet", "vhet", "maf", "hwe", "dens"):
        hist_path = os.path.join(stats_dir, metric)
        with open(hist_path, "r") as hist_in:
            hist_data[metric] = json.load(hist_in)

    with open(os.path.join(stats_dir, "cnts"), "r") as cnts_in:
        samp_num, var_num = json.load(cnts_in)

    markers_root = os.path.join(data_dir, genome, "assoc")
    saved_sets = list()
    for d in os.listdir(markers_root):
        init_path = os.path.join(markers_root, d, "stats", "init")
        if os.path.exists(init_path):
            saved_sets.append(d.casefold())
            hap_fn = os.path.join(markers_root, d, "stats", "is_hap")
            if not os.path.exists(hap_fn):
                form.filt_from.choices.append((d, d))

    form.coord.choices = [
        [json.dumps(group[0]), group[1]] for group in gdata["group"].get("coord", [])
    ]

    form.chr.choices = gdata["chr_names"]

    form.sample_group.choices = [
        (json.dumps(group[0]), group[1]) for group in gdata["group"].get("sample", [])
    ]

    return render_template(
        "assoc/qc.html",
        genome=genome,
        form=form,
        hist_data=hist_data,
        samp_num=samp_num,
        var_num=var_num,
        mdbin=mdbin,
        chr_names=gdata["chr_names"],
        chr_sizes=gdata["chr_sizes"],
        chr_offsets=gdata["chr_offsets"],
        concat=json.loads(config["concat"]),
        saved_sets=saved_sets,
    )


###############################################################################
### Apply filters to create new marker set, redirect
###############################################################################


@app.route("/<genome>/Association_analysis/Filter", methods=["POST"])
def assoc_filt(genome):
    form = AssocQCform(request.form)
    markers = form.data["name"]
    task_id = str(uuid.uuid4())
    task_spec = ("filt_qc", genome, markers)

    task_params = [("Task:", "Filtering a marker set"), ("Genome:", genome)]
    param_keys = [
        "name",
        "filt_from",
        "sample_group",
        "coord",
        "smiss",
        "vmiss",
        "shet",
        "vhet",
        "maf",
        "hwe",
        "prune",
        "window",
        "imp",
        "hap",
    ]
    if form.data["hap"]:
        if form.data["btype"] == "ld":
            param_keys += ["btype", "bsize"]
        else:
            param_keys += ["btype", "bspan", "bstep"]
    for key in param_keys:
        task_params.append((form._fields[key].label.text, str(form.data[key])))
    task_params.insert(6, ("Missingness per:",))
    task_params.insert(9, ("Heterozygosity per:",))
    task_params.insert(14, ("Pruning:",))

    task = [
        "python3",
        os.path.join("modules", "assoc", "filt_qc.py"),
        genome,
        json.dumps(form.data),
        task_id,
    ]

    set_dir = os.path.join(data_dir, genome, "assoc", markers)
    if os.path.exists(os.path.join(set_dir, "stats", "init")):
        return render_template(
            "error.html", error="Set name is taken. Choose another name."
        )
    elif os.path.exists(set_dir):
        shutil.rmtree(set_dir)

    # Make sure the same task isn't running or queued
    with master_lock:
        status_lock = get_lock("__task_status__")
    with status_lock:
        if task_spec in task_status["specs"] or task_spec in [
            t[1] for t in list(task_queue.queue)
        ]:
            return render_template(
                "error.html", error="Set name is taken. Choose another name."
            )
        else:
            task_queue.put((task_id, task_spec, task_params, task))
            return render_template(
                "loading.html",
                genome=genome,
                task="filt",
                task_id=task_id,
                task_spec=task_spec,
                warning="",
                form="",
            )


###############################################################################
### Choose parameters for association analysis
###############################################################################


@app.route("/<genome>/Association_analysis/Options", methods=["GET"])
def assoc_opts(genome):
    gdata = load_gdata(genome)
    config = load_config(genome)
    form = AssocResForm()

    phe_path = os.path.join(data_dir, genome, config["phe"])
    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        with open(phe_path, "r") as phe_in:
            phe = [line.strip("\n").split("\t") for line in phe_in if line.strip()]

    form.trait.choices = phe[0][1:]
    form.by_phe.choices = phe[0][1:]

    markers_root = os.path.join(data_dir, genome, "assoc")
    haps = dict()
    for m in os.listdir(markers_root):
        init_path = os.path.join(markers_root, m, "stats", "init")
        if os.path.exists(init_path):
            hap_fn = os.path.join(markers_root, m, "stats", "is_hap")
            is_hap = os.path.exists(hap_fn)
            haps[m] = is_hap

            form.markers_sgl.choices.append(m)
            form.kq_from.choices.append(m)

    form.cov.choices = [
        fn for fn in os.listdir(os.path.join(data_dir, genome, "covariates"))
    ]

    form.trait_group.choices = [
        [json.dumps(group[0]), group[1]] for group in gdata["group"].get("trait", [])
    ]
    form.sample_group.choices = [
        (json.dumps(group[0]), group[1]) for group in gdata["group"].get("sample", [])
    ]
    form.plot.choices += [("kin", "Kinship matrix"), ("pca", "Principal components")]
    concat = json.loads(config["concat"])
    form.chr.choices = [c for c in gdata["chr_names"] if c not in concat]

    return render_template(
        "assoc/opts.html",
        genome=genome,
        form=form,
        phe_upl=True,
        rm_markers=True,
        haps=haps,
        mdbin=mdbin,
        chr_names=gdata["chr_names"],
        chr_sizes=gdata["chr_sizes"],
        chr_offsets=gdata["chr_offsets"],
        phe=phe,
        concat=concat,
        errors=None,
    )


###############################################################################
### Run association analysis
###############################################################################


@app.route("/<genome>/Association_analysis/Run", methods=["POST"])
def assoc_run(genome):
    errors = list()
    form = AssocResForm(request.form)
    data = form.data

    res_dir = os.path.join(data_dir, genome, "assoc", data["markers_sgl"], "res")
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    task_params = [("Task:", "Running association analysis"), ("Genome:", genome)]
    param_keys = [
        "markers_sgl",
        "model_sgl",
        "trait",
        "trait_group",
        "k",
        "q",
        "q_num",
        "kq_from",
        "cov",
        "overwrite",
    ]
    if data["model_sgl"] in ("rf", "rfc"):
        param_keys += ["mtry", "ntree"]
    elif data["model_sgl"] in ("B", "C"):
        param_keys += ["b_burn", "b_chain", "b_thin"]
    for key in param_keys:
        task_params.append((form._fields[key].label.text, str(data[key])))
    task_params[3] = ("Model:", models_dict[data["model_sgl"]])

    if data["model_sgl"] != "mtmlm":
        for group in data["trait_group"]:
            group, group_name = json.loads(group)
            for trait in group:
                if trait not in data["trait"]:
                    data["trait"].append(trait)

        if not data["overwrite"]:
            existing_traits = list()
            for trait in data["trait"][:]:
                # Skip existing runs
                trait_out = os.path.join(
                    data_dir,
                    genome,
                    "assoc",
                    data["markers_sgl"],
                    "res",
                    f"{data['model_sgl']}.{trait}.feather",
                )
                if os.path.exists(trait_out):
                    data["trait"].remove(trait)
                    existing_traits.append(trait)
            if existing_traits:
                errors.append("Results already exist for:")
                errors.append(", ".join(existing_traits))

            # Skip queued runs
            queued_traits = list()
            for task_spec in task_status["specs"]:
                if (
                    task_spec[0] == "run_assoc"
                    and task_spec[1] == genome
                    and task_spec[2] == data["markers_sgl"]
                    and task_spec[3] == data["model_sgl"]
                ):
                    for trait in task_spec[4]:
                        if trait in data["trait"]:
                            data["trait"].remove(trait)
                            queued_traits.append(trait)
            for t in list(task_queue.queue):
                task_spec = t[1]
                if (
                    task_spec[0] == "run_assoc"
                    and task_spec[1] == genome
                    and task_spec[2] == data["markers_sgl"]
                    and task_spec[3] == data["model_sgl"]
                ):
                    for trait in task_spec[4]:
                        if trait in data["trait"]:
                            data["trait"].remove(trait)
                            queued_traits.append(trait)

            if queued_traits:
                errors.append("Results are already being generated for:")
                errors.append(", ".join(queued_traits))

        if data["trait"]:
            data["trait"] = sorted(data["trait"])

            task_id = str(uuid.uuid4())
            task_spec = (
                "run_assoc",
                genome,
                data["markers_sgl"],
                data["model_sgl"],
                data["trait"],
            )

            task = [
                "python3",
                os.path.join("modules", "assoc", "run_assoc.py"),
                genome,
                json.dumps(data),
                task_id,
            ]

            task_queue.put((task_id, task_spec, task_params, task))
            return render_template(
                "loading.html",
                genome=genome,
                task="assoc",
                task_id=task_id,
                task_spec=task_spec,
                errors=errors,
                form=data,
            )

        else:
            errors = [
                "Results already exist or are being generated for all selected traits."
            ]
            flash(errors)
            return redirect(url_for("assoc_res", genome=genome))
    else:
        if data["trait"]:
            pass
            # errors = ["To use multi-trait MLM, traits must be grouped."]
        else:
            data["trait_group"] = [json.loads(g) for g in data["trait_group"]]
            existing_groups = list()
            for group, group_name in data["trait_group"][:]:
                # Skip existing runs
                trait_out = os.path.join(
                    data_dir,
                    genome,
                    "assoc",
                    data["markers_sgl"],
                    "res",
                    f"{data['model_sgl']}.{group_name}.feather",
                )
                if not data["overwrite"]:
                    if os.path.exists(trait_out):
                        data["trait_group"].remove([group, group_name])
                        existing_groups.append(group_name)
            if existing_groups:
                errors.append("Results already exist for:")
                errors.append(", ".join(existing_groups))
            # Skip queued runs
            queued_groups = list()
            for task_spec in task_status["specs"]:
                if (
                    task_spec[0] == "run_assoc"
                    and task_spec[1] == genome
                    and task_spec[2] == data["markers_sgl"]
                    and task_spec[3] == data["model_sgl"]
                ):
                    for group in task_spec[4]:
                        if group in data["trait_group"]:
                            data["trait_group"].remove(group)
                            queued_groups.append(group[1])
            if queued_groups:
                errors.append("Results are already being generated for:")
                errors.append(", ".join(queued_groups))

            if data["trait_group"]:
                task_id = str(uuid.uuid4())
                task_spec = (
                    "run_assoc",
                    genome,
                    data["markers_sgl"],
                    data["model_sgl"],
                    data["trait_group"],
                )
                task = [
                    "python3",
                    os.path.join("modules", "assoc", "run_assoc.py"),
                    genome,
                    json.dumps(data),
                    task_id,
                ]

                task_queue.put((task_id, task_spec, task))
                return render_template(
                    "loading.html",
                    genome=genome,
                    task="assoc",
                    task_id=task_id,
                    task_spec=task_spec,
                    errors=errors,
                    form=data,
                )
            else:
                errors = ["Results already exist for all selected groups."]
                flash(errors)
                return redirect(url_for("assoc_res", genome=genome))


###############################################################################
### Association analysis results
###############################################################################


@app.route("/<genome>/Association_analysis/Results", methods=["GET", "POST"])
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
        form["trait"].data = form["trait"].data
        form["trait_group"].data = [
            json.dumps(json.loads(group.replace("'", '"'))[0])
            for group in form["trait_group"].data
        ]

    markers_root = os.path.join(data_dir, genome, "assoc")

    form.markers.choices = [
        (d, d)
        for d in os.listdir(markers_root)
        if os.path.isdir(os.path.join(markers_root, d))
        and os.path.exists(os.path.join(markers_root, d, "stats", "init"))
    ]
    form.markers_sgl.choices = [("", "Raw VCF")] + form.markers.choices[:]

    analysed_traits = list()
    analysed_models = list()
    for marker_set in form.markers.choices:
        marker_set = marker_set[1]
        assoc_root = os.path.join(data_dir, genome, "assoc", marker_set, "res")
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

    phe_path = os.path.join(data_dir, genome, config["phe"])
    with master_lock:
        lock = get_lock(phe_path)
    with lock:
        with open(phe_path, "r") as phe_in:
            phe_columns = next(phe_in).rstrip().split("\t")[1:]
    form.trait.choices = [t for t in analysed_traits if t in phe_columns]

    form.models.choices = [(m, models_dict[m]) for m in analysed_models]

    form.trait_group.choices = [
        [json.dumps(group[0]), group[1]] for group in gdata["group"].get("trait", [])
    ]

    form.hgl_coord_group.choices = [
        [json.dumps(group[0]), group[1]] for group in gdata["group"].get("coord", [])
    ]

    form.hgl_gene_group.choices = [
        (json.dumps(group[0]), group[1]) for group in gdata["group"].get("gene", [])
    ]

    chroms = gdata["chr_names"]
    sizes = gdata["chr_sizes"]
    offsets = gdata["chr_offsets"]
    concat = json.loads(config["concat"])

    is_hap = dict()
    for marker_set in form.markers.choices:
        marker_set = marker_set[1]
        hap_fn = os.path.join(markers_root, marker_set, "stats", "is_hap")
        is_hap[marker_set] = os.path.exists(hap_fn)

    return render_template(
        "assoc/res.html",
        genome=genome,
        form=form,
        errors=errors,
        chroms=chroms,
        sizes=sizes,
        offsets=offsets,
        models=models_dict,
        concat=concat,
        is_hap=is_hap,
    )


@app.route("/<genome>/Association_analysis/QTLs", methods=["POST"])
def assoc_qtls(genome):
    form = AssocResForm(request.form)
    data = form.data
    errors = list()

    task_id = str(uuid.uuid4())
    task_spec = ("extract_qtls", genome, task_id)

    task_params = [("Task:", "Extracting significant QTLs"), ("Genome:", genome)]
    param_keys = ["qtl_runs", "markers_sgl", "plot_thr", "window_ld", "r2_ld"]
    for key in param_keys:
        task_params.append((form._fields[key].label.text, str(data[key])))
    task_params[2] = (
        "Results to process:",
        json.dumps(
            [
                ", ".join((r[0], models_dict[r[1]], r[2]))
                for r in json.loads(task_params[2][1])
            ]
        ),
    )

    task = [
        "python3",
        os.path.join("modules", "assoc", "extract_qtls.py"),
        genome,
        json.dumps(data),
        task_id,
    ]

    task_queue.put((task_id, task_spec, task_params, task))
    return render_template(
        "loading.html",
        genome=genome,
        task="qtl",
        task_id=task_id,
        task_spec=task_spec,
        errors=errors,
        form="",
    )
