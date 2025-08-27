#!/usr/bin/env python

import os
import sys
import gzip
import json
import pysam
import shutil
import gffutils
import numpy as np
import configparser
import pandas as pd
from Bio import SeqIO
from functools import partial

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from app.dirs import data_dir
from modules.safe_col import safe_col
from modules.assoc.assoc_qc import run_qc
from modules.load_input import load_config
from modules.snpeff import snpeff_annotate

from modules.sev import (
    #make_db,
    run_blastp,
)

def genome_setup(genome, params, task_id):
    stats_dir = os.path.join(data_dir, genome, "stats")
    if not os.path.exists(stats_dir):
        os.mkdir(stats_dir)

    if os.path.exists(os.path.join(data_dir, genome, "stats", "init")):
        raise OSError("<b>The selected dataset is already set up.</b>")

    if params["setup_from"]:
        genome_from = params["setup_from"]
        config_from = load_config(genome_from)

        # Process VCF file
        file_path = os.path.join(data_dir, genome, params["vcf"])
        if not file_path.endswith(".gz"):
            if os.path.islink(file_path):
                print(f"Copying vcf file")
                link_path = os.readlink(file_path)
                os.unlink(file_path)
                shutil.copy2(link_path, file_path)

            print(f"Compressing vcf file")
            pysam.tabix_compress(file_path, file_path + ".gz")
            os.remove(file_path)
            params["vcf"] += ".gz"

        file_path = os.path.join(data_dir, genome, params["vcf"])
        if not os.path.exists(file_path + ".tbi"):
            print(f"Indexing vcf file")
            pysam.tabix_index(
                file_path,
                preset="vcf",
                force=True,
            )

        # Link other files
        for key in (
            "seq", "ann", "exp", "run", "info", "bam",
            "cds", "pep", "db_prefix"):
            params[key] = config_from[key]
        params["concat"] = json.loads(config_from["concat"])
        for key in ("seq", "cds", "pep"):
            for ext in ("", ".fai", ".gzi"):
                if os.path.exists(os.path.join(data_dir, genome, config_from[key]+ext)):
                    continue
                os.link(
                    os.path.join(data_dir, genome_from, config_from[key]+ext),
                    os.path.join(data_dir, genome, config_from[key]+ext)
                )
        for key in ("ann", "exp"):
            for ext in ("", ".tbi"):
                if os.path.exists(os.path.join(data_dir, genome, config_from[key]+ext)):
                    continue
                os.link(
                    os.path.join(data_dir, genome_from, config_from[key]+ext),
                    os.path.join(data_dir, genome, config_from[key]+ext)
                )
        for bam in json.loads(config_from["bam"]):
            for ext in ("", ".bai"):
                if os.path.exists(os.path.join(data_dir, genome, config_from[key]+ext)):
                    continue
                os.link(
                    os.path.join(data_dir, genome_from, config_from[key]+ext),
                    os.path.join(data_dir, genome, config_from[key]+ext)
                )
        for key in ("run", "info"):
            if os.path.exists(os.path.join(data_dir, genome, config_from[key])):
                continue
            os.link(
                os.path.join(data_dir, genome_from, config_from[key]),
                os.path.join(data_dir, genome, config_from[key])
            )
        if not os.path.exists(os.path.join(data_dir, genome, "stats", "gff.db")):
            os.link(
                os.path.join(data_dir, genome_from, "stats", "gff.db"),
                os.path.join(data_dir, genome, "stats", "gff.db")
            )

        try:
            os.symlink(
                os.path.join(data_dir, genome_from, "prot_db"),
                os.path.join(data_dir, genome, "prot_db")
            )
        except FileExistsError:
            pass

        db = gffutils.FeatureDB(os.path.join(data_dir, genome, "stats", "gff.db"))

    else:

###############################################################################
### Add a dummy pos column to expression file
###############################################################################

        if params["exp"]:
            file_path = os.path.join(data_dir, genome, params["exp"])
            if os.path.islink(file_path):
                print(f"Copying exp file")
                link_path = os.readlink(file_path)
                os.unlink(file_path)
                shutil.copy2(link_path, file_path)

            is_gz = file_path.endswith(".gz")
            if not is_gz:
                params["exp"] += ".gz"

            out_path = file_path[:-3] if is_gz else file_path
            exp = pd.read_csv(
                file_path,
                sep="\t",
                index_col=0,
            )
            exp["Pos"] = 1
            exp[["Pos"] + exp.columns.tolist()[:-1]].to_csv(
                out_path,
                sep="\t",
            )
            pysam.tabix_compress(out_path, out_path + ".gz", force=True)
            file_path = out_path + ".gz"
            
            os.remove(out_path)
            pysam.tabix_index(
                file_path,
                seq_col=0,
                start_col=1,
                end_col=1,
                line_skip=1,
                force=True,
            )

###############################################################################
### BGZip and tabix files if necessary
###############################################################################

        for file in ("vcf", "ann", "seq", "pep", "cds"):
            file_path = os.path.join(data_dir, genome, params[file])
            if params[file] and not file_path.endswith(".gz"):
                if os.path.islink(file_path):
                    print(f"Copying {file} file")
                    link_path = os.readlink(file_path)
                    os.unlink(file_path)
                    shutil.copy2(link_path, file_path)

                print(f"Compressing {file} file")
                pysam.tabix_compress(file_path, file_path + ".gz")
                os.remove(file_path)
                params[file] += ".gz"

            file_path = os.path.join(data_dir, genome, params[file])
            if params[file] and file in ("seq", "pep", "cds"):
                if not os.path.exists(file_path + ".fai") or not os.path.exists(
                    file_path + ".gzi"
                ):
                    print(f"Indexing {file} file")
                    _ = pysam.FastaFile(file_path)
            else:
                idx_path = file_path + (".tbi" if not params[file].endswith(".bcf.gz") else ".csi")
                if params[file] and not os.path.exists(idx_path):
                    print(f"Indexing {file} file")
                    pysam.tabix_index(
                        file_path,
                        preset=("vcf" if file == "vcf" else "gff"),
                        force=True,
                    )
        for file in params["bam"]:
            file_path = os.path.join(data_dir, genome, file)
            if not os.path.exists(file_path + ".bai"):
                print("Indexing a bam file")
                pysam.index(file)

        if params["exp"] and not params["run"]:
            file_path = os.path.join(data_dir, genome, params["exp"])
            exp_runs = pd.read_csv(file_path, sep="\t", nrows=0)
            (
                pd.DataFrame({"Run_ID": exp_runs.columns.tolist()[4:]}).to_csv(
                    os.path.join(data_dir, genome, "run_info.txt"), sep="\t", index=False
                )
            )
            params["run"] = "run_info.txt"

###############################################################################
### Set up GFFUtils database
###############################################################################

        ann_path = os.path.join(data_dir, genome, params["ann"])
        db_path = os.path.join(data_dir, genome, "stats", "gff.db")
        if not os.path.exists(db_path):
            print("<b>Creating an annotation database</b>")
            gffutils.create_db(ann_path, dbfn=db_path, merge_strategy="merge")
        #db = gffutils.FeatureDB(db_path)

###############################################################################
### Add gene attributes to GFFUtils database
###############################################################################

        if params["info"]:
            print("<b>Updating database information</b>")
            annotate_db(genome, db, params["info"], task_id)
        # if params["enr"]:
        #    print("Calculating gene enrichment")
        #    num_genes = db.count_features_of_type("gene")
        #    calc_enr(
        #      genome,
        #      params["info"],
        #      params["enr"],
        #      num_genes,
        #      task_id)

###############################################################################
### Annotate variants and setup BLAST / severity prediction
###############################################################################

        if params["sev"]:
            # Some commands cannot handle spaces in paths
            # (which are allowed in genome names)
            old_wd = os.getcwd()
            os.chdir(os.path.join(data_dir, genome))
            # if params["db_fasta"]:
            #    make_db(
            #      task_id,
            #      params["db_fasta"])
            # params["db_prefix"] = params["db_fasta"]
            if params["db_prefix"]:
                run_blastp(task_id, params["pep"], params["db_prefix"])
                # process_blastp()
            os.chdir(old_wd)

#    else:
#        if params["pep"]:
#            prot_dir = os.path.join(data_dir, genome, "prot_db", "prot")
#            if not os.path.exists(prot_dir):
#                os.makedirs(prot_dir, exist_ok=True)
#            prot_file = os.path.join(data_dir, genome, params["pep"])

#            _open = partial(gzip.open, mode="rt") if prot_file.endswith(".gz") else open
#            with _open(prot_file) as f:
#                for prot in SeqIO.parse(f, "fasta"):
#                    prot_path = os.path.join(prot_dir, f"{prot.id}.fasta")
#                    SeqIO.write([prot], prot_path, "fasta")

#        if params["cds"]:
#            cds_dir = os.path.join(data_dir, genome, "prot_db", "cds")
#            if not os.path.exists(cds_dir):
#                os.makedirs(cds_dir, exist_ok=True)
#            cds_file = os.path.join(data_dir, genome, params["cds"])

#            _open = partial(gzip.open, mode="rt") if cds_file.endswith(".gz") else open
#            with _open(cds_file) as f:
#                for cds in SeqIO.parse(f, "fasta"):
#                    cds_path = os.path.join(cds_dir, f"{cds.id}.fasta")
#                    SeqIO.write([cds], cds_path, "fasta")

###############################################################################
### Get genome and sample data
###############################################################################

    print("<b>Parsing VCF file</b>")
    vcf_path = os.path.join(data_dir, genome, params["vcf"])

    # Create a genome data file if none exists
    gdata_path = os.path.join(data_dir, genome, "stats", ".gdata")
    #    if not os.path.exists(gdata_path):
    gdata = dict()

    # Get chromosome names, sizes and offsets
    idx_path = os.path.join(data_dir, genome, params["seq"] + ".fai")
    idx_df = pd.read_csv(idx_path, sep="\t", header=None)
    chroms = idx_df.iloc[:, 0].astype(str).tolist()
    sizes = idx_df.iloc[:, 1].tolist()

    if params["concat"]:
        chroms = [c for c in chroms if c not in params["concat"]] + [
            c for c in chroms if c in params["concat"]
        ]

    offsets = [sum(sizes[:i]) for i in range(len(sizes))]
    gdata["chr_names"] = chroms
    gdata["chr_sizes"] = sizes
    gdata["chr_offsets"] = offsets

    # Get first gene ID in the file
    #for rec in db.region(region=str(chroms[0]), featuretype="gene"):
    #    gdata["id_fmt"] = rec.id
    #    break
    gdata["id_fmt"] = "PGSC0003DMG400015133"

    # Check that samples match between VCF and SEQ files
    chroms = idx_df.iloc[:, 0].astype(str).tolist()
    vcf = pysam.VariantFile(vcf_path)
    vcf_chroms = list(vcf.header.contigs)
    vcf_sizes = [vcf.header.contigs[c].length for c in vcf_chroms]

    if not all(vcf_chroms) or not all(vcf_sizes):
        print(
            "<b>Cannot verify that chromosome names and sizes match "
            "between vcf and seq files.</b>"
        )

    elif vcf_chroms != chroms or vcf_sizes != sizes:
        print(
            "<b>Chromosomes do not match in size or names between "
            "vcf and seq files.</b>"
        )
        raise OSError

    if params["snpeff"]:
        print("<b>Running VCF annotation</b>")
        snpeff_annotate(
            genome,
            params["vcf"],
            params["ann"],
            params["seq"],
            params["cds"],
            params["pep"],
            task_id,
        )
        params["vcf"] = f"ann.{params['vcf']}"

    print("<b>Parsing phenotypes</b>")
    # Save sample list, combine with phenotyped samples
    samples = list(vcf.header.samples)
    if not params["phe"]:
        params["phe"] = "phe.txt"
        phe_path = os.path.join(data_dir, genome, params["phe"])
        with open(phe_path, "w") as phe_out:
            phe_out.write("\n".join(["Acc"] + samples))
    else:
        phe_path = os.path.join(data_dir, genome, params["phe"])
        if os.path.islink(phe_path):
            link_path = os.readlink(phe_path)
            os.unlink(phe_path)
            shutil.copy2(link_path, phe_path)
        phe_df = pd.read_csv(phe_path, sep="\t", header=0, index_col=0)
        # Make trait names filename-friendly
        safe_cols = [safe_col(col) for col in phe_df.columns]
        phe_df = phe_df.set_axis(safe_cols, axis=1)

        # De-duplicate trait names
        cols = pd.Series(phe_df.columns)
        for dup in phe_df.columns[phe_df.columns.duplicated(keep=False)]:
            cols[phe_df.columns.get_loc(dup)] = [
                dup + "." + str(d_idx) if d_idx != 0 else dup
                for d_idx in range(phe_df.columns.get_loc(dup).sum())
            ]
        phe_df.columns = cols

        for sample in phe_df.index:
            if sample not in samples:
                samples.append(sample)
        miss_df = pd.DataFrame(
            np.nan,
            index=[s for s in samples if s not in phe_df.index],
            columns=phe_df.columns,
        )

        (
            pd.concat([phe_df, miss_df])
            .loc[samples]
            .rename_axis("Acc")
            .to_csv(phe_path, sep="\t")
        )

    # Rename samples
    if not params["samp"]:
        params["samp"] = "samp.txt"
        samp_path = os.path.join(data_dir, genome, params["samp"])
        (
            pd.DataFrame(samples, index=samples, columns=["Name"])
            .rename_axis("Acc")
            .to_csv(samp_path, sep="\t")
        )
    else:
        samp_path = os.path.join(data_dir, genome, params["samp"])
        if os.path.islink(samp_path):
            link_path = os.readlink(samp_path)
            os.unlink(samp_path)
            shutil.copy2(link_path, samp_path)
        samp_df = pd.read_csv(samp_path, sep="\t", header=0, index_col=0)
        samp_df.index = samp_df.index.map(str)

        miss_samples = [s for s in samples if s not in samp_df.index]
        miss_df = pd.DataFrame(
            miss_samples, index=miss_samples, columns=[samp_df.columns[0]]
        )
        (
            pd.concat([samp_df, miss_df])
            .loc[samples]
            .rename_axis("Acc")
            .to_csv(samp_path, sep="\t")
        )

    # Initialise empty list for groups
    gdata["group"] = dict()

    # Save genome data to JSON file
    with open(gdata_path, "w") as gdata_out:
        json.dump(gdata, gdata_out)

    # Save parameters as config.ini
    config_params = {k: v for k, v in params.items()}
    _ = [config_params.pop(v) for v in ("name", "csrf_token")]
    config_params["bam"] = json.dumps(config_params["bam"])
    # config_params["enr"] = json.dumps(config_params["enr"])
    config_params["concat"] = json.dumps(config_params["concat"])
    config_params = {k: str(v) for k, v in config_params.items()}
    config_params["rm_markers"] = True

    config = configparser.ConfigParser()
    config["setup"] = config_params
    config_path = os.path.join(data_dir, genome, "stats", "config.ini")
    with open(config_path, "w") as config_out:
        config.write(config_out)

    # Run raw marker set QC
    assoc_dir = os.path.join(data_dir, genome, "assoc")
    if not os.path.exists(assoc_dir):
        os.mkdir(assoc_dir)
    cof_dir = os.path.join(data_dir, genome, "covariates")
    if not os.path.exists(cof_dir):
        os.mkdir(cof_dir)

    if not params["run_ext"]:
        run_qc(genome, task_id)

        init_path = os.path.join(stats_dir, "init")
        open(init_path, "w").close()


def annotate_db(genome, db, info, task_id):
    info_path = os.path.join(data_dir, genome, info)
    info_df = pd.read_csv(info_path, sep="\t", header=0, index_col=0)
    info_df = info_df.fillna("")
    info_df = info_df.astype(str)
    cols = info_df.columns

    new_entries = list()
    for tr in db.features_of_type("mRNA"):
        # GFFutils splits attributes by comma
        # Rejoin them for uniformity
        for attr, val in tr.attributes.items():
            tr.attributes[attr] = ",".join(val)
        # Save new attributes
        if tr.id in info_df.index:
            data = info_df.loc[[tr.id]].reset_index(drop=True).iloc[0]
            for col in cols:
                val = data[col]
                tr.attributes[col] = val  # if not pd.isna(val) else ""
        else:
            for col in cols:
                tr.attributes[col] = ""
        new_entries.append(tr)
    db = db.update(new_entries, merge_strategy="replace", make_backup=False)


# def calc_enr(
#  genome,
#  info,
#  enr_cols,
#  num_genes,
#  task_id):
#    info_path = os.path.join(
#      data_dir,
#      genome,
#      info)
#    info_df = pd.read_csv(
#      info_path,
#      sep="\t",
#      header=0,
#      index_col=0)
#
#    enr = dict()
#    for col in enr_cols:
#        attrs = []
#        for row_attrs in info_df[col].dropna():
#            attrs += list(
#              set(
#                [attr.strip() for attr
#                 in row_attrs.strip().split(";") if attr]))
#        enr[col] = Counter(attrs)
#
#    enr_path = os.path.join(
#      data_dir,
#      genome,
#      ".enr")
#    with open(enr_path, "w") as enr_out:
#        json.dump([enr, num_genes], enr_out)

if __name__ == "__main__":
    genome, params, task_id = sys.argv[1:4]
    params = json.loads(params)
    genome_setup(genome, params, task_id)
