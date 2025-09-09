#!/usr/bin/env python

import io
import os
import re
import sys
import gzip
import json
import time
import pysam
import statistics
import subprocess
import portalocker
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from app.dirs import tmp_dir, data_dir, tool_dir, prot_db_dir

from modules.vcf_reader import VcfReader
from modules.load_input import load_config

from Bio.Data.IUPACData import protein_letters_3to1 as convertAA

# BLAST options
blastp_evalue = "0.1"
blastp_max_seqs = "5000"

blastp_columns = [
    "qseqid",
    "sseqid",
    "qlen",
]

cdhit_path = "cdhit"
blastp_path = os.path.join(tool_dir, "diamond")
provean_path = os.path.join("..", "..", "..", tool_dir, "provean")
blastdbcmd_path = os.path.join("..", "..", "..", tool_dir,"blast", "bin", "blastdbcmd")

# Add the conversion of "Ter" to "*" to the dictionary
convertAA["Ter"] = "*"

sev_helpstring = {
    "Number": ".",
    "Type": "String",
    "Description": (
        "\"Variant severity information. Format: '"
        + "Allele | "
        + "Transcript | "
        + "GeneId | "
        + "GeneName | "
        + "Region | "
        + "VariantType | "
        + "Variant | "
        + "SEV_score | "
#        + "SEV_stdev | "
        + "SEV_num_seqs | "
        + "SEV_num_clusters | "
        + "SEV_prediction'\""
    ),
}


def run_blastp(genome, prot_file, db_prefix, task_id):

    if any(
        not os.path.exists(os.path.join(prot_db_dir, re.sub("\.phr$", ".acc", file)))
        for file in os.listdir(prot_db_dir) if file.endswith(".phr")):
        print("<b>Preparing DIAMOND database</b>")
        pr = subprocess.run(
            [
                blastp_path,
                "prepdb",
                "--db",
                os.path.join(prot_db_dir, db_prefix),
            ],
            check=True)

    genome_prot_dir = os.path.join(data_dir, genome, "prot_db")
    if not os.path.exists(genome_prot_dir):
        os.mkdir(genome_prot_dir)

    prot_path = os.path.join(data_dir, genome, prot_file)
    tmp_prot_path = os.path.join(tmp_dir, "prot.fasta")
    with gzip.open(prot_path, "rt") as f_in, open(tmp_prot_path, "w") as f_out:
        f_out.write(f_in.read())

    out_path = os.path.join(tmp_dir, "blastp_out.txt")
    print("<b>Running BLAST</b>")
    pr = subprocess.run(
        [
            blastp_path,
            "blastp",
            "--query",
            tmp_prot_path,
            "--db",
            os.path.join(prot_db_dir, db_prefix),
            "--out",
            out_path,
            "--evalue",
            blastp_evalue,
            "--outfmt",
            "6", *blastp_columns,
            "--max-target-seqs",
            blastp_max_seqs,
        ],
        check=True)

    out_gz_path = os.path.join(genome_prot_dir, "blastp_out.txt.gz")
    pysam.tabix_compress(out_path, out_gz_path, force=True)
    os.remove(out_path)
    os.remove(tmp_prot_path)

    pysam.tabix_index(
        out_gz_path, seq_col=0, start_col=2, end_col=2, line_skip=0, force=True
    )

def annotate_sev(genome, markers, task_id):
    print("<b>Annotating variants</b>")
    config = load_config(genome)
    prot_dir = os.path.join(data_dir, genome, "prot_db")
    vcf_path = os.path.join(data_dir, genome, config["vcf"])
    sev_path = os.path.join(data_dir, genome, f"sev.{config['vcf']}")
    vcf_reader = VcfReader(vcf_path)

    if not os.path.exists(sev_path):
        vcf_reader.header.add_header_field(
            name="SEV", field_type="INFO", dictionary=sev_helpstring
        )
        sev_header = (
            vcf_reader.header.__str__().rsplit("\n", 2)[0]
            + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        sev_info = pd.DataFrame()
    else:
        sev_header = VcfReader(sev_path).header.__str__()
        sev_info = pd.read_csv(
            sev_path, sep="\t", comment="#", header=None, index_col=[0, 1, 3, 4]
        )
        sev_info.index = sev_info.index.set_levels(
            sev_info.index.levels[0].astype(str), level=0
        )

    sev_lines = list()
    all_anns = list()

    joined_markers = list()
    for marker in markers:

        if not joined_markers:
            joined_markers.append(list(marker) + [marker[1]])
            continue
        if (
            joined_markers[-1][0] == marker[0]
            and marker[1] - joined_markers[-1][2] < 10000
        ):
            joined_markers[-1][2] = marker[1]
        else:
            joined_markers.append(list(marker) + [marker[1]])

    pep_fa = pysam.FastaFile(os.path.join(data_dir, genome, config["pep"]))

    print("<b>Fetching variant data</b>")
    skip_var = list()
    for marker in joined_markers:
        vcf_reader.set_region(marker[0], marker[1] - 1, marker[2])
        for rec in vcf_reader:
            if [rec.CHROM, rec.POS] not in markers:
                continue
            if (rec.CHROM, rec.POS, rec.REF, ",".join(rec.ALT)) in sev_info.index:
                continue
            if not hasattr(rec, "CSQ"):
                continue

            anns = list()
            for ann in rec.CSQ:
                # Skip variants not affecting the protein sequence
                if not ann["HGVS.p"]:
                    continue
                if "synonymous_variant" in ann["Annotation"]:
                    continue

                #t_id = ann["Feature_ID"]
                #try:
                #    t_seq = pep_fa[t_id]
                #except KeyError:
                #    continue

                # HGVS.p has format p.VAR
                var = ann["HGVS.p"][2:]

                # Convert 3-letter amino acid codes to 1-letter
                var = re.sub("([A-Z][a-z]{2})", lambda x: convertAA[x.group()], var)

                pattern = re.compile(r"""^(
                    [A-Z]\d+[A-Z] |                      # SNP
                    [A-Z]+\d+[A-Z]+ |                    # MNP
                    [A-Z]\d+del |                        # simple deletion
                    [A-Z]\d+_[A-Z]\d+del |               # range deletion
                    [A-Z]\d+_(?:[A-Z]\d+|\*)ins[A-Z]+ |  # insertion
                    [A-Z]\d+delins[A-Z]+ |               # single indel
                    [A-Z]\d+_[A-Z]\d+delins[A-Z]+ |      # range indel
                    [A-Z]\d+dup |                        # single duplication
                    [A-Z]\d+_[A-Z]\d+dup                 # range duplication
                )$""", re.VERBOSE)

                # Stop gains (treat as deletion of residues after the new stop)
#                elif re.match("^[A-Z]\d+\*$", var):
#                    start_var = var[:-1]
#                    start_ref, start_pos = start_var[0], int(start_var[1:])
#                    end_ref = ""
#                    end_pos = len(t_seq) - 1
#                    alt = ""

                # Cannot predict effects of start_gain/start_loss/stop_loss/frameshift
#                elif any(s in var for s in ("fs", "*", "?")):
#                    continue

#                else:
#                    print(f"<b>Unrecognised variant: {var}</b>")
#                    continue

                if re.match(pattern, var):
                    anns.append(
                        [
                            rec.CHROM,
                            rec.POS,
                            rec.REF,
                            ",".join(rec.ALT),
                            ann["Allele"],
                            ann["Feature_ID"],
                            ann["Gene_ID"],
                            ann["Gene_Name"],
                            ann["Annotation"],
                            var,
                        ]
                    )
                else:
                    var_str = f"{rec.CHROM}:{rec.POS} {var}"
                    if  var_str not in skip_var:
                        skip_var.append(var_str)

            if not anns:
                continue
            all_anns += anns

    if skip_var:
        print(
          "Skipping variants that cannot be annotated by PROVEAN:\n" +\
          "\n".join(skip_var))

    print()
    ann_df = pd.DataFrame(
        all_anns,
        columns=[
            "Chrom",
            "Pos",
            "Ref",
            "Alt",
            "Allele",
            "Feature_ID",
            "Gene_ID",
            "Gene_Name",
            "Annotation",
            "Var",
        ],
    )

    blastp_out_fn = os.path.join(prot_dir, "blastp_out.txt.gz")
    blastp_out = pysam.TabixFile(blastp_out_fn)

    for t_id, anns in ann_df.groupby("Feature_ID", sort=False):
        print(f"<b>{t_id}</b>")
        anns = anns.reset_index(drop=True)

        print("Fetching homologous sequences")
        # Get transcript sequence
        try:
            t_seq = pep_fa[t_id]
        except KeyError:
            continue
        record = SeqRecord(Seq(t_seq), id=t_id)
        seq_path = os.path.join(prot_dir, "seq.fasta")
        SeqIO.write(record, seq_path, "fasta")

        # Extract sequences of aligned hits
        try:
            t_blast_out = blastp_out.fetch(t_id)
        except ValueError:
            print(f"No homologs found for {t_id}.")
            continue

        blastp_out_df = pd.DataFrame(
            [line.split("\t") for line in t_blast_out],
            columns=blastp_columns,
        ).convert_dtypes()

        hits_out = os.path.join(prot_dir, "aligned_hits.txt")
        (
            pd.Series(blastp_out_df["sseqid"].unique()).to_csv(
                hits_out, header=False, index=False
            )
        )

        # Save variants
        vars_path = os.path.join(prot_dir, "vars.txt")
        anns["Var"].drop_duplicates().to_csv(vars_path, header=False, index=False)

        old_wd = os.getcwd()
        os.chdir(prot_dir)
        subprocess.run(
            [
                blastdbcmd_path,
                "-db",
                os.path.join("..", "..", "..", "prot_db", config["db_prefix"]),
                "-entry_batch",
                "aligned_hits.txt",
                "-out",
                "aligned_hits.fasta",
                "-outfmt",
                "%f",
            ],
            check=True,
        )

        print("Annotating variants")
        pr = subprocess.run(
            [provean_path,
             "-q", "seq.fasta",
             "-v", "vars.txt",
             "-d", os.path.join("..", "..", "..", prot_db_dir, config["db_prefix" ]),
             "--psiblast", "x",
             "--cdhit", "cdhit",
             "--blastdbcmd", "x",
             "--subject_sequences", "aligned_hits.fasta",
            ],
            capture_output=True,
            text=True)

        provean_out = pr.stdout

        match = re.search(r"#\sNumber of clusters:\s(\d+)", provean_out)
        if match:
            num_clusters = int(match.group(1))
        match = re.search(r"#\sNumber of supporting sequences used:\s(\d+)", provean_out)
        if match:
            num_seqs = int(match.group(1))

        scores = provean_out.rsplit("#", 1)[1].strip()
        scores_df = pd.read_csv(io.StringIO(scores), sep="\s+", index_col=0)

        pos_groups = anns.groupby(["Chrom", "Pos"])
        for marker, pos in pos_groups:
            sev_anns = list()
            for var in pos.to_dict(orient="records"):

                if var["Var"] not in scores_df.index:
                    sev_score = "-"
                else:
                    sev_score = scores_df.at[var["Var"], "SCORE"]
                sev_anns.append(
                    "|".join(
                        map(
                            str,
                            (
                                var["Allele"],
                                var["Feature_ID"],
                                var["Gene_ID"],
                                var["Gene_Name"],
                                "CDS",
                                var["Annotation"],
                                var["Var"],
                                sev_score,
                                num_seqs,
                                num_clusters,
                                (
                                    "-"
                                    if sev_score == "-"
                                    else (
                                        "DELETEROUS_ALT"
                                        if sev_score < -2.5
                                        else (
                                            "TOLERATED"
                                            if sev_score <= 2.5
                                            else "DELETERIOUS_REF"
                                        )
                                    )
                                ),
                            ),
                        )
                    )
                )

            if sev_anns:
                sev_lines.append(
                    (
                        var["Chrom"],
                        var["Pos"],
                        ".",
                        var["Ref"],
                        var["Alt"],
                        ".",
                        ".",
                        "SEV=" + ",".join(sev_anns),
                    )
                )

        os.chdir(old_wd)

    new_sev = pd.DataFrame(sev_lines)
    if not sev_info.empty:
        sev_info = sev_info.reset_index()[[0, 1, 2, 3, 4, 5, 6, 7]]
    if not new_sev.empty:
        sev_info = pd.concat([sev_info, new_sev])
        sev_info[0] = pd.Categorical(
            sev_info[0],
            categories=list(pysam.VariantFile(vcf_path).header.contigs),
            ordered=True,
        )
        sev_info = sev_info.sort_values(by=[0, 1])
        sev_output = io.StringIO()
        sev_info.to_csv(sev_output, sep="\t", header=False, index=False)

        with portalocker.Lock(sev_path, mode="wb", flags=portalocker.LOCK_EX):
            with pysam.BGZFile(sev_path, "wb") as sev_out:
                _ = sev_out.write(sev_header.encode("utf-8"))
                _ = sev_out.write(sev_output.getvalue().encode("utf-8"))
            pysam.tabix_index(sev_path, preset="vcf", force=True)


if __name__ == "__main__":
    genome, markers, task_id = sys.argv[1:4]
    markers = json.loads(markers)
    annotate_sev(genome, markers, task_id)
