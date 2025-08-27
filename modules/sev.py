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
from Bio import SeqIO, Seq, SeqRecord
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from app.dirs import data_dir, tool_dir

from modules.vcf_reader import VcfReader
from modules.load_input import load_config

from Bio.Data.IUPACData import protein_letters_3to1 as convertAA

# BLAST options
blastp_evalue = "0.1"
blastp_max_seqs = "5000"

blastp_columns = [
    "qseqid",
    "sacc",
    "sallseqid",
    "qstart",
    "qend",
    "evalue",
    "qlen",
    "sstart",
    "send",
    "slen",
    "pident",
    # "salltitles",
    # "length",
    # "nindent",
    # "score",
]

cdhit_path = "cdhit"
blastp_path = os.path.join(tool_dir,"blast", "bin", "blastp")
blastdbcmd_path = os.path.join("..", "..", "..", tool_dir,"blast", "bin", "blastdbcmd")

seq_identity = 0.8
cluster_to_most_similar = 1


# Add the conversion of "Ter" to "*" to the dictionary
convertAA["Ter"] = "*"

aln_gapopen = -10
aln_gapextend = -1
aln_endscore = 0

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
        + "SEV_stdev | "
        + "SEV_num_seqs | "
        + "SEV_num_clusters | "
        + "SEV_prediction'\""
    ),
}


def run_blastp(task_id, prot_file, db_prefix):
    # Split input sequences file into single-sequence FASTA files
    # Initialize a progress counter based on number of proteins finished
    # print("<b>Processing input proteins</b>")
    # split_dir = os.path.join("prot_db", "prot")
    # if not os.path.exists(split_dir):
    #    os.mkdir(split_dir)

    with gzip.open(prot_file, "rt") as f_in, open("prot.fasta", "w") as f_out:
        f_out.write(f_in.read())

    prot_order = dict()
    for idx, prot in enumerate(SeqIO.parse("prot.fasta", "fasta"), start=1):
        #    prot_path = os.path.join(
        #      split_dir,
        #      f"{prot.id}.fasta")
        #    SeqIO.write(
        #      [prot],
        #      prot_path,
        #      "fasta")
        prot_order[prot.id] = idx
    prot_cnt = len(prot_order)

    out_path = os.path.join("prot_db", "blastp_out.txt")
    print("<b>Running BLAST</b>")
    pr = subprocess.Popen(
        [
            blastp_path,
            "-query",
            "prot.fasta",
            "-db",
            os.path.join("prot_db", db_prefix),
            "-out",
            out_path,
            "-evalue",
            blastp_evalue,
            "-num_threads",
            str(os.cpu_count()),
            "-outfmt",
            "6 " + " ".join(blastp_columns),
            "-max_target_seqs",
            blastp_max_seqs,
        ],
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    while pr.poll() is None:
        tail = subprocess.Popen(
            ["tail", "-n", "1", out_path],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        stdout, stderr = tail.communicate()
        if not stderr.decode():
            last_prot = stdout.decode().split("\t", 1)[0]
            if last_prot:
                pct = int(prot_order[last_prot] / prot_cnt * 100)
                print(f"\rProgress: {pct}%", end="")
            else:
                print("\rProgress: 0%", end="")
        time.sleep(5)

    stdout, stderr = pr.communicate()
    if pr.returncode == 0:
        print(f"\rProgress: 100%")
    else:
        raise OSError("BLAST failed.")

    # Index blastp file
    pysam.tabix_compress(out_path, out_path + ".gz", force=True)
    os.remove(out_path)
    os.remove("prot.fasta")

    pysam.tabix_index(
        out_path + ".gz", seq_col=0, start_col=6, end_col=6, line_skip=0, force=True
    )


# def process_blastp():
#    split_dir = os.path.join("prot_db", "blastp_out")
#    if not os.path.exists(split_dir):
#        os.mkdir(split_dir)
#
#    blastp_reader = pd.read_csv(
#      os.path.join("prot_db", "blastp_out.txt"),
#      sep="\t",
#      header=None,
#      chunksize=100000)
#    prev_data = None
#
#    for chunk in blastp_reader:
#        chunk = chunk.set_axis(
#          blastp_columns,
#          axis=1)
#        if prev_data != None:
#            chunk = pd.concat([prev_data, chunk])
#
#        for prot_id, hits in chunk.groupby("qseqid", sort=False):
#            # Save hits, evalues and aligned regions
#            out_file = os.path.join(
#              split_dir,
#              f"{prot_id}.blastp_out.txt")
#            (hits[
#               ["sacc",
#                "sallseqid",
#                "qstart",
#                "qend",
#                "evalue"]]
#              .to_csv(
#                out_file,
#                index=False,
#                header=False,
#                sep="\t"))
#        prev_data = hits


def annotate_sev(genome, markers, task_id):
    # print("<b>Annotating variants</b>")
    config = load_config(genome)
    vcf_path = os.path.join(data_dir, genome, config["vcf"])
    # blast_db_root = os.path.join(
    #  data_dir,
    #  genome,
    #  "prot_db",
    #  config["db_fasta"])
    prot_dir = os.path.join(data_dir, genome, "prot_db")
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
            and marker[1] - joined_markers[-1][2] < 1000
        ):
            joined_markers[-1][2] = marker[1]
        else:
            joined_markers.append(list(marker) + [marker[1]])

    pep_fa = pysam.FastaFile(os.path.join(data_dir, genome, config["pep"]))

    print("<b>Fetching variant data</b>")
    for marker in joined_markers:
        marker = tuple(marker)
        vcf_reader.set_region(marker[0], marker[1] - 1, marker[2])
        for rec in vcf_reader:
            if (rec.CHROM, rec.POS, rec.REF, ",".join(rec.ALT)) in sev_info.index:
                continue
            try:
                _ = rec.CSQ
            except AttributeError:
                continue

            anns = list()
            for ann in rec.CSQ:
                # Skip variants not affecting the protein sequence
                if not ann["HGVS.p"]:
                    continue
                if "synonymous_variant" in ann["Annotation"]:
                    continue

                t_id = ann["Feature_ID"]
                try:
                    t_seq = pep_fa[t_id]
                except KeyError:
                    continue

                # HGVS.p has format p.VAR
                var = ann["HGVS.p"][2:]

                # Convert 3-letter amino acid codes to 1-letter
                var = re.sub("([A-Z][a-z]{2})", lambda x: convertAA[x.group()], var)

                # SNPs
                if re.match("^[A-Z]\d+[A-Z]$", var):
                    ref, pos, alt = var[0], var[1:-1], var[-1]
                    start_ref = ref
                    end_ref = ref
                    start_pos = int(pos)
                    end_pos = start_pos

                # MNPs
                elif re.match("^[A-Z]+\d+[A-Z]+$", var):
                    ref = re.search("^[A-Z]+", var).group()
                    pos = re.search("\d+", var).group()
                    alt = re.search("[A-Z]+$", var).group()
                    start_ref = ref[0]
                    end_ref = ref[-1]
                    start_pos = int(pos)
                    end_pos = start_pos + (len(ref) - 1)

                # Deletions
                elif re.match("^[A-Z]\d+del$", var):
                    start_var = var[:-3]
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref, end_pos = start_ref, start_pos
                    alt = ""

                elif re.match("^[A-Z]\d+\_[A-Z]\d+del$", var):
                    start_var, end_var = var[:-3].split("_")
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref, end_pos = end_var[0], int(end_var[1:])
                    alt = ""

                # Insertions (only if no stop codon was inserted)
                elif re.match("^[A-Z]\d+\_([A-Z]\d+|\*)ins[A-Z]+$", var):
                    start_var, end_var, alt = var.replace("ins", "_").split("_")
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref = end_var[0]
                    end_pos = int(end_var[1:]) if end_var[1:] else len(t_seq)

                # Del-insertions (only if no stop codon was deleted or inserted)
                elif re.match("^[A-Z]\d+delins[A-Z]+$", var):
                    start_var, alt = var.replace("delins", "_").split("_")
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref, end_pos = start_ref, start_pos

                elif re.match("^[A-Z]\d+\_[A-Z]\d+delins[A-Z]+$", var):
                    start_var, end_var, alt = var.replace("delins", "_").split("_")
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref, end_pos = end_var[0], int(end_var[1:])

                # Duplications
                elif re.match("^[A-Z]\d+dup$", var):
                    start_var = var[:-3]
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref, end_pos = start_ref, start_pos
                    alt = str(t_seq[start_pos - 1 : end_pos]) * 2

                elif re.match("^[A-Z]\d+\_[A-Z]\d+dup$", var):
                    start_var, end_var = var[:-3].split("_")
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref, end_pos = end_var[0], int(end_var[1:])
                    alt = str(t_seq[start_pos - 1 : end_pos]) * 2

                # Stop gains (treat as deletion of residues after the new stop)
                elif re.match("^[A-Z]\d+\*$", var):
                    start_var = var[:-1]
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref = ""
                    end_pos = len(t_seq) - 1
                    alt = ""

                # Cannot predict effects of start_gain/start_loss/stop_loss/frameshift
                elif any(s in var for s in ("fs", "*", "?")):
                    # print("Skipping", var)
                    continue

                else:
                    print(f"<b>Unrecognised variant: {var}</b>")
                    continue

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
                        start_ref,
                        start_pos,
                        end_ref,
                        end_pos,
                        alt,
                    ]
                )

            if not anns:
                continue
            all_anns += anns

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
            "Start_ref",
            "Start_pos",
            "End_ref",
            "End_pos",
            "Subs",
        ],
    )
    ann_df["Supporting_set"] = ""

    blastp_out_fn = os.path.join(prot_dir, "blastp_out.txt.gz")
    blastp_out = pysam.TabixFile(blastp_out_fn)

    aligner = PairwiseAligner(
        mode="global",
        internal_open_gap_score=aln_gapopen,
        internal_extend_gap_score=aln_gapextend,
        end_gap_score=aln_endscore,
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
    )

    for t_id, anns in ann_df.groupby("Feature_ID", sort=False):
        print(f"<b>{t_id}</b>")
        anns = anns.reset_index(drop=True)
        full_set = set()

        print("Fetching homologous sequences")
        # Get transcript sequence
        try:
            t_seq = pep_fa[t_id]
        except KeyError:
            continue

        blastp_out_df = pd.DataFrame(
            [line.split("\t") for line in blastp_out.fetch(t_id)],
            columns=blastp_columns,
        ).convert_dtypes()
        blastp_out_df[["qstart", "qend", "evalue"]] = blastp_out_df[
            ["qstart", "qend", "evalue"]
        ].apply(pd.to_numeric, errors="coerce")

        # Extract sequences of aligned hits
        hits_out = os.path.join(prot_dir, "aligned_hits.txt")
        (
            pd.Series(blastp_out_df["sallseqid"].unique()).to_csv(
                hits_out, header=False, index=False
            )
        )

        old_wd = os.getcwd()
        os.chdir(os.path.join(data_dir, genome, "prot_db"))
        subprocess.run(
            [
                blastdbcmd_path,
                "-db",
                config["db_prefix"],
                "-entry_batch",
                "aligned_hits.txt",
                "-out",
                "aligned_hits.fasta",
                "-outfmt",
                "%f",
            ],
            check=True,
        )

        # Cluster by similarity
        subprocess.run(
            [
                cdhit_path,
                "-i",
                "aligned_hits.fasta",
                "-o",
                "clust",
                "-c",
                str(seq_identity),
                "-g",
                str(cluster_to_most_similar),
            ],
            stdout=subprocess.DEVNULL,
            check=True,
        )
        with open("clust.clstr") as clust_in:
            clust_raw = clust_in.read().strip().split("\n>")

        hits = {rec.name: rec for rec in SeqIO.parse("aligned_hits.fasta", "fasta")}

        clust_seqs = list()
        clust_evalues = list()
        for cluster in clust_raw:
            # Format: >Cluster n
            #         <cluster sequences>
            #         >Cluster n+1
            clust_lines = [line for line in cluster.strip().split("\n")[1:]]
            # Format: <whitespace> >seq_id...at x% <identity>
            seq_names = [line.split(">")[1].split("...")[0] for line in clust_lines]
            clust_seqs.append(seq_names)
            mean_evalue = blastp_out_df[blastp_out_df["sacc"].isin(seq_names)][
                "evalue"
            ].mean()
            clust_evalues.append((seq_names, mean_evalue))
        clusters = [cluster[0] for cluster in sorted(clust_evalues, key=lambda x: x[1])]

        os.chdir(old_wd)

        for idx in range(anns.shape[0]):
            start_pos, end_pos = anns.loc[idx, ["Start_pos", "End_pos"]]
            supp_set = list()

            matching_seqs = set(
                blastp_out_df[
                    (blastp_out_df["qstart"] <= start_pos)
                    & (blastp_out_df["qend"] >= end_pos)
                ]["sacc"]
            )
            for cluster in clusters:
                filt_cluster = [seq for seq in cluster if seq in matching_seqs]
                if filt_cluster:
                    supp_set.append(filt_cluster)
                    if len(supp_set) == 45:
                        break
            anns.at[idx, "Supporting_set"] = supp_set
            full_set.update(sum(supp_set, []))

        if len(full_set) == 0:
            continue
        full_set = [hits[seq_name] for seq_name in full_set]

        wt_scores = list()
        try:
            for seq in full_set:
                alignments = aligner.align(t_seq, seq.seq)
                wt_scores.append(statistics.mean([a.score for a in alignments]))
        except Exception:
            print(f"<b>Alignment for ref-{t_id} failed.</b>")
            continue

        pos_groups = anns.groupby(["Chrom", "Pos"])
        num_markers = len(pos_groups)
        for idx, (marker, pos) in enumerate(pos_groups, start=1):
            print(f"\rVariant {idx} of {num_markers}", end="")
            sev_anns = list()
            for var in pos.to_dict(orient="records"):
                if len(var["Supporting_set"]) == 0:
                    continue

                mutated_seq = Seq.MutableSeq(t_seq)
                if mutated_seq[var["Start_pos"] - 1] != var["Start_ref"] or (
                    var["End_ref"] and mutated_seq[var["End_pos"] - 1] != var["End_ref"]
                ):
                    print()
                    print(
                        "<b>Mismatch between variant and sequence "
                        + f"in {t_id}, {':'.join(list(map(str, marker)))}, {var['Var']}</b>"
                    )
                    continue

                mutated_seq[var["Start_pos"] - 1 : var["End_pos"]] = var["Subs"]

                # Get alignment for the mutated sequence
                mut_scores = list()
                try:
                    for seq in full_set:
                        alignments = aligner.align(mutated_seq, seq.seq)
                        mut_scores.append(statistics.mean([a.score for a in alignments]))
                except Exception:
                    print()
                    print(f"<b>Alignment for {t_id}-{var['Subs']} failed.</b>")
                    continue

                if len(wt_scores) != len(mut_scores):
                    print()
                    print(
                        f"<b>Alignment output mismatched in length: {t_id}-{var['Subs']}.</b>"
                    )
                    continue

                score_df = pd.DataFrame(
                  {"wt": wt_scores, "mut": mut_scores}).set_axis([seq.id for seq in full_set])
                score_df["delta"] = score_df["mut"] - score_df["wt"]

                avg_deltas = list()
                for cluster in var["Supporting_set"]:
                    deltas = score_df[score_df.index.isin(cluster)]["delta"]
                    avg_deltas.append(deltas.mean())
                    sev_score = (
                        round(statistics.mean(avg_deltas), 4)
                        if len(avg_deltas) > 0
                        else "-"
                    )
                    sev_stdev = (
                        round(statistics.stdev(avg_deltas), 4)
                        if len(avg_deltas) > 1
                        else "-"
                    )

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
                                sev_stdev,
                                len(sum(var["Supporting_set"], [])),
                                len(var["Supporting_set"]),
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
        print()

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
