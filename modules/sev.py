#!/usr/bin/env python

import io
import os
import re
import time
import pysam
import statistics
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO, Seq, SeqRecord

from flask_app import (
  get_lock,
  master_lock)

from app.dirs import data_dir

from app.task_worker import progress

from modules.assoc import popen
from modules.vcf_reader import VcfReader
from modules.load_input import load_config

from Bio.Data.IUPACData import protein_letters_3to1 as convertAA

# BLAST options
blastp_evalue = "0.1"
blastp_max_seqs = "5000"

blastp_columns = [
  "qseqid",
  "sseqid",
  "sallseqid",
  "qstart",
  "qend",
  "evalue",
  "qlen",
  "sstart",
  "send",
  "slen",
  "pident",
  #"salltitles",
  #"length",
  #"nindent",
  #"score",
]

cdhit_path = "cdhit"
blastdbcmd_path = "blastdbcmd"

seq_identity = 0.8
cluster_to_most_similar = 1


# Add the conversion of "Ter" to "*" to the dictionary
convertAA["Ter"] = "*"

# Needle options
needle_path = "needle"

needle_gapopen = "10"
needle_gapextend = "1"

sev_helpstring = {
  "Number": ".",
  "Type": "String",
  "Description": (
    '"Variant severity information. Format: \''
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
    + 'SEV_prediction\'"')}

def make_db(task_id, db_fasta):
    progress(task_id, "Making a BLAST database")
    popen(
      task_id,
      ["makeblastdb",
       "-in", os.path.join("prot_db", db_fasta),
       "-out", os.path.join("prot_db", db_fasta),
       "-dbtype", "prot",
       "-parse_seqids"])


def run_blastp(task_id, prot_file, db_prefix):
    # Split input sequences file into single-sequence FASTA files
    # Initialize a progress counter based on number of proteins finished
    progress(task_id, "Processing input proteins")
    split_dir = os.path.join("prot_db", "prot")
    if not os.path.exists(split_dir):
        os.mkdir(split_dir)

    prot_order = dict()
    for idx, prot in enumerate(
      SeqIO.parse(
        prot_file,
        "fasta"),
      start=1):
        prot_path = os.path.join(
          split_dir,
          f"{prot.id}.fasta")
        SeqIO.write(
          [prot],
          prot_path,
          "fasta")
        prot_order[prot.id] = idx
    prot_cnt = len(prot_order)

    progress(task_id, "Running BLAST")
    pr = subprocess.Popen(
      ["blastp",
       "-query", prot_file,
       "-db", os.path.join("prot_db", db_prefix),
       "-out", os.path.join("prot_db", "blastp_out.txt"),
       "-evalue", blastp_evalue,
       "-num_threads", str(os.cpu_count()),
       "-outfmt", "6 " + " ".join(blastp_columns),
       "-max_target_seqs", blastp_max_seqs],
      stderr=subprocess.PIPE,
      stdout=subprocess.PIPE)

    while pr.poll() is None:
        tail = subprocess.Popen(
          ["tail", 
           "-n", "1",
           os.path.join("prot_db", "blastp_out.txt")],
          stderr=subprocess.PIPE,
          stdout=subprocess.PIPE)
        stdout, stderr = tail.communicate()
        if not stderr.decode():
            last_prot = stdout.decode().split("\t", 1)[0]
            pct = int(prot_order[last_prot] / prot_cnt * 100)
            progress(
              task_id,
              f"\rProgress: {pct}%",
              bold=False,
              newline=False)
        time.sleep(1)

    stdout, stderr = pr.communicate()
    if pr.returncode == 0:
        progress(task_id, f"\rProgress: 100%", bold=False)
    else:
        raise OSError("Process failed.")

def process_blastp():
    split_dir = os.path.join("prot_db", "blastp_out")
    if not os.path.exists(split_dir):
        os.mkdir(split_dir)

    blastp_reader = pd.read_csv(
      os.path.join("prot_db", "blastp_out.txt"),
      sep="\t",
      header=None,
      chunksize=100000)
    prev_data = None

    for chunk in blastp_reader:
        chunk = chunk.set_axis(
          blastp_columns,
          axis=1)
        if prev_data != None:
            chunk = pd.concat([prev_data, chunk])

        for prot_id, hits in chunk.groupby("qseqid", sort=False):
            # Save hits, evalues and aligned regions
            out_file = os.path.join(
              split_dir,
              f"{prot_id}.blastp_out.txt")
            (hits[
               ["sseqid",
                "sallseqid",
                "qstart",
                "qend",
                "evalue"]]
              .to_csv(
                out_file,
                index=False,
                header=False,
                sep="\t"))
        prev_data = hits


def annotate_sev(genome, markers, task_id):
    progress(task_id, "Annotating variants")
    config = load_config(genome)
    vcf_path = os.path.join(
      data_dir,
      genome,
      config["vcf"])
    blast_db_root = os.path.join(
      data_dir,
      genome,
      "prot_db",
      config["db_prefix"])
    prot_dir = os.path.join(
      data_dir,
      genome,
      "prot_db")
    sev_path = os.path.join(
      data_dir,
      genome,
      f"sev.{config['vcf']}")
    vcf_reader = VcfReader(vcf_path)

    if not os.path.exists(sev_path):
        vcf_reader.header.add_header_field(
          name="SEV",
          field_type="INFO",
          dictionary=sev_helpstring)
        sev_header = (
          vcf_reader.header.__str__()
          .rsplit("\n", 2)[0]
          + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        sev_info = pd.DataFrame()
    else:
        sev_header = VcfReader(sev_path).header.__str__()
        sev_info = (pd
          .read_csv(
            sev_path,
            sep="\t",
            comment="#",
            header=None,
            index_col=[0,1,3,4]))
        sev_info.index = sev_info.index.set_levels(
          sev_info.index.levels[0].astype(str),
          level=0)
    sev_lines = list()
    all_anns = list()

    joined_markers = list()
    for marker in markers:
        if not joined_markers:
            joined_markers.append(list(marker) + [marker[1]])
            continue
        if (
          joined_markers[-1][0] == marker[0] and
          marker[1] - joined_markers[-1][2] < 1000):
            joined_markers[-1][2] = marker[1]
        else:
            joined_markers.append(list(marker) + [marker[1]])

    progress(task_id, "Fetching marker data", bold=False)
    for marker in joined_markers:
        marker = tuple(marker)
        vcf_reader.set_region(
          marker[0],
          marker[1] - 1,
          marker[2])
        for rec in vcf_reader:
            if (rec.CHROM, rec.POS, rec.REF, ",".join(rec.ALT)) in sev_info.index:
                continue
            try:
                _ = rec.CSQ
            except:
                continue

            anns = list()
            for ann in rec.CSQ:
                # Skip variants not affecting the protein sequence
                if not ann["HGVS.p"]:
                    continue
                if "synonymous_variant" in ann["Annotation"]:
                    continue

                t_id = ann["Feature_ID"]
                seq_path = os.path.join(
                  prot_dir,
                  "prot",
                f"{t_id}.fasta")
                t_seq = SeqIO.read(
                  seq_path,
                  "fasta")

                # HGVS.p has format p.VAR
                var = ann["HGVS.p"][2:]

                # Convert 3-letter amino acid codes to 1-letter
                var = re.sub(
                  "([A-Z][a-z]{2})",
                  lambda x: convertAA[x.group()],
                  var)

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
                    end_pos = int(end_var[1:]) if end_var[1:] else len(t_seq.seq)

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
                    alt = str(t_seq.seq[start_pos - 1:end_pos]) * 2

                elif re.match("^[A-Z]\d+\_[A-Z]\d+dup$", var):
                    start_var, end_var = var[:-3].split("_")
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref, end_pos = end_var[0], int(end_var[1:])
                    alt = str(t_seq.seq[start_pos - 1:end_pos]) * 2

                # Stop gains (treat as deletion of residues after the new stop)
                elif re.match("^[A-Z]\d+\*$", var):
                    start_var = var[:-1]
                    start_ref, start_pos = start_var[0], int(start_var[1:])
                    end_ref = ""
                    end_pos = len(t_seq.seq) - 1
                    alt = ""

                # Cannot predict effects of start_gain/start_loss/stop_loss/frameshift
                elif any (s in var for s in ("fs", "*", "?")):
                    #print("Skipping", var)
                    continue

                else:
                    progress(task_id, f"Unrecognised variant: {var}")
                    continue

                anns.append(
                  [rec.CHROM,
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
                   alt])

            if not anns:
                continue
            all_anns += anns

    ann_df = pd.DataFrame(
      all_anns,
      columns=[
        "Chrom", "Pos", "Ref", "Alt", "Allele", "Feature_ID",
        "Gene_ID", "Gene_Name", "Annotation", "Var", "Start_ref",
        "Start_pos", "End_ref", "End_pos", "Subs"])
    ann_df["Supporting_set"] = ""

    for t_id, anns in ann_df.groupby("Feature_ID", sort=False):
        progress(task_id, t_id, bold=False)
        anns = anns.reset_index(drop=True)
        full_set = set()

        progress(task_id, "Fetching homologous sequences", newline=False, bold=False)
        # Get transcript sequence
        seq_path = os.path.join(
          prot_dir,
          "prot",
          f"{t_id}.fasta")
        if not os.path.exists(seq_path):
            continue
        t_seq = SeqIO.read(
          seq_path,
          "fasta")

        blastp_out_fn = os.path.join(
          prot_dir,
          "blastp_out",
          f"{t_id}.blastp_out.txt")
        blastp_out_df = (pd
          .read_csv(
            blastp_out_fn,
            header=None,
            sep="\t")
          .set_axis(
            ["sseqid", "sallseqid", "qstart", "qend", "evalue"],
            axis=1))

        # Extract sequences of aligned hits
        hits_out = os.path.join(prot_dir, "aligned_hits.txt")
        hits_fasta = os.path.join(prot_dir, "aligned_hits.fasta")
        (pd
          .Series(
            blastp_out_df["sallseqid"].unique())
          .to_csv(
            hits_out,
            header=False,
            index=False))

        old_wd = os.getcwd()
        os.chdir(os.path.join(data_dir, genome, "prot_db"))
        popen(
          task_id,
          [blastdbcmd_path,
           "-db", config["db_prefix"],
           "-entry_batch", "aligned_hits.txt",
           "-out", "aligned_hits.fasta",
           "-outfmt", "%f"])

        # Cluster by similarity
        popen(
          task_id,
          [cdhit_path,
          "-i", "aligned_hits.fasta",
          "-o", "clust",
          "-c", str(seq_identity),
          "-g", str(cluster_to_most_similar)],
          show_stdout=False)
        with open("clust.clstr") as clust_in:
            clust_raw = clust_in.read().strip().split("\n>")

        hits = {
          rec.name: rec
          for rec in SeqIO.parse("aligned_hits.fasta", "fasta")}

        clust_seqs = list()
        clust_evalues = list()
        for cluster in clust_raw:
            # Format: >Cluster n
            #         <cluster sequences> 
            #         >Cluster n+1
            clust_lines = [
              line for line
              in cluster.strip().split("\n")[1:]]
            # Format: <whitespace> >seq_id...at x% <identity>
            seq_names = [
              line.split(">")[1].split("...")[0]
              for line in clust_lines]
            clust_seqs.append(seq_names)
            mean_evalue = blastp_out_df[
              blastp_out_df["sseqid"].isin(seq_names)]["evalue"].mean()
            clust_evalues.append(
              (seq_names, 
               mean_evalue))
        clusters = [
          cluster[0] for cluster
          in sorted(clust_evalues, key=lambda x: x[1])]

        os.chdir(old_wd)
        progress(task_id, "\r", newline=False, bold=False)

        for idx in range(anns.shape[0]):
            start_pos, end_pos = anns.loc[idx, ["Start_pos", "End_pos"]]
            supp_set = list()

            matching_seqs = set(
              blastp_out_df[
                (blastp_out_df["qstart"] <= start_pos) &
                (blastp_out_df["qend"] >= end_pos)]["sseqid"])
            for cluster in clusters:
                filt_cluster = [
                  seq for seq in cluster
                  if seq in matching_seqs]
                if filt_cluster:
                    supp_set.append(filt_cluster)
                    if len(supp_set) == 45:
                        break
            anns.at[idx, "Supporting_set"] = supp_set
            full_set.update(sum(supp_set, []))

        if len(full_set) == 0:
            continue

        supp_set_fn = os.path.join(
          prot_dir,
          "supporting_set.fasta")
        SeqIO.write(
          [hits[seq_name] for seq_name in full_set],
          supp_set_fn,
          "fasta")

        try:
            p = subprocess.run(
              [needle_path,
               "-asequence", seq_path,
               "-bsequence", supp_set_fn,
               "-gapopen", needle_gapopen,
               "-gapextend", needle_gapextend,
               "-endweight", "false",
               "-aformat", "score",
               "-stdout",
               "-auto"],
              check=True,
              capture_output=True,
              text=True)
        except subprocess.CalledProcessError as e:
            progress(
              task_id,
              f"Alignment for ref-{t_id} failed: {e}")
            continue
        else:
            query_needle_out = p.stdout

            query_needle_df = (pd
              .read_csv(
                io.StringIO(query_needle_out),
                header=None,
                sep=" ",
                comment="#")
              .set_axis(
                ["qseqid", "sseqid", "length", "score"],
                axis=1)
              .drop(
                ["qseqid", "length"],
                axis=1))
            query_needle_df["score"] = (pd
              .to_numeric(query_needle_df["score"].str.strip("()")))

        pos_groups = anns.groupby(["Chrom", "Pos"])
        num_markers = len(pos_groups)
        for idx, (marker, pos) in enumerate(pos_groups, start=1):
            progress(
              task_id,
              f"\rVariant {idx} of {num_markers}",
              newline=False)
            sev_anns = list()
            for var in pos.to_dict(orient="records"):
                if len(var["Supporting_set"]) == 0:
                    continue

                mutated_seq = Seq.MutableSeq(t_seq.seq)
                if (
                  mutated_seq[var["Start_pos"] - 1] != var["Start_ref"] or
                  (var["End_ref"] and mutated_seq[var["End_pos"] - 1] != var["End_ref"])):
                    progress(
                      task_id,
                      "Mismatch between variant and sequence " +\
                      f"in {t_id}, {':'.join(list(map(str, marker)))}, {var['Var']}")
                    continue

                mutated_seq[var["Start_pos"] - 1:var["End_pos"]] = var["Subs"]
                mutated_seq_fn = os.path.join(
                  prot_dir,
                  "mutated.fasta")
                SeqIO.write(
                  [SeqRecord.SeqRecord(mutated_seq)],
                  mutated_seq_fn,
                  "fasta")

                # Get alignment for the mutated sequence
                try:
                    p = subprocess.run(
                      [needle_path,
                       "-asequence", mutated_seq_fn,
                       "-bsequence", supp_set_fn,
                       "-gapopen", needle_gapopen,
                       "-gapextend", needle_gapextend,
                       "-endweight", "false",
                       "-aformat", "score",
                       "-stdout",
                       "-auto"],
                      check=True,
                      capture_output=True,
                      text=True)
                except subprocess.CalledProcessError as e:
                    progress(
                      task_id,
                      f"Alignment for {t_id}-{var['Subs']} failed: {e}")
                    continue
                else:
                    mut_needle_out = p.stdout
                    mut_needle_df = (pd
                      .read_csv(
                        io.StringIO(mut_needle_out),
                        header=None,
                        sep=" ",
                        comment="#")
                      .set_axis(
                        ["qseqid", "sseqid", "length", "score"],
                        axis=1)
                      .drop(
                        ["qseqid", "length"],
                        axis=1))
                    mut_needle_df["score"] = (pd
                      .to_numeric(mut_needle_df["score"].str.strip("()")))

                if not query_needle_df.shape[0] == mut_needle_df.shape[0]:
                    progress(
                      task_id,
                      f"Needle output mismatched in length: {t_id}-{var['Subs']}.")
                    continue

                mut_needle_df["delta"] = (
                  mut_needle_df["score"]
                  - query_needle_df["score"])

                avg_deltas = list()
                for cluster in var["Supporting_set"]:
                    deltas = mut_needle_df[mut_needle_df["sseqid"].isin(cluster)]["delta"]
                    avg_deltas.append(deltas.mean())
                    sev_score = (
                      round(statistics.mean(avg_deltas), 4)
                      if len(avg_deltas) > 0 else "-")
                    sev_stdev = (
                      round(statistics.stdev(avg_deltas), 4)
                      if len(avg_deltas) > 1 else "-")

                sev_anns.append(
                  "|".join(map(str,
                    (var["Allele"],
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
                     ("-" if sev_score == "-" 
                      else "DELETEROUS_ALT" if sev_score < -2.5
                      else "TOLERATED" if sev_score <= 2.5
                      else "DELETERIOUS_REF")))))

            if sev_anns:
                sev_lines.append(
                  (var["Chrom"],
                   var["Pos"],
                   ".",
                   var["Ref"],
                   var["Alt"],
                   ".",
                   ".",
                   "SEV=" + ",".join(sev_anns)))
        progress(task_id, "")

    new_sev = pd.DataFrame(sev_lines)
    if not sev_info.empty:
        sev_info = sev_info.reset_index()[[0,1,2,3,4,5,6,7]]
    if not new_sev.empty:
        sev_info = (pd
          .concat([
            sev_info,
            new_sev]))
        sev_info[0] = pd.Categorical(
          sev_info[0],
          categories=list(pysam.VariantFile(vcf_path).header.contigs),
          ordered=True)
        sev_info = sev_info.sort_values(by=[0,1])
        sev_output = io.StringIO()
        sev_info.to_csv(
          sev_output,
          sep="\t",
          header=False,
          index=False)
        with master_lock:
            lock = get_lock(sev_path)
        with lock:
            with pysam.BGZFile(sev_path, "wb") as sev_out:
                _ = sev_out.write(sev_header.encode("utf-8"))
                _ = sev_out.write(sev_output.getvalue().encode("utf-8"))
            pysam.tabix_index(
              sev_path,
              preset="vcf",
              force=True)


def fetch_sev():
    pass