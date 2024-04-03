import os
import sys
import copy
import gffutils
import subprocess
import pandas as pd
from Bio import SeqIO
#from io import StringIO
from pyfaidx import Fasta
from parseVcf import VcfReader
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation

from pymemesuite.fimo import FIMO
from pymemesuite.common import MotifFile, Sequence


# Step 1: Read cluster information
clusters_file = "../../WGCNA/gene_modules.txt"
clusters = pd.read_csv(clusters_file, sep="\t", header=None, names=["gene_id", "cluster_id"])
clusters["gene_id"] = clusters["gene_id"].str[:-2]

# Step 2: Access gene information from the GFF database
# Assuming you have a GFF database named 'gff_db'
# and a genomic fasta file named 'genome.fa'
gff_db_path = "../ITAG2.4_gene_models.db"
gff_db = gffutils.FeatureDB(gff_db_path)

seq_path = "../../../Genome_data/Genomes/Solyc2.4/S_lycopersicum_chromosomes.2.40.fa"
seq = Fasta(seq_path)

#vcf_path = "../../TFBS/var_demo.vcf.gz"
vcf_path = "../../TFBS/snpEff.filtered.SL2.40.info_only.vcf.gz"
vcf_reader = VcfReader(vcf_path)

# No. of variants: 137256

motifs_fn = "../../TFBS/Sly_TF_binding_motifs.meme"
ann_fn = "../../TFBS/promoter_annotations.ann"

promoter_size = (1000, -500)

# Load motifs
motifs = dict()
with MotifFile(motifs_fn) as motif_file:
    while True:
        motif = motif_file.read()
        if motif:
            motifs[motif.accession.decode()] = motif
        else:
            break
    

with open(ann_fn, "w") as ann_file:
    for cluster_id, cluster_data in clusters.groupby("cluster_id"):
        gene_ids = cluster_data["gene_id"]
        promoter_sequences = dict()

        for gene_id in gene_ids:
            try:
                record = gff_db[gene_id]
            except:
                continue
            if record.strand == "-":
                continue

            # Get promoter upstream of transcription start
            # Get promoter upstream of transcription start
            # (for genes with multiple transcripts, get each promoter)
            transcripts = gff_db.children(record["ID"][0],
                                          featuretype="mRNA")
            for t in transcripts:
                transcript_id = t["ID"][0]
                # Extract promoter sequence
                promoter_start, promoter_end = (
                                                (max(1, t.start - promoter_size[0]),
                                                 t.start - promoter_size[1] - 1
                                                 )
                                                if t.strand == "+" else
                                                (max(1, t.end + promoter_size[1]),
                                                 t.end + promoter_size[0] - 1
                                                )
                                               )
                promoter_seq = str(seq[record.chrom][promoter_start - 1:promoter_end - 1])

                # Create a SeqRecord for the promoter sequence
                seq_record = SeqRecord(MutableSeq(promoter_seq), id=transcript_id, description=f"{record.chrom}:{promoter_start}..{promoter_end};parent={gene_id}")
                promoter_sequences[transcript_id] = seq_record

        # Save promoter sequences to a FASTA file
        fasta_fn = f"../../TFBS/promoters_{cluster_id}.fasta"
        SeqIO.write(promoter_sequences.values(), fasta_fn, "fasta")

        # Step 4: Run MEME-SEA
        out_dir = f"../../TFBS/sea_out_{cluster_id}"
        thresh = 0.001/len(promoter_sequences)
        if 0:
            command = ["sea", "--p", fasta_fn, "-m", motifs_fn, "--oc", out_dir, "--qvalue", "--thresh", str(thresh)]
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Wait for the process to complete and get the output
            stdout, stderr = process.communicate()

            # Print the output and return code
            print("Output:", stdout.decode())
            print("Error:", stderr.decode())
            print("Return Code:", process.returncode)

        sea_result = pd.read_csv(f"{out_dir}/sequences.tsv", sep="\t")
        # Get only true positives
        sea_result = sea_result[sea_result["seq_Class"] == "tp"]

    #    out_dir = f"../../TFBS/fimo_out_{cluster_id}"
    #    if 0:
    #        command = ["fimo", "--oc", out_dir]
    #        for motif_id in set(sea_result["motif_ID"]):
    #            command += ["--motif", motif_id]
    #        command += [motifs_fn, fasta_fn]
    #        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Wait for the process to complete and get the output
    #        stdout, stderr = process.communicate()

            # Print the output and return code
    #        print("Output:", stdout.decode())
    #        print("Error:", stderr.decode())
    #        print("Return Code:", process.returncode)

    #    fimo_result = pd.read_csv(f"{out_dir}/sequences.tsv", sep="\t")

        fimo = FIMO(both_strands=True)

        # Step 5: Run FIMO for each motif
        for motif_id, t_id in zip(sea_result["motif_ID"],
                                  sea_result["seq_ID"]):

            motif = motifs[motif_id]
            promoter = promoter_sequences[t_id]
            meme_sequence = Sequence(str(promoter.seq), name=promoter.id.encode(), description=promoter.description.encode())

            region, parent = promoter.description.split(";")
            parent_gene = parent.split("parent=")[1]
            chrom, span = region.split(":")
            start, end = map(int, span.split(".."))
    #        print(start, end)

            pattern = fimo.score_motif(motif, [meme_sequence], motif_file.background)
    #        print(len([m for m in pattern.matched_elements if m.qvalue < 0.001 and m.strand == "+"]),
    #              len([m for m in pattern.matched_elements if m.qvalue < 0.001 and m.strand == "-"]),
    #              len([m for m in pattern.matched_elements if m.qvalue >= 0.001]))
            for m in pattern.matched_elements:
                if m.qvalue >= thresh:
                    continue
                if m.strand == "-":
                    continue

                # matched_element numbering starts from 1,
                # python numbering (used in Bio.Seq indexing) starts from 0
                if m.strand == "+":
                    matched_start, matched_stop = m.start, m.stop
                else:
                    matched_start, matched_stop = m.stop, m.start
                matched_sequence = promoter[matched_start-1:matched_stop]

    #            print()
    #            for i, pos in enumerate(range(start + m.start - 1, start + m.stop)):
    #                print(pos, matched_sequence[i])

                # VcfReader numbering starts from 1 and doesn't consider first position
                vcf_reader.set_region(chrom, start + matched_start - 1, start + matched_stop - 1)
                for var_record in vcf_reader:
                    annotations = list()
                    var_sequence = copy.deepcopy(matched_sequence).seq

                    offset_start = abs(min(0, var_record.POS - (start + matched_start - 1)))
                    offset_end = - max(0, var_record.POS + len(var_record.REF) - 1 - (start + matched_start - 1))
    #                print(offset_start, offset_end)

    #                print(str(var_sequence[var_record.POS - (start + m.start - 1) + offset_start : var_record.POS - (start + m.start - 1) + len(var_record.REF) + offset_end]))
    #                print(var_record.REF[offset_start: len(var_record.REF) + offset_end])
                    assert(str(var_sequence[var_record.POS - (start + matched_start - 1) + offset_start : var_record.POS - (start + matched_stop - 1) + len(var_record.REF) + offset_end]).upper() == var_record.REF[offset_start: len(var_record.REF) + offset_end])
    #                continue
                    for alt_allele in var_record.ALT:
                        var_sequence[var_record.POS - (start + matched_start - 1) + offset_start : var_record.POS - (start + matched_start - 1) + len(var_record.REF) + offset_end] = alt_allele

                        meme_sequence = Sequence(str(var_sequence))
                        var_pattern = fimo.score_motif(motif, [meme_sequence], motif_file.background)

                        scores = [vm.score for vm in var_pattern.matched_elements if vm.qvalue < thresh]
    #                continue

                        score = max(scores) - m.score if scores else "Activity_lost"

                        annotations.append("|".join((alt_allele,
                                                     "promoter_region",
                                                     "MODIFIER",
                                                     parent_gene,
                                                     t_id,
                                                     motif_id,
                                                     motif.consensus,
                                                     m.strand,
                                                     str(score),
                                                     str(matched_sequence.seq),
                                                     str(var_sequence),
                                                     "",
                                                     "",
                                                     "",
                                                     "",
                                                     ""
                                                    )))
                    ann_file.write("\t".join((chrom,
                                              str(var_record.POS),
                                              ".",
                                              var_record.REF,
                                              ",".join(var_record.ALT),
                                              ".",
                                              ".",
                                              "ANN=")))
                    ann_file.write(",".join(annotations))
                    ann_file.write("\n")
                                                 
        break
    #            print(
    #                m.source.accession.decode(),
    #                m.start,
    #                m.stop,
    #                m.strand,
    #                m.score,
    #                m.pvalue,
    #                m.qvalue
    #            )
                #print(annotations)
                #sys.exit()

    #        for record in vcf:
    #            chrom = record.CHROM
    #            pos = record.POS

    #            # Check if the variant is within the motif
    #            motif_location = fimo_result.get_location(chrom, pos)
    #            if motif_location:
    #                motif_seq = str(Fasta("promoters.fasta")[motif_location.sequence_name][motif_location.start:motif_location.end])

    #                # Create mutated sequence
    #                mutated_seq = promoter_seq[:motif_location.start] + variant["ALT"] + promoter_seq[motif_location.end:]

    #                # Run FIMO for the mutated sequence
    #                mutated_fimo_result = fimo(motif_file, {"query": mutated_seq, "source": "text"}, output_dir="mutated_fimo_output")

    #                # Compare scores
    #                wild_type_score = fimo_result.get_scores()[0]
    #                mutated_score = mutated_fimo_result.get_scores()[0]

    #                # Report the difference in scores
    #                score_difference = mutated_score - wild_type_score
    #                print(f"Cluster {cluster_id}, Gene {gene_id}, Motif {motif_file}: Score Difference = {score_difference}")
