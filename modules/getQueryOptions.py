#!/usr/bin/env python
import os
import re
import json
import gffutils
import numpy as np
import pandas as pd
from .parseVcf import VcfReader
#from gffutils.inspect import inspect

"""
Haplosmasher v2.0 input processing functions
Author: Peter Poláček

Usage: To be called from within HaploSmasherApp.py or other scripts
"""

# Define working directories
root_dir = os.getcwd()
data_dir = os.path.join(root_dir, "data")


def getDbDir():
    """Returns path to gffutils DB files"""
    with open("res/db_dir.txt", "r") as db_file:
        lines = db_file.read().split("\n")

    lines = [line for line in lines
             if line and
             not line.startswith("#")]

    # Return error unless exactly one line found and path is valid
    if len(lines) != 1 or not os.path.isdir(lines[0]):
        raise OSError("Invalid path to GFFutils databases.")

    return lines[0]


def getFile(genome, extension):
    """Parses directory of a given genome for a query file

    genome: name of folder in the data/ directory
    extension: str, one of "vcf", "gff", "pep", "seq", "phe"

    returns: path of file
    """
    genome_dir = f"data/{genome}/"

    if extension == "gff":
        file = [f for f in os.listdir(genome_dir)
                if f.casefold().endswith((".gff",
                                          ".gff3",
                                          ".gtf",
                                          ".gff.gz",
                                          ".gff3.gz",
                                          ".gtf.gz"))]

    elif extension == "vcf":
        file = [f for f in os.listdir(genome_dir)
                if f.casefold().endswith((".vcf.gz",
                                          ".bcf.gz"))]

    elif extension == "pep":
        file = [f for f in os.listdir(genome_dir)
                if f.casefold().endswith(("pep.fasta",
                                          "pep.fa",
                                          "pep.fna"))]

    elif extension == "seq":
        file = [f for f in os.listdir(genome_dir)
                if f.casefold().endswith(("seq.fasta",
                                          "seq.fa",
                                          "seq.fna"))]

    elif extension == "phe":
        file = [f for f in os.listdir(genome_dir)
                if f.casefold().endswith((".phe"))]
        

    if len(file) != 1:
        raise OSError(f"{extension.upper()} file for {genome} " +\
                       "missing or multiple files found.")

    return genome_dir + file[0]


def getVariantConsequences(impact_tag=True):
    """Creates a list of tuples of variant consequences

    impact_tag: bool (default=True), add impact to display name

    returns: list, [(id_name, display_name)]
    """
    csq_list = list()

    with open("static/variant_consequences.tsv", "r") as csq_in:
        for line in csq_in.readlines():
            id_name, display_name, impact = line.split("\t")[:3]

            if impact_tag == True:
                # Add a tag e.g. "[HIGH] start loss"
                csq_dict = {"HIGH": "H", "MODERATE": "M",
                            "LOW": "L", "MODIFIER": "-"}
                display_name = f"{csq_dict[impact]} {display_name}"

            csq_list.append((id_name, display_name))

    return csq_list


def convertSampleNames(genome):
    """Creates a dictionary of sample names and IDs for query genome

    genome: name of folder in the data/ directory

    returns: dict(), {sample_ID: sample_name}
    """
    sample_dict = dict()
    samplenames_file = f"data/{genome}/sample_names.tsv"
    if os.path.exists(samplenames_file):
        with open(samplenames_file, "r") as samples_in:
            for line in samples_in.readlines():
                (sample_id,
                 sample_name) = line.split("\t")[:2]
                sample_dict[sample_id] = sample_name

    return sample_dict


def getGenomeData(genome):
    """Returns information a given genome

    returns: dict, {"Chromosomes": list(),
                    "Samples": list(),
                    "Groups": [([samples], group_name)],
                    "Sample_ID": str}

    Keeps data in genome_data.json file in genome source directory in data/
    """
    if not os.path.exists("res/gffutils_db/"):
        os.mkdir("res/gffutils_db/")

    try:
        db_dir = getDbDir()
    except OSError as error:
        raise OSError(error.args[0])

    # Create a gffutils database if none exists
    if not os.path.exists(f"{db_dir}/{genome}.db"):
        try:
            gff_file = getFile(genome, "gff")
        except OSError as error:
            raise OSError(error)
        gffutils.create_db(gff_file,
                           dbfn=f"{db_dir}/{genome}.db",
                           merge_strategy="merge")

    # Create a genome data file if none exists
    if not os.path.exists(f"data/{genome}/genome_data.json"):
        genome_data = dict()

        db = gffutils.FeatureDB(f"{db_dir}/{genome}.db")
        # Create list of chromosomes, sorted by trailing digits
#        chromosomes = list(inspect(db, look_for=["chrom"])["chrom"].keys())
#        chromosomes.sort(key=lambda x: int(re.search('\d+$', x)[0]))
        seq_file = getFile(genome, "seq")
        seq_df = pd.read_csv(seq_file + ".fai", sep="\t", header=None)
        chromosomes = seq_df.iloc[:, 0].tolist()
        sizes = seq_df.iloc[:, 1].tolist()
        #chromosomes.sort(key=lambda x: int(re.search('\d+$', x)[0]))

        genome_data["Chromosomes"] = chromosomes

        # Get information about chromosome sizes,
        # offsets for calculating plot marker positions
        # and ticks for calculating plot label positions
        offset = 0
        genome_data["Chromosome_data"] = {}
        for chrom, size in zip(chromosomes, sizes):
            tick = offset + size // 2
            genome_data["Chromosome_data"][chrom] = {"size": size,
                                                     "offset": offset,
                                                     "tick": tick}
            offset += size

        # Get first gene ID in the file
        for record in db.region(region=chromosomes[0],
                                featuretype="gene"):
            genome_data["Sample_ID"] = record.id
            break

        # Get a list of samples in the VCF file
        try:
            vcf_file = getFile(genome, "vcf")
        except OSError as error:
            raise OSError(error)
        vcf_reader = VcfReader(vcf_file)
        samples = vcf_reader.header.samples

        sample_dict = convertSampleNames(genome)
        genome_data["Samples"] = [(sample, sample_dict.get(sample, sample))
                                  for sample in samples]

        # Initialise empty list for groups
        genome_data["Groups"] = dict()

        # Save genome data to JSON file
        with open(f"data/{genome}/genome_data.json", "w") as data_out:
            json.dump(genome_data, data_out)#, indent=4)

    else:
        # Load genome data from file
        with open(f"data/{genome}/genome_data.json", "r") as data_in:
            genome_data = json.load(data_in)

    return genome_data


if __name__ == "__main__":
    pass
