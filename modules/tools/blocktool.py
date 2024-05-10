#!/usr/bin/env python
import os
import sys
from ast import literal_eval
import pathlib

from argparse import (ArgumentParser,
                      RawTextHelpFormatter,
                      ArgumentTypeError,
                      SUPPRESS)

from subprocess import run

from collections import Counter
import pandas as pd

#sys.tracebacklimit = 0
#sys.stdout = open(os.devnull, "w")

"""
Unified tool for haplotype construction and genotype calling

Usage: python3 blocktool.py [-h]
"""

plink_path = os.path.join("modules", "tools", "plink1.9")

##############################################################################
### Functions for processing input arguments
##############################################################################

def is_dir(path):
    """Checks if input is a valid path, creates if doesn't exist"""
    if os.path.isdir(path):
        return path
    else:
        try:
            pathlib.Path(path).mkdir(parents=True)
            print(f"Created directory {path} for output.")
            return path
        except OSError:
            raise ArgumentTypeError(f"{path} is not a valid directory.")

def above_zero(value):
    """Checks if value is a valid integer higher than 0"""
    if type(literal_eval(value)) == int and int(value) > 0:
        return int(value)
    else:
        raise ArgumentTypeError("Value must be an integer higher than 0.")

def above_one(value):
    """Checks if value is a valid integer higher than 1"""
    if type(literal_eval(value)) == int and int(value) > 1:
        return int(value)
    else:
        raise ArgumentTypeError("Value must be an integer higher than 1.")

def zero_one(value):
    """Checks if value is a valid float >= 0 and < 1"""
    if (type(literal_eval(value)) in (int, float) and
        (0 <= float(value) < 1)):
        return float(value)
    else:
        raise ArgumentTypeError("Value must be <= 0 < 1.")

##############################################################################
### Input argument parser
##############################################################################

class HelpfulParser(ArgumentParser):
   def error(self, message):
      self.print_help()
      sys.stderr.write('error: %s\n' % message)
      sys.exit(2)

parser = HelpfulParser(description="Unified utility for haplotype block "
                                    "construction and genotype calling",
                        epilog="Author: Peter Poláček",
                        formatter_class=RawTextHelpFormatter,
                        usage=SUPPRESS)

subprs = parser.add_subparsers(title="Commands",
                               dest="command",
                               required=True)

###
### Temporary parser, checks if command is 'block'
### Checks the chosen type to make sub-options required
###

block_decider = ArgumentParser("block decider",
                               add_help=False,
                               usage=SUPPRESS)

block_arg = block_decider.add_argument("--type",
                           choices=["gab", "4gam", "bpwin", "snpwin"],
                           help="Type of block construction method to use:\n"
                                "gab     groups markers "
                                "with LD above a threshold\n"
                                "4gam    groups markers with with "
                                "no evidence of recombination\n"
                                "bpwin   groups markers based on "
                                "a sliding window of a given bp size\n"
                                "snpwin  groups markers based on a sliding "
                                "window with a given number of markers")

deciding_args = block_decider.parse_known_args()[0]

block_pr = subprs.add_parser("block",
                             help="Consctruct haplotype blocks",
                             parents=[block_decider],
                             formatter_class=RawTextHelpFormatter,
                             usage=SUPPRESS)

call_pr = subprs.add_parser("call",
                            help="Genotype based on blocks",
                            formatter_class=RawTextHelpFormatter,
                            epilog="Note: root filename of each genotype file "
                                   "and the corresponding block file must "
                                   "be identical.",
                            usage=SUPPRESS)

block_arg.required = True

###
### 'block' command arguments
###

block_pr.add_argument("--gab-max-kb",
                     type=above_zero,
                     required=(deciding_args.type == "gab"),
                     help="Maximum span of gabriel blocks, in kb")

block_pr.add_argument("--bp-span",
                      type=above_one,
                      required=(deciding_args.type == "bpwin"),
                      help="Sliding window span for bp blocking, in bases")
block_pr.add_argument("--snp-span",
                      type=above_one,
                      required=(deciding_args.type == "snpwin"),
                      help="Sliding window span for marker-based blocking")

block_pr.add_argument("--bp-step",
                      type=above_zero,
                      required=(deciding_args.type == "bpwin"),
                      help="Sliding window step for bp blocking, in bases")
block_pr.add_argument("--snp-step",
                      type=above_zero,
                      required=(deciding_args.type == "snpwin"),
                      help="Sliding window step for marker-based blocking")

block_pr.add_argument("--gt-dir",
                     type=is_dir,
                     required=True,
                     help="Directory containing input genotype files\n"
                          "gab    .bed, .bim, .fam files (plink --make-bed)\n"
                          "4gam   .tped file (plink --recode transpose)\n"
                          "bpwin  .map file (plink --recode)\n"
                          "snpwin .map file (plink --recode)")

block_pr.add_argument("--out",
                      type=is_dir,
                      help="Directory for block file output\n"
                           "Leave empty to save to --gt-dir")

###
### 'call' command arguments
###

call_pr.add_argument("--format",
                     choices=["biallelic",
                              "per-hap",
                              "per-loc"],
                     default="biallelic",
                     nargs="+",
                     help="GWAS model for genotype output, space-separated\n"
                          "biallelic  One call per locus, allele_1 = major "
                          "haplotype, allele_2 = all others,\n           "
                          "call represents zygosity at locus for alelles\n"
                          "per-hap    One call per haplotype, call represents "
                          "allele dosage for given haplotype\n"
                          "per-loc    Two calls per haplotype block, one per "
                          "chromatid; call represents haplotype number\n"
                          "Default = biallelic")

call_pr.add_argument("--block-dir",
                     type=is_dir,
                     required=True,
                     help="Directory containing input block files")

call_pr.add_argument("--gt-dir",
                     type=is_dir,
                     required=True,
                     help="Directory containing input .tped and .tfam files")

call_pr.add_argument("--out",
                     type=is_dir,
                     help="Directory for genotype and info file output\n"
                          "Leave empty to save to --gt-dir")

call_pr.add_argument("--transpose",
                     choices=["biallelic",
                              "per-hap",
                              "per-loc"],
                     default="",
                     nargs="+",
                     help="Transpose output (accessions in rows)"
                          " for chosen formats, space-separated")

call_pr.add_argument("--maf",
                     type=zero_one,
                     default=0.05,
                     help="Filter haplotypes in output tables "
                          "on minor allele frequency\n"
                          "0 >= maf < 1, default = 0.05"
                          )

call_pr.add_argument("--min-hap",
                     type=int,
                     default=1,
                     help="Minimum number of haplotypes above the "
                          "MAF threshold for a block to be retained\n"
                          "Default = 1"
                          )

call_pr.add_argument("--group-minor",
                     action="store_true",
                     help="Group haplotypes below MAF threshold into one\n"
                          "Group still checked against the MAF threshold")

call_pr.add_argument("--include-singles",
                     action="store_true",
                     help="Include markers not placed in a haplotype")

call_pr.add_argument("--tab",
                     action="store_true",
                     help="Output tab-separated tables, "
                          "instead of comma-separated")

##############################################################################
### Process input arguments
##############################################################################

# Print help if no arguments provided
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()
elif len(sys.argv) == 2 and sys.argv[1] == "block":
    block_pr.print_help()
    parser.exit()
elif len(sys.argv) == 2 and sys.argv[1] == "call":
    call_pr.print_help()
    parser.exit()

args = parser.parse_args()

try:
    if args.block_dir:
        args.block_dir = args.block_dir.strip("/")
# When using 'block' option, this argument will be missing
except AttributeError:
    pass

if not args.out:
    args.out = args.gt_dir

args.gt_dir = args.gt_dir.strip("/")
args.out = args.out.strip("/")

# Get a list of input genotype (and block) files
if args.command == "block":
    ext_dict = {"gab":    ".bed",
                "4gam":   ".tped",
                "bpwin":  ".map",
                "snpwin": ".map"}

    gt_files = [file.rsplit(".", 1)[0] for file in os.listdir(args.gt_dir)
                if file.endswith(ext_dict[args.type])]
    if not gt_files:
        raise ArgumentTypeError(f"No {ext_dict[args.type]} files " +\
                                f"found in '{args.gt_dir}'.")

    print(f"Blocking these files in '{args.gt_dir}':")
    print("\n".join(gt_files))

elif args.command == "call":
    block_files = [file for file in os.listdir(args.block_dir)
                   if file.endswith((".blocks",
                                     ".4gamblocks",
                                     ".bpwinblocks",
                                     ".snpwinblocks"))]
    if not block_files:
        raise ArgumentTypeError(f"No block files found " +\
                                f"in '{args.block_dir}'.")

    gt_files = [file.rsplit(".", 1)[0] for file in os.listdir(args.gt_dir)
                 if file.endswith(".tped") and
                 os.path.isfile(args.gt_dir + "/" +\
                                file.rsplit(".", 1)[0] + ".tfam")]

    file_pairs = list()
    for gt_file in gt_files:
        for block_file in block_files:
            if block_file.startswith(gt_file):
                file_pairs.append((gt_file, block_file))

    if not file_pairs:
        raise ArgumentTypeError(f"No .tped, .tfam and block files with " +\
                                f"matching names found in '{args.gt_dir}'.")

    print(f"Genotyping these files:")
    for (gt_file, block_file) in file_pairs:
        print(f"{args.gt_dir}/{gt_file} using " +\
              f"{args.block_dir}/{block_file}")

##############################################################################
### Check for unused arguments
##############################################################################

unused_args = list()
if args.command == "block":
    if args.type != "gab" and args.gab_max_kb:
        unused_args.append("--gab-max-kb")

    if args.type != "snpwin" and args.snp_span:
        unused_args.append("--snp-span")
    if args.type != "snpwin" and args.snp_step:
        unused_args.append("--snp-step")

    if args.type != "bpwin" and args.bp_span:
        unused_args.append("--bp-span")
    if args.type != "bpwin" and args.bp_step:
        unused_args.append("--bp-step")

if unused_args:
    print("Unused arguments found:")
    print("\n".join(unused_args))
    print("Check if your input was correct.")

##############################################################################
### Four-gamete block construction method
##############################################################################

def fourgam_blocking(gt_fn,
                     out_fn):
    """Hapotypes acc. to Hudson and Kaplan (Genetics, 111(1), 147-164)"""
    gt_fn += ".tped"
    out_fn += ".4gamblocks"

    gt_file = open(gt_fn, "r")
    out_file = open(out_fn, "w")
    gt_size = os.path.getsize(gt_fn)

    # Read first marker
    line_1 = next(gt_file)
    line_size = sys.getsizeof(line_1)

    line_1 = line_1.rstrip().split()
    marker_1 = line_1[1]
    alleles_1 = line_1[4:]

    block = [marker_1]

    print("Blocking progress:")

    for lc, line in enumerate(gt_file):

        print(f"\r{(line_size * lc) * 100 // gt_size}%", end="")

        # Read second marker to compare with
        line_2 = line.rstrip().split()
        marker_2 = line_2[1]
        alleles_2 = line_2[4:]

        # Check the four-gamete rule
        if len(Counter(zip(alleles_1, alleles_2))) < 4:
            block.append(marker_2)

        else:
            if len(block) > 1:
                out_file.write(" ".join(["*"] + block))
                out_file.write("\n")
            block = [marker_2]

        # Shift markers to compare by 1
        marker_1 = marker_2
        alleles_1 = alleles_2

    print("\r100%")

    gt_file.close()
    out_file.close()

    return out_fn

##############################################################################
### Sliding window block construction method
##############################################################################

def win_blocking(haptype,
                 marker_fn,
                 out_fn,
                 span,
                 step):
    """Haplotypes acc. to the sliding window approach"""
    marker_fn = marker_fn + ".map"
    with open(marker_fn, "r") as markers_in:

        markers = [line.split() for line in
                   markers_in.read().strip().split("\n")]
        markers = [(marker[1], int(marker[3]))
                   for marker in markers]

    out_fn = out_fn + (".bpwinblocks" if haptype == "bpwin"
                       else ".snpwinblocks")
    block_file = open(out_fn, "w")
    block = list()
    # Last block kept to make sure duplicate blocks are not saved
    last_block = ""

    print("Blocking progress:")
    marker_count = len(markers)
    marker_pos = 0

    if haptype == "bpwin":

        chr_span = markers[-1][1] + step
        for block_start in range(0,
                                 chr_span,
                                 step):

            print(f"\r{marker_pos * 100 // marker_count}%", end="")

            block_end = block_start + span
            # Remove markers outside current block from list
            block = [marker for marker in block
                     if marker[1] >= block_start]

            for marker in markers:
                # Keep marker if within current block, otherwise save block
                if marker[1] <= block_end:
                    block.append(marker)

                else:
                    # Make sure to not include duplicate blocks
                    if len(block) > 1 and block != last_block:

                        block_file.write(" ".join(["*"] + [marker[0] for marker
                                                           in block]))
                        block_file.write("\n")

                        last_block = block

                    break

            while markers and markers[0][1] <= block_end:
                markers.pop(0)
                marker_pos += 1

    elif haptype == "snpwin":
        while markers:

            print(f"\r{marker_pos * 100 // marker_count}%", end="")

            block_file.write(" ".join(["*"] + [marker[0] for marker in
                                               markers[:span]]))
            block_file.write("\n")

            for i in range(step):
                if markers:
                    markers.pop(0)
                    marker_pos += 1

    print("\r100%")

    block_file.close()

    return out_fn

##############################################################################
### Block genotyping pipeline
##############################################################################

def block_genotyping(gt_file,
                     block_file,
                     out_dir,
                     out_format,
                     transpose,
                     maf,
                     min_haps,
                     group_minor,
                     include_singles,
                     tab):
    """Genotypes given accessions based on input haplotype blocks"""
    sep = "\t" if tab else ","
    # Get list of blocks and accessions
    with open(block_file, "r") as in_b:
        # Skip leading "*", split on space
        blocks = [l.split()[1:] for l in
                  in_b.read().strip().split("\n")]
        block_count = len(blocks)

    with open(gt_file + ".tfam", "r") as in_acc:
        accessions = [l.split()[1] for l in
                      in_acc.read().strip().split("\n")]
        acc_chromatids = []
        for acc in accessions:
            acc_chromatids.append(acc + "_1")
            acc_chromatids.append(acc + "_2")

    # Create output files
    biallelic_out, perhap_out, perloc_out = list(), list(), list()
    biallelic, per_hap, per_loc = ("biallelic" in out_format,
                                   "per-hap" in out_format,
                                   "per-loc" in out_format)

    block_fn = block_file.rsplit("/", 1)[1].rsplit(".", 1)[0]
    block_info = open(f"{out_dir}/{block_fn}.blockinfo", "w")
    hap_info = open(f"{out_dir}/{block_fn}.hapinfo", "w")

    #hapmap = open(f"{out_dir}/{block_fn}.hapmap", "w")
    hapmap = open(f"{out_dir}/{block_fn}.map", "w")
    locmap = open(f"{out_dir}/{block_fn}.locmap", "w")

    # Write header line of .blockinfo file
    block_info.write(sep.join(["Block_ID",
                                "Span",
                                "Marker_count",
                                "Hap_count"]))
    if maf > 0:
        block_info.write(f"\tHap_count(MAF>{maf})")

    block_info.write("\n")

    # Write header line of .hapinfo file
    hap_info.write(sep.join(["Haplotype_ID",
                              "Alleles",
                              "Frequency"]))
    hap_info.write("\n")

    # Write header line of MAP files
#    hapmap.write(sep.join(["Marker",
#                           "Chromosome",
#                           "Start",
#                           "End"]))
#    hapmap.write("\n")

    locmap.write(sep.join(["Marker",
                            "Chromosome",
                            "Start",
                            "End"]))
    locmap.write("\n")

    df_gt = pd.DataFrame()
    gt_list = list()
    positions = list()
    in_gt = open(gt_file + ".tped", "r")

    print("Genotyping progress:")
    for block_index, block in enumerate(blocks, 1):

        # Print progess status
        print(f"\r{block_index * 100 // block_count}%", end="")

        # Filter data left from previous block
        for series_index in range(len(gt_list)):
            if gt_list[0].name not in block:
                gt_list.pop(0)
                positions.pop(0)
            else:
                break

        # Assemble a list of markers corresponding to the block
        while len(gt_list) < len(block):

            line = next(in_gt).split()
            marker = line[1]
            alleles = line[4:]

            # If marker is not part of any block
            if marker not in block:
                if include_singles:
                    if biallelic:
                        biallelic_out.append(pd.Series(data=alelles,
                                                       name=marker))
                    if per_loc:
                        perloc_out.append(pd.Series(data=alleles,
                                                    name=marker))
                    if per_hap:
                        perhap_out.append(pd.Series(data=alleles,
                                                    name=marker))
                continue

            chrom = line[0]
            positions.append(line[3])
            gt_list.append(pd.Series(data=alleles, name=marker))

        # Create a block ID in the format chrom:start-end
        first_pos, last_pos = positions[0], positions[-1]
        block_id = chrom + ":" + first_pos + "-" + last_pos

        df_gt = pd.concat(gt_list, axis=1)
        # Count occurrences of haplotypes
        haps = Counter(df_gt.to_records(index=False).tolist()).most_common()
        filtered_haps = dict()

        for alleles, count in haps:

            # Filter on MAF
            freq = count / len(acc_chromatids)
            if freq >= maf:
                filtered_haps[alleles] = (len(filtered_haps) + 1,
                                          freq)
            else:
                if group_minor:
                    # Keep all haplotypes below MAF threshold
                    # Add frequencies together
                    alleles = (f"all_at_MAF<{maf}",)
                    filtered_haps.setdefault(alleles,
                                             (len(filtered_haps) + 1,
                                              0)
                                            )
                    filtered_haps[alleles] = (filtered_haps[alleles][0],
                                              filtered_haps[alleles][1] + freq)

        # If some haplotypes were grouped, check if group passes MAF threshold
        if (f"all_at_MAF<{maf}",) in filtered_haps.keys():
            if filtered_haps[(f"all_at_MAF<{maf}",)][1] < maf:
                filtered_haps.pop((f"all_at_MAF<{maf}",))

        # Skip block if not enough haplotypes pass MAF threshold
        if len(filtered_haps) < min_haps:
            continue


        # Make a new block entry to .blockinfo file
        span = int(last_pos) - (int(first_pos) - 1)

        block_entry = [block_id,
                       str(span),
                       str(len(block)),
                       str(len(haps))]
        if maf > 0:
            block_entry.append(str(len(filtered_haps)))

        block_info.write(sep.join(block_entry))
        block_info.write("\n")

        # Make a new entry to .locmap file
        locmap_entry = [block_id,
                        chrom,
                        first_pos,
                        last_pos]
        locmap.write(sep.join(locmap_entry))
        locmap.write("\n")

        # Make a new haplotype entry to .hapinfo file
        # Sort haplotype entries by index (first element in value of each key)
        # Make a new entry to .hapmap file
        for alleles, (hap_index, freq) in sorted(filtered_haps.items(),
                                                 key=lambda i: i[1][0]):

            hap_id = block_id + "_" + str(hap_index)
            allele_string = "".join(str(a) for a in alleles)

            hap_info.write(sep.join([hap_id,
                                     allele_string,
                                     str(freq)]))
            hap_info.write("\n")

#            hapmap_entry = [hap_id,
#                            chrom,
#                            first_pos,
#                            last_pos]
            hapmap_entry = [chrom,
                            hap_id,
                            "0",
                            first_pos]

            hapmap.write(sep.join(hapmap_entry))
            hapmap.write("\n")


        # Get haplotype number for each haplotype
        minor_val = len(filtered_haps) if group_minor else 0
        alleles = [filtered_haps.get(hap, (minor_val,))[0] for hap in
                   [tuple(hap) for hap in df_gt.to_numpy()]]
        alleles = pd.Series(data=alleles, name=block_id)

        if per_loc:
            perloc_out.append(alleles)

        if per_hap:
            # Count number of occurrences per haplotype in each accession
            for index in range(1, len(filtered_haps) + 1):
                hap_id = block_id + "_" + str(index)

                # Get dosage of each haplotype per accession
                alleles_1 = (alleles[0::2] == index).astype(int)
                alleles_1 = alleles_1.reset_index(drop=True)
                alleles_2 = (alleles[1::2] == index).astype(int)
                alleles_2 = alleles_2.reset_index(drop=True)

                perhap_out.append((alleles_1 + alleles_2).rename(hap_id))

        if biallelic:
            # Get zygosity at a locus per accession
            alleles_1 = (alleles[0::2] != 1).astype(int)
            alleles_1 = alleles_1.reset_index(drop=True)
            alleles_2 = (alleles[1::2] != 1).astype(int)
            alleles_2 = alleles_2.reset_index(drop=True)

            biallelic_out.append((alleles_1 + alleles_2).rename(block_id))

    if include_singles:

        while in_gt:
            line = next(in_gt).split()
            marker = line[1]
            alleles = line[4:]

            if biallelic:
                biallelic_out.append(pd.Series(data=alelles,
                                               name=marker))
            if per_loc:
                perloc_out.append(pd.Series(data=alleles,
                                            name=marker))
            if per_hap:
                perhap_out.append(pd.Series(data=alleles,
                                            name=marker))

    print("\r100%")

    if biallelic:
        if biallelic_out:
            biallelic_out = pd.concat(biallelic_out,
                                      axis=1).set_axis(accessions)
            if "biallelic" not in transpose:
                biallelic_out = biallelic_out.T
            biallelic_out.to_csv(f"{out_dir}/{block_fn}.biallelic",
                                 sep=sep)
        else:
            print("No haplotypes to include in .biallelic file.")

    if per_hap:
        if perhap_out:
            perhap_out = pd.concat(perhap_out,
                                   axis=1)
            perhap_markers = perhap_out.columns.tolist()
            perhap_out["IID"] = accessions
            perhap_out[["FID", "PAT", "MAT", "SEX", "PHENOTYPE"]] = ""
            perhap_out = perhap_out[["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] + perhap_markers]
            if "per-hap" not in transpose:
                perhap_out = perhap_out.T
            perhap_out.to_csv(f"{out_dir}/{block_fn}.raw",
                              sep=sep, index=False)
        else:
            print("No haplotypes to include in .per_hap file.")

    if per_loc:
        if perloc_out:
            perloc_out = pd.concat(perloc_out,
                                   axis=1).set_axis(acc_chromatids)
            if "per-loc" not in transpose:
                perloc_out = perloc_out.T
            perloc_out.to_csv(f"{out_dir}/{block_fn}.per_loc",
                              sep=sep)
        else:
            print("No haplotypes to include in .per_loc file.")

    in_gt.close()
    block_info.close()
    hap_info.close()
    hapmap.close()
    locmap.close()

##############################################################################
### MAIN FUNCTION
##############################################################################

def main():

###
### 'block' pipeline
###

    if args.command == "block":

        for file in gt_files:
            print(f"Constructing blocks from {file}")

            if args.type == "gab":
                run([plink_path,
                     "--bfile", os.path.join(args.gt_dir, file),
                     "--allow-extra-chr",
                     "--blocks", "no-pheno-req",
                     "--blocks-max-kb", str(args.gab_max_kb),
                     "--out", os.path.join(args.out, file)],
                    check=True)
                out_fn = os.path.join(args.out, f"{file}.blocks")

            elif args.type == "4gam":
                out_fn = fourgam_blocking(f"{args.gt_dir}/{file}",
                                          f"{args.out}/{file}")

            elif args.type == "bpwin":
                out_fn = win_blocking("bpwin",
                                      f"{args.gt_dir}/{file}",
                                      f"{args.out}/{file}",
                                      args.bp_span,
                                      args.bp_step)

            elif args.type == "snpwin":
                out_fn = win_blocking("snpwin",
                                      f"{args.gt_dir}/{file}",
                                      f"{args.out}/{file}",
                                      args.snp_span,
                                      args.snp_step)

            print(f"Blocks saved in {out_fn}")

###
### 'call' pipeline
###

    elif args.command == "call":
        for (gt_file, block_file) in file_pairs:
            print(f"Genotyping accessions in {gt_file} using {block_file}")

            block_genotyping(f"{args.gt_dir}/{gt_file}",
                             f"{args.block_dir}/{block_file}",
                             args.out,
                             args.format,
                             args.transpose,
                             args.maf,
                             args.min_hap,
                             args.group_minor,
                             args.include_singles,
                             args.tab)

            print(f"Data and info files saved in '{args.out}'.")


if __name__ == "__main__":

    main()