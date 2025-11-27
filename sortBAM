#!/usr/bin/env python3
# sortBAM.py
# Generate samtools sort commands for all BAM files in a directory

import argparse
import os
import sys


############
# argparse #
############

parser = argparse.ArgumentParser(
	description="Generate samtools sort commands for all BAM files in a directory")
parser.add_argument("-d", "--bamdir", default=".",
	help="Directory containing BAM files [%(default)s]")
parser.add_argument("-o", "--outdir", default="./Sorted",
	help="Directory to store sorted BAMs [%(default)s]")
parser.add_argument("-f", "--outfile", default="sortCMDs.txt",
	help="Output file to write sort commands [%(default)s]")
args = parser.parse_args()


########
# main #
########

if not os.path.exists(args.bamdir):
	sys.exit(f"[sortBAM] Error: Directory {args.bamdir} does not exist")

os.makedirs(args.outdir, exist_ok=True)

c = 0

with open(args.outfile, "w") as out:
	for fname in os.listdir(args.bamdir):
		if fname.endswith(".bam"):
			bam_path = os.path.join(args.bamdir, fname)
			base = fname.replace(".bam", "")
			sorted_bam = os.path.join(args.outdir, f"{base}.sorted.bam")
			sort_cmd = f"samtools sort -@ 4 -o {sorted_bam} {bam_path}"
			out.write(sort_cmd + "\n")
			c += 1

print(f"[sortBAM] Generated sort commands for {c} BAM files in {args.outfile}")
