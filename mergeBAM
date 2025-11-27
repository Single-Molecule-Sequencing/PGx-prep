#!/usr/bin/env python3
# mergeBAM.py
# Generate samtools merge commands per sample from a SampleSheet CSV file

import argparse
import csv
from pathlib import Path
import sys


#########
# funcs #
#########

def find_bams(bamdir: Path, source_prefixes: list[str], barcode: str) -> list[str]:
	barcode_pattern = f"barcode{int(barcode):02d}"
	matches = []

	# Iterate through files using pathlib
	for fpath in bamdir.iterdir():
		if fpath.suffix == ".bam" and barcode_pattern in fpath.name:
			# Check if filename contains one of the required source prefixes
			for prefix in source_prefixes:
				if prefix in fpath.name:
					# fpath is an absolute Path object here because bamdir is absolute
					matches.append(str(fpath))
					break

	return sorted(matches)


############
# argparse #
############

parser = argparse.ArgumentParser(
	description="Generate samtools merge commands per sample from a SampleSheet CSV file")
parser.add_argument("csv",
	help="Path to the SampleSheet CSV file")
parser.add_argument("-d", "--bamdir", default=".",
	help="Directory containing sorted&indexed BAM files [%(default)s]")
parser.add_argument("-o", "--outdir", default="./Merged",
	help="Directory to store merged BAMs [%(default)s]")
parser.add_argument("-f", "--outfile", default="mergeCMDs.txt",
	help="Output file to write merge commands [%(default)s]")
args = parser.parse_args()


########
# main #
########

if not Path(args.csv).exists():
	sys.exit(f"[mergeBAM] Error: CSV file {args.csv} does not exist")

# Resolve directories to absolute paths
bam_dir = Path(args.bamdir).resolve()
out_dir = Path(args.outdir).resolve()

# Create output directory
out_dir.mkdir(parents=True, exist_ok=True)

# Read CSV and generate commands
c = 0

with open(args.csv, newline='', encoding='utf-8-sig') as csvfile, open(args.outfile, "w") as out:
	reader = csv.DictReader(csvfile)

	for row in reader:
		pid = row["Patient_ID"]
		experiment = row["Experiment_ID"]
		barcode = row["Barcode_Number"]

		# Get list of valid filename components
		source_prefixes = row["Source_Folders"].split(';')

		# Find matching files
		bam_files = find_bams(bam_dir, source_prefixes, barcode)

		# Generate command
		bam_list = " ".join(bam_files)

		# Construct absolute path for the output file
		output_bam = out_dir / f"{pid}-{experiment}-{barcode}.bam"

		cmd = f"samtools merge -o {output_bam} {bam_list}"
		out.write(cmd + "\n")
		c += 1

print(f"[mergeBAM] Generated merge commands for {c} BAM files in {args.outfile}")
