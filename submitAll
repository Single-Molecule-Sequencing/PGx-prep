#!/usr/bin/env python3
# submitAll.py
# Batch submit sbatch files from a directory

import os
import sys
import glob
import time
import argparse
import subprocess


############
# argparse #
############

parser = argparse.ArgumentParser(description="Batch submit sbatch files from a directory")
parser.add_argument("-i", "--input", default="./Sbatch",
	help="Path to the directory containing .sbatch files [%(default)s]")
args = parser.parse_args()


########
# main #
########

sbatch_dir = os.path.abspath(os.path.expanduser(args.input))

if not os.path.exists(sbatch_dir):
	sys.exit(f"[submitAll] Error: Directory {sbatch_dir} does not exist")

# Find files
search_path = os.path.join(sbatch_dir, "*.sbatch")
sbatch_files = glob.glob(search_path)
num_files = len(sbatch_files)

if num_files == 0:
	sys.exit(f"[submitAll] Error: No .sbatch files found in {sbatch_dir}")

print(f"[submitAll] Found {num_files} jobs in {sbatch_dir}. Submitting now...")

# --- Submission Loop ---
c = 0

for file_path in sorted(sbatch_files):
	subprocess.run(["sbatch", file_path])
	print(f"[submitAll] Submitted script: {os.path.basename(file_path)}")
	time.sleep(0.05)
	c += 1

print(f"[submitAll] All {c} jobs submitted")
