#!/usr/bin/env python3
# cmdtxt2sbatch.py
# Generate Slurm .sbatch files from a list of shell commands in a text file

import argparse
from pathlib import Path
import shlex
import sys


#########
# funcs #
#########

def derive_job_name(cmd: str, job_prefix: str, fallback_idx: int) -> str:

	TRIM_EXTS = {'bam', 'sam', 'cram', 'bai', 'csi', 'tbi', 'sorted', 'gz', 'bed', 'bw'}

	parts = cmd.split()
	base = None

	if "-o" in parts:
		flag_index = parts.index("-o")
		if flag_index + 1 < len(parts):
			output_file = parts[flag_index + 1]
			base = output_file.rsplit("/", 1)[-1]

	if not base:
		base = cmd.rsplit("/", 1)[-1].split("_")[-1]

	while True:
		split_parts = base.rsplit('.', 1)
		if len(split_parts) > 1 and split_parts[-1] in TRIM_EXTS:
			base = split_parts[0]
		else:
			break

	# base = base.replace(" ", "_").replace("-", "_")

	if base.endswith("_R1") or base.endswith("_R2"):
		base = base[:-3]

	return f"{job_prefix}_{base}" if base else f"{job_prefix}_{fallback_idx:04d}"


def build_header(partition: str, job_name: str, account: str, gres: str, cpus: int, mem: str,
				 walltime: str, email: str, logs_dir: Path, module: str) -> str:
	lines = [
		"#!/bin/bash",
		f"#SBATCH --partition={partition}",
		f"#SBATCH --account={account}",
		f"#SBATCH --job-name={job_name}",
		f"#SBATCH --gres={gres}" if gres else "",
		f"#SBATCH --nodes=1",
		f"#SBATCH --ntasks=1",
		f"#SBATCH --cpus-per-task={cpus}",
		f"#SBATCH --mem={mem}",
		f"#SBATCH --time={walltime}",
		f"#SBATCH --output={logs_dir}/{job_name}_%j.out",
		f"#SBATCH --error={logs_dir}/{job_name}_%j.err",
	]
	if email:
		lines += [f"#SBATCH --mail-user={email}",
				  "#SBATCH --mail-type=BEGIN,END,FAIL"]
	lines += [
		"",
		"set -euo pipefail",
		'echo "[$(date)] Node: $(hostname)"',
		'echo "[$(date)] CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"' if gres else "",
	]
	if module:
		lines += ["", f"module load {module}"]
	lines += ["", "# run command"]
	return "\n".join(lines)


############
# argparse #
############

parser = argparse.ArgumentParser(description="Generate Slurm .sbatch files from a list of shell commands in a text file.")
parser.add_argument("-i", "--input", required=True,
	help="Path to input text file with one command per line")
parser.add_argument("-o", "--outdir", default=".",
	help="Directory to write .sbatch files [%(default)s]")
parser.add_argument("-l", "--logs", default=None,
	help="Directory for Slurm logs [None]")
parser.add_argument("-p", "--partition", default="standard",
	help="Slurm partition [%(default)s]")
parser.add_argument("-a", "--account", required=True,
	help="Slurm account")
parser.add_argument("-g", "--gres", default="",
	help="Slurm GRES string (e.g., 'gpu:gpuname:1') [%(default)s]")
parser.add_argument("-@", "-c", "--cpus", type=int, default=8,
	help="CPUs per task [%(default)i]; note: '-@' may need quoting in some shells")
parser.add_argument("--mem", default="32G",
	help="Memory per task [%(default)s]")
parser.add_argument("--time", default="12:00:00",
	help="Walltime (HH:MM:SS) [%(default)s]")
parser.add_argument("--email",
	help="Email for Slurm notifications")
parser.add_argument("--job-prefix", default="job",
	help="Prefix for job names [%(default)s]")
parser.add_argument("--module", default="",
	help="Optional module to load")
args = parser.parse_args()


########
# main #
########

inpath = Path(args.input)
if not inpath.exists():
	sys.exit(f"[cmdtxt2sbatch] Input file not found: {inpath}")

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

logs_dir = Path(args.logs) if args.logs else Path("./Logs")
logs_dir.mkdir(parents=True, exist_ok=True)

cmds = [ln.strip() for ln in inpath.read_text(encoding="utf-8").splitlines() if ln.strip()]
if not cmds:
	sys.exit(f"[cmdtxt2sbatch] No commands found in input file: {inpath}")

written = 0
for idx, cmd in enumerate(cmds, start=1):
	job_name = derive_job_name(cmd, args.job_prefix, idx)
	header = build_header(
		partition = args.partition,
		account   = args.account,
		job_name  = job_name,
		gres      = args.gres,
		cpus      = args.cpus,
		mem       = args.mem,
		walltime  = args.time,
		email     = args.email,
		logs_dir  = logs_dir,
		module    = args.module.strip()
	)

	cmd_tokens = shlex.split(cmd)
	grouped_tokens = []
	i = 0
	while i < len(cmd_tokens):
		if cmd_tokens[i].startswith("-") and i + 1 < len(cmd_tokens) and not cmd_tokens[i + 1].startswith("-"):
			grouped_tokens.append(f"{cmd_tokens[i]} {cmd_tokens[i + 1]}")
			i += 2
		else:
			grouped_tokens.append(cmd_tokens[i])
			i += 1
	cmd_multiline = " \\\n  ".join(grouped_tokens)

	script = f"{header}\n{cmd_multiline}\n"

	outpath = outdir / f"{job_name}.sbatch"
	if outpath.exists():
		outpath = outdir / f"{job_name}_{idx:03d}.sbatch"
	outpath.write_text(script, encoding="utf-8")
	written += 1

print(f"[cmdtxt2sbatch] Wrote {written} .sbatch files to {outdir}/")
print(f"[cmdtxt2sbatch] Logs will go to {logs_dir}/")
