# PGx-prep

Preparatory algorithms and HPC+Slurm solutions for processing BAM files ahead of the ONT's PGx workflow.

This repository provides a suite of Python scripts to automate the sorting, indexing, and merging of BAM files, tailored for high-performance computing (HPC) environments using the Slurm workload manager.

## Manifest

-   `sortBAM.py`: Generates `samtools sort` commands for all BAM files in a given directory.
-   `indexBAM.py`: Indexes all BAM files in a directory using parallel processing.
-   `mergeBAM.py`: Generates `samtools merge` commands for grouping BAM files by sample, based on a `SampleSheet.csv` file.
-   `cmdtxt2sbatch.py`: Converts a text file of shell commands into individual Slurm `.sbatch` scripts for job submission.
-   `submitAll.py`: Submits all `.sbatch` files in a specified directory to the Slurm scheduler.

## Usage

### `sortBAM.py`

Generates `samtools sort` commands.

```bash
python3 sortBAM.py -d <bam_directory> -o <output_directory> -f <output_file>
```

-   `-d, --bamdir`: Directory containing input BAM files (default: `.`).
-   `-o, --outdir`: Directory to store sorted BAMs (default: `./Sorted`).
-   `-f, --outfile`: Output file to write sort commands (default: `sortCMDs.txt`).

### `indexBAM.py`

Indexes BAM files in parallel.

```bash
python3 indexBAM.py -i <bam_directory> -j <num_jobs>
```

-   `-i, --input`: Path to the directory containing BAM files (required).
-   `-j, --jobs`: Number of parallel jobs to run (default: 4).

### `mergeBAM.py`

Generates `samtools merge` commands from a SampleSheet obtained at [SSS](https://single-molecule-sequencing.github.io/sss/).

```bash
python3 mergeBAM.py <samplesheet_csv> -d <bam_directory> -o <output_directory> -f <output_file>
```

-   `csv`: Path to the SampleSheet CSV file (required).
-   `-d, --bamdir`: Directory containing sorted & indexed BAM files (default: `.`).
-   `-o, --outdir`: Directory to store merged BAMs (default: `./Merged`).
-   `-f, --outfile`: Output file to write merge commands (default: `mergeCMDs.txt`).

### `cmdtxt2sbatch.py`

Generates `.sbatch` files from a command list.

```bash
python3 cmdtxt2sbatch.py -i <input_file> -a <slurm_account> [options]
```

-   `-i, --input`: Path to input text file with one command per line (required).
-   `-a, --account`: Slurm account (required).
-   `-o, --outdir`: Directory to write `.sbatch` files (default: `.`).
-   `-l, --logs`: Directory for Slurm logs (default: `./Logs`).
-   `-p, --partition`: Slurm partition (default: `standard`).
-   `-g, --gres`: Slurm GRES string (e.g., 'gpu:1') (default: "").
-   `-@, --cpus`: CPUs per task (default: 8).
-   `--mem`: Memory per task (default: `32G`).
-   `--time`: Walltime (HH:MM:SS) (default: `12:00:00`).
-   `--email`: Email for Slurm notifications.
-   `--job-prefix`: Prefix for job names (default: `job`).
-   `--module`: Optional module to load.

### `submitAll.py`

Submits all `.sbatch` files in a directory.

```bash
python3 submitAll.py -i <sbatch_directory>
```

-   `-i, --input`: Directory containing `.sbatch` files (default: `./Sbatch`).

## Workflow

The workflow expects a directory of raw BAM files as the initial input.

```
Input: Directory of BAM files (e.g., 'BAMs/')
  │
  ▼
┌───────────────────┐
│   1. Sort BAMs    │
└───────────────────┘
  │
  │ Input: 'BAMs/'
  │ Output: 'Sorted/' directory with sorted BAMs
  │
  ├─► python3 sortBAM.py -d BAMs -o Sorted -f sortCMDs.txt
  │
  ├─► python3 cmdtxt2sbatch.py -i sortCMDs.txt -o Sbatch/sort --account bleu99 --partition standard \
  │     --cpus 4 --mem 16G --time 03:00:00 --email theatheylab@gmail.com --job-prefix sort
  │
  └─► python3 submitAll.py -i Sbatch/sort
  │
  ▼
┌───────────────────┐
│  2. Index BAMs    │
└───────────────────┘
  │
  │ Input: 'Sorted/'
  │ Output: '.bai' index files in 'Sorted/'
  │
  └─► nohup python3 indexBAM.py -i Sorted -j 8 > index.log 2>&1 &
  │
  ▼
┌───────────────────┐
│   3. Merge BAMs   │
└───────────────────┘
  │
  │ Input: 'Sorted/' directory, 'SampleSheet.csv'
  │ Output: 'Merged/' directory with merged BAMs
  │
  ├─► python3 mergeBAM.py SampleSheet.csv -d Sorted -o Merged -f mergeCMDs.txt
  │
  ├─► python3 cmdtxt2sbatch.py -i mergeCMDs.txt -o Sbatch/merge --account bleu99 --partition standard \
  │     --cpus 2 --mem 8G --time 03:00:00 --email theatheylab@gmail.com --job-prefix merge
  │
  └─► python3 submitAll.py -i Sbatch/merge
  │
  ▼
Final Output: Directory of sorted, indexed, and merged BAM files ready for the PGx workflow.
```

### Step 1: Sort BAM Files

First, generate the `samtools sort` commands. This example assumes the raw BAMs are in a directory named `BAMs/`.

```bash
python3 sortBAM.py -d BAMs -o Sorted -f sortCMDs.txt
```

Next, convert the commands into Slurm jobs. The following command creates `.sbatch` files in `Sbatch/sort/`.

```bash
python3 cmdtxt2sbatch.py \
  -i sortCMDs.txt \
  -o Sbatch/sort \
  -a bleu99 \
  -p standard \
  -@ 4 \
  --mem 16G \
  --time "03:00:00" \
  --email "theatheylab@gmail.com" \
  --job-prefix "sort"
```

Finally, submit all sorting jobs to Slurm.

```bash
python3 submitAll.py -i Sbatch/sort
```

### Step 2: Index Sorted BAM Files

After the sorting jobs are complete, index the sorted BAM files in the `Sorted/` directory. It is recommended to run this process with `nohup` to prevent interruptions.

```bash
nohup python3 indexBAM.py -i Sorted -j 8 > index.log 2>&1 &
```

### Step 3: Merge BAM Files by Sample

Generate merge commands using a `SampleSheet.csv` file.

```bash
python3 mergeBAM.py SampleSheet.csv -d Sorted -o Merged -f mergeCMDs.txt
```

Convert the merge commands into Slurm jobs, writing the `.sbatch` files to `Sbatch/merge/`.

```bash
python3 cmdtxt2sbatch.py \
  -i mergeCMDs.txt \
  -o Sbatch/merge \
  -a bleu99 \
  -p standard \
  -@ 2 \
  --mem 8G \
  --time "03:00:00" \
  --email "theatheylab@gmail.com" \
  --job-prefix "merge"
```

Submit the merge jobs.

```bash
python3 submitAll.py -i Sbatch/merge
```

## Future Plans

-   A master script to execute the entire workflow sequentially.
-   Potentially turning this into a NextFlow workflow.
-   Pray to the Automation God that this all will work properly.
