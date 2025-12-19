# PGx-prep

Preparatory post-basecalling demultiplex algorithms and HPC+Slurm solutions for processing BAM files ahead of the ONT's PGx workflow.

This repository provides a suite of Python scripts to automate the sorting, merging, and indexing of BAM files, tailored for high-performance computing (HPC) environments using the Slurm workload manager.

## Manifest

-   `sortBAM`: Generates `samtools sort` commands for all BAM files in a given directory.
-   `mergeBAM`: Generates `samtools merge` commands for grouping BAM files by sample, based on a `SampleSheet.csv` file.
-   `indexBAM`: Indexes all BAM files in a directory using parallel processing.
-   `cmdtxt2sbatch`: Converts a text file of shell commands into individual Slurm `.sbatch` scripts for job submission.
-   `submitAll`: Submits all `.sbatch` files in a specified directory to the Slurm scheduler.

## Usage

### `sortBAM`

Generates `samtools sort` commands.

```bash
sortBAM -d <bam_directory> -o <output_directory> -f <output_file>
```

-   `-d, --bamdir`: Directory containing input BAM files (default: `.`).
-   `-o, --outdir`: Directory to store sorted BAMs (default: `./Sorted`).
-   `-f, --outfile`: Output file to write sort commands (default: `sortCMDs.txt`).

### `mergeBAM`

Generates `samtools merge` commands from a SampleSheet obtained at [SSS](https://single-molecule-sequencing.github.io/sss/).

```bash
mergeBAM <samplesheet_csv> -d <bam_directory> -o <output_directory> -f <output_file>
```

-   `csv`: Path to the SampleSheet CSV file (required).
-   `-d, --bamdir`: Directory containing sorted & indexed BAM files (default: `.`).
-   `-o, --outdir`: Directory to store merged BAMs (default: `./Merged`).
-   `-f, --outfile`: Output file to write merge commands (default: `mergeCMDs.txt`).

### `indexBAM`

Indexes BAM files in parallel.

```bash
indexBAM -i <bam_directory> -j <num_jobs>
```

-   `-i, --input`: Path to the directory containing BAM files (required).
-   `-j, --jobs`: Number of parallel jobs to run (default: 4).

### `cmdtxt2sbatch`

Generates `.sbatch` files from a command list.

```bash
cmdtxt2sbatch -i <input_file> -a <slurm_account> [options]
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

### `submitAll`

Submits all `.sbatch` files in a directory.

```bash
submitAll -i <sbatch_directory>
```

-   `-i, --input`: Directory containing `.sbatch` files (default: `./Sbatch`).

## Workflow

Make sure to set up and activate the conda environment using the provided `env.yml`.

```
conda env create -f env.yml
conda activate pgxprep
```

The workflow expects a directory of raw BAM files as the initial input.

```
Input: Directory of BAM files (e.g., 'BAMs/')
  │
  ▼
┌────────────────┐
│  1. Sort BAMs  │
└────────────────┘
  │
  │ Input: 'BAMs/'
  │ Output: 'Sorted/' directory with sorted BAMs
  │
  ├─► sortBAM -d BAMs -o Sorted -f sortCMDs.txt
  │
  ├─► cmdtxt2sbatch -i sortCMDs.txt -o Sbatch/sort --account bleu99 --partition standard \
  │     --cpus 4 --mem 16G --time 03:00:00 --email theatheylab@gmail.com --job-prefix sort
  │
  └─► submitAll -i Sbatch/sort
  │
  ▼
┌─────────────────┐
│  2. Merge BAMs  │
└─────────────────┘
  │
  │ Input: 'Sorted/' directory, 'SampleSheet.csv'
  │ Output: 'Merged/' directory with merged BAMs
  │
  ├─► mergeBAM SampleSheet.csv -d Sorted -o Merged -f mergeCMDs.txt
  │
  ├─► cmdtxt2sbatch -i mergeCMDs.txt -o Sbatch/merge --account bleu99 --partition standard \
  │     --cpus 2 --mem 8G --time 03:00:00 --email theatheylab@gmail.com --job-prefix merge
  │
  └─► submitAll -i Sbatch/merge
  │
  ▼
┌─────────────────┐
│  3. Index BAMs  │
└─────────────────┘
  │
  │ Input: 'Sorted/'
  │ Output: '.bai' index files in 'Sorted/'
  │
  └─► nohup indexBAM -i Sorted -j 8 > index.log 2>&1 &
  │
  ▼
Final Output: Directory of sorted, merged, and indexed BAM files ready for the PGx workflow.
```

### Step 1: Sort BAM Files

First, generate the `samtools sort` commands. This example assumes the raw BAMs are in a directory named `BAMs/`.

```bash
sortBAM -d BAMs -o Sorted -f sortCMDs.txt
```

Next, convert the commands into Slurm jobs. The following command creates `.sbatch` files in `Sbatch/sort/`.

```bash
cmdtxt2sbatch \
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
submitAll -i Sbatch/sort
```

### Step 2: Merge BAM Files by Sample

Generate merge commands using a `SampleSheet.csv` file.

```bash
mergeBAM SampleSheet.csv -d Sorted -o Merged -f mergeCMDs.txt
```

Convert the merge commands into Slurm jobs, writing the `.sbatch` files to `Sbatch/merge/`.

```bash
cmdtxt2sbatch \
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
submitAll -i Sbatch/merge
```

### Step 3: Index Merged BAM Files

After the merging jobs are complete, index the merged BAM files in the `Merged/` directory. It is recommended to run this process with `nohup` to prevent interruptions.

```bash
nohup indexBAM -i Merged -j 8 > index.log 2>&1 &
```

## Future Plans

-   A master script to execute the entire workflow sequentially.
-   Potentially turning this into a NextFlow workflow.
-   Pray to the Automation God that this all will work properly.
