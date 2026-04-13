#!/bin/bash

# Directory to delete
OUT="outputs"

# Clean the directory
rm -rf "$OUT"
mkdir -p "$OUT"

# Remove all old slurm outputs
rm slurm-*

# Configure the slurm file and submit the job
dos2unix job.slurm
sbatch job.slurm

# Configure the process file
dos2unix process.sh
chmod +x process.sh