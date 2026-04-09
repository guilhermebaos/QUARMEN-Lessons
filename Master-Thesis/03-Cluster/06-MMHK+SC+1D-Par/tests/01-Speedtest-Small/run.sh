#!/bin/bash

# Directory to delete
OUT="outputs"

# Check if the directory exists
if [ -d "$OUT" ]; then
    rm -rf "$OUT"
fi


# Remove all old slurm outputs
rm slurm-*


# Configure the slurm file and submit the job
dos2unix job.slurm
sbatch job.slurm