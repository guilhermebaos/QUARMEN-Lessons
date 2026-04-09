#!/bin/bash

#SBATCH --job-name="06-MMHK+SC+1D-Par"
#SBATCH --account=def-vcrepel

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=0
#SBATCH --exclusive


# Load Julia
module load julia

# Path to the project files
export JULIA_PROJECT="/project/6108152/gbotelho/"

# Set number of threads
# export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the code
julia speedtest.jl