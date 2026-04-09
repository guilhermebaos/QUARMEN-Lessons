#!/bin/bash

#SBATCH --job-name="01-Test-Run"
#SBATCH --account=rrg-vcrepel

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G

#SBATCH --array=0-49


# Sweep a single variable linearly
varMin=0.01
varMax=1.99

# Bash cannot do floating-point math natively, so we pipe the equation to 'bc'
varVal=$(echo "$varMin + ($varMax - $varMin) * ($SLURM_ARRAY_TASK_ID / ($SLURM_ARRAY_TASK_COUNT - 1))" | bc -l)


# Variables to run the code with
runID=$SLURM_ARRAY_TASK_ID
nMMHK=2
L=16000
mu=1
U=0
g=0
T=0
nTarget=$varVal


# Load Julia
module load julia

# Switch to new project!
# Path to the project files (.toml files)
export JULIA_PROJECT="/project/6108152/gbotelho/"

# Set number of threads
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the code
julia code.jl $runID $nMMHK $L $mu $U $g $T $nTarget