#!/bin/bash

#SBATCH --job-name="02-Speed-Test"
#SBATCH --account=rrg-vcrepel

#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=0

#SBATCH --exclusive


# Load Julia
module load julia

# Path to the project files (.toml files)
export JULIA_PROJECT="/project/6114430/gbotelho/"

# Run the code
julia speedtest.jl