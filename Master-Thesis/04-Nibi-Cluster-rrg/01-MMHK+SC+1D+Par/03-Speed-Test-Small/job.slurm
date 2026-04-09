#!/bin/bash

#SBATCH --job-name="03-Speed-Test-Small"
#SBATCH --account=rrg-vcrepel

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G


# Load Julia
module load julia

# Path to the project files (.toml files)
export JULIA_PROJECT="/project/6114430/gbotelho/"

# Run the code
julia speedtest.jl