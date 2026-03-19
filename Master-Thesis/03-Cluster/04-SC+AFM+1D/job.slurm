#!/bin/bash

#SBATCH --job-name="04-SC+AFM+1D"
#SBATCH --account=def-vcrepel

#SBATCH --time=03:00:00
#SBATCH --mem=8000M

#SBATCH --array=0-14


# Values of g (later multiplied by W)
VALUES_G=(0 0.05 0.15)

# Values of U (Later multiplied by W)
VALUES_U=(0 0.2 0.8 1.2 2.0)

# Lengths
LEN_G=${#VALUES_G[@]}
LEN_U=${#VALUES_U[@]}

# Pick the value based on the current Task ID
SELECT_G=${VALUES_G[$SLURM_ARRAY_TASK_ID / ${LEN_U}]}
SELECT_U=${VALUES_U[$SLURM_ARRAY_TASK_ID % ${LEN_U}]}


# Run the code
module load python
source /home/gbotelho/projects/def-vcrepel/gbotelho/gbvenv/bin/activate
python code.py ${SELECT_G} ${SELECT_U}