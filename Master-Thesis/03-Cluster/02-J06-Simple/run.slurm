#!/bin/bash

#SBATCH --job-name="01-J07-001"
#SBATCH --account=def-vcrepel

#SBATCH --time=12:00:00
#SBATCH --mem=8000M

module load python
source /home/gbotelho/projects/def-vcrepel/gbotelho/gbvenv/bin/activate
python code-0.py out-0.py
python code-5.py out-5.py
python code-15.py out-15.py