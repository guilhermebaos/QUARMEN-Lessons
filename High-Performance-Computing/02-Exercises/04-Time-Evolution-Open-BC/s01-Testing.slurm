#!/bin/bash

#SBATCH --job-name="KPM-G-P"                    # job name

#SBATCH --account=epor50                        # account name
#SBATCH --qos=gp_debug                          # quality of service / queue

#SBATCH --ntasks=1                              # number of processor cores
#SBATCH --cpus-per-task=112                     # number of threads for each process

#SBATCH --time=00:30:00                         # walltime


./03-Time-Evolution.sh


# chmod +x 03-Time-Evolution.sh ; perl -pi -e 's/\r\n/\n/' 03-Time-Evolution.sh ; perl -pi -e 's/\r\n/\n/' s01*.slurm ; sbatch s01*.slurm
