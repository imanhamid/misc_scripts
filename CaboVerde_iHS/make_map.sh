#!/bin/bash
#SBATCH -a 1-22
#SBATCH --mem=5G

module load R/4.0.0

c=$SLURM_ARRAY_TASK_ID

./make_map.R $c

