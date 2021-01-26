#!/bin/bash
#SBATCH -p scavenger
#SBATCH -a 1-22

module load OpenMPI/4.0.1

c=$SLURM_ARRAY_TASK_ID

~/hapbin/build/ihsbin --hap ~/work/CV_localancestry/iHS_calculation/NWCluster_chr"$c".hap --map ~/work/CV_localancestry/iHS_calculation/chr"$c".map --out NWCluster_chr"$c"_iHS --minmaf 0.1 --cutoff 0.1
