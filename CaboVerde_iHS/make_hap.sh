#!/bin/bash
#SBATCH -a 1-22

c=$SLURM_ARRAY_TASK_ID

infile=merged_auto_chr"$c".phased.CVsubset.haps
outfile=chr"$c".hap

cut -d" " --complement -f1-5 $infile > $outfile

