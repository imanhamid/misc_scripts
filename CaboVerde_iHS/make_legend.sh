#!/bin/bash
#SBATCH -a 1-22

c=$SLURM_ARRAY_TASK_ID

infile=merged_auto_chr"$c".phased.CVsubset.haps

outfile=CV_chr"$c".legend

cut -f2,3 -d" " $infile > $outfile
