#!/bin/bash
#SBATCH --mem=2G
#SBATCH -p scavenger

~/hapbin/build/ehhbin --hap ~/work/CV_localancestry/iHS_calculation/NWCluster_chr1.hap --map ~/work/CV_localancestry/iHS_calculation/chr1.map --locus rs2814778 > NWCluster_Duffy_EHH.txt
