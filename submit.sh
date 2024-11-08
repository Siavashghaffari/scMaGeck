#!/bin/bash

#SBATCH -p himem
#SBATCH --mem=512G
#SBATCH -c 8


ml R/r421-bioc316-20221122_prd
Rscript EFF_GeneList.R
