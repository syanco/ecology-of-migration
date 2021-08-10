#!/bin/bash

#SBATCH -t 4- #4 days
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH -c 8
#SBATCH -J flow_ndvi_mod

module load miniconda
conda activate parallel_r

scriptsP=~/project/

Rscript $scriptsP/hpc_fit_ndvi.R