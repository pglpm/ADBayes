#!/bin/bash

#SBATCH --account=nn0000k
#SBATCH --job-name=MCMC3
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=1G
#SBATCH --time=40:00:00

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

## module restore
module --quiet purge
module load R/4.1.0-foss-2021a
## module load intel/2020b

cd ~/ADBayes

Rscript new_posterior_probabilities3.R > new_posterior_probabilities3.Rout 2>&1
