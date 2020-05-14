#!/bin/bash
#SBATCH --job-name=variant-vis-HAV          # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=16gb                          # Job memory request
#SBATCH --time=120:00:00                    # Time limit hrs:min:sec
#SBATCH --output=variant-vis-HAV_%j.log     # Standard output and error log

pwd; hostname; date

module load R/3.6

echo "Generating visualization of variants for havanesis"

Rscript /ufrc/lee/braskey/R_scripts/variant-vis-HAV.R