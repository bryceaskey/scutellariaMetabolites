#!/bin/bash
#SBATCH --job-name=combine-wgs-data        # Job name
#SBATCH --mail-type=END,FAIL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu        # Where to send mail	
#SBATCH --ntasks=1                         # Run on a single CPU
#SBATCH --mem=4gb                          # Job memory request
#SBATCH --time=12:00:00                    # Time limit hrs:min:sec
#SBATCH --output=combine-wgs-data_%j.log   # Standard output and error log

pwd; hostname; date

module load python/3.8

python /ufrc/lee/braskey/Python_scripts/combine-wgs-data.py