#!/bin/bash
#SBATCH --job-name=wgs-getOrgName-RAC          # Job name
#SBATCH --account=lee                          # Account name
#SBATCH --qos=lee                              # QOS name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu            # Where to send mail	
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --mem=4gb                              # Job memory request
#SBATCH --time=744:00:00                       # Time limit hrs:min:sec
#SBATCH --output=wgs-getOrgName-RAC_%j.log     # Standard output and error log

pwd; hostname; date
module load R/4.0
echo "Finding organism names for BLAST search results"

Rscript /blue/lee/braskey/R_scripts/get_sci_name_RAC.R