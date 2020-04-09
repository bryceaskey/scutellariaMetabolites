#!/bin/bash
#SBATCH --job-name=wgs-alignment        # Job name
#SBATCH --account=lee                   # Account name
#SBATCH --qos=lee                       # QOS name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu     # Where to send mail	
#SBATCH --ntasks=1                      # Run on a single CPU
#SBATCH --mem=4gb                       # Job memory request
#SBATCH --time=720                      # Time limit hrs:min:sec
#SBATCH --output=wgs-alignment_%j.log   # Standard output and error log

pwd; hostname; date

module load bwa/0.7.17

echo "Aligning altissima WGS data to reference genome"
# bwa mem reference reads > output
bwa mem /ufrc/lee/braskey/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/... /ufrc/lee/braskey/WGS/ATM_1.fq.gz > /ufrc/lee/braskey/WGS/ATM_1-aln.sam
