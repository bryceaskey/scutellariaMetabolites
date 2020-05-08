#!/bin/bash
#SBATCH --job-name=bwa-test-HAV-stats       # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=120:00:00                    # Time limit hrs:min:sec
#SBATCH --output=bwa-test-HAV-stats_%j.log  # Standard output and error log

pwd; hostname; date

module load samtools/1.10 bcftools/1.10.2

wgs=/ufrc/lee/braskey/Data/WGS/bwa_test/

samtools stats ${wgs}HAV_aln_sorted.bam
bcftools stats ${wgs}HAV_aln_filtered.vcf