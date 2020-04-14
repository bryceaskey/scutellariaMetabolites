#!/bin/bash
#SBATCH --job-name=get-scutellaria-data # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu     # Where to send mail	
#SBATCH --ntasks=1                      # Run on a single CPU
#SBATCH --mem=1gb                       # Job memory request
#SBATCH --time=12:00:00                 # Time limit hrs:min:sec
#SBATCH --output=download_data_%j.log   # Standard output and error log
pwd; hostname; date

module load sra/2.10.4

echo "Downloading scutellaria RNA seq data"

# SRP179837
fastq-dump SRR8449931
fastq-dump SRR8449932
fastq-dump SRR8449933
fastq-dump SRR8449934
fastq-dump SRR8449935
fastq-dump SRR8449936
fastq-dump SRR8449937
fastq-dump SRR8449938
fastq-dump SRR8449939

mv SRR8449931.fastq SRP179837/
mv SRR8449932.fastq SRP179837/
mv SRR8449933.fastq SRP179837/
mv SRR8449934.fastq SRP179837/
mv SRR8449935.fastq SRP179837/
mv SRR8449936.fastq SRP179837/
mv SRR8449937.fastq SRP179837/
mv SRR8449938.fastq SRP179837/
mv SRR8449939.fastq SRP179837/