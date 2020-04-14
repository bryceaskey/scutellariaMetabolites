#!/bin/bash
#SBATCH --job-name=get-scutellaria-data # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu     # Where to send mail	
#SBATCH --ntasks=1                      # Run on a single CPU
#SBATCH --mem=1gb                       # Job memory request
#SBATCH --time=24:00:00                 # Time limit hrs:min:sec
#SBATCH --output=download_data_%j.log   # Standard output and error log
pwd; hostname; date

module load sra/2.10.4

echo "Downloading scutellaria RNA seq data"

# SRP068883
fastq-dump SRR3130396
fastq-dump SRR3123399

# SRP156996
fastq-dump SRR7689119
fastq-dump SRR7689120
fastq-dump SRR7689121
fastq-dump SRR7689122
fastq-dump SRR7689123
fastq-dump SRR7665600

# SRP048707
fastq-dump SRR1605127
fastq-dump SRR3367955
fastq-dump SRR3367956
fastq-dump SRR3367957

# SRP096180
fastq-dump SRR6940088 # WGS
fastq-dump SRR5150730 # RNA-seq

# SRP068819
fastq-dump SRR3727161
fastq-dump SRR3114960
