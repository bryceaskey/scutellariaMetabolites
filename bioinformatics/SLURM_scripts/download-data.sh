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

dest=/ufrc/lee/braskey/Data/

# SRP179837
fastq-dump --split-files SRR8449931 -O ${dest}SRP179837/
fastq-dump --split-files SRR8449932 -O ${dest}SRP179837/
fastq-dump --split-files SRR8449933 -O ${dest}SRP179837/
fastq-dump --split-files SRR8449934 -O ${dest}SRP179837/
fastq-dump --split-files SRR8449935 -O ${dest}SRP179837/
fastq-dump --split-files SRR8449936 -O ${dest}SRP179837/
fastq-dump --split-files SRR8449937 -O ${dest}SRP179837/
fastq-dump --split-files SRR8449938 -O ${dest}SRP179837/
fastq-dump --split-files SRR8449939 -O ${dest}SRP179837/

# SRP068883
fastq-dump SRR3130396 -O ${dest}SRP068883/
fastq-dump SRR3123399 -O ${dest}SRP068883/

# SRP156996
fastq-dump SRR7689119 -O ${dest}SRP156996/
fastq-dump SRR7689120 -O ${dest}SRP156996/
fastq-dump SRR7689121 -O ${dest}SRP156996/
fastq-dump SRR7689122 -O ${dest}SRP156996/
fastq-dump SRR7689123 -O ${dest}SRP156996/
fastq-dump SRR7665600 -O ${dest}SRP156996/

# SRP048707
fastq-dump SRR1605127 -O ${dest}SRP048707/
fastq-dump SRR3367955 -O ${dest}SRP048707/
fastq-dump SRR3367956 -O ${dest}SRP048707/
fastq-dump SRR3367957 -O ${dest}SRP048707/

# SRP096180
fastq-dump --split-files SRR6940088 -O ${dest}SRP096180/ # WGS
fastq-dump SRR5150730 -O ${dest}SRP096180/ # RNA-seq

# SRP068819
fastq-dump SRR3727161 -O ${dest}SRP068819/
fastq-dump SRR3114960 -O ${dest}SRP068819/