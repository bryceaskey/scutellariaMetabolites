#!/bin/bash
#SBATCH --job-name=wgs-blastUnmapped-HAV       # Job name
#SBATCH --account=lee                          # Account name
#SBATCH --qos=lee                              # QOS name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu            # Where to send mail	
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --mem=24gb                              # Job memory request
#SBATCH --time=120:00:00                       # Time limit hrs:min:sec
#SBATCH --output=wgs-blastUnmapped-HAV_%j.log  # Standard output and error log

pwd; hostname; date

module load ncbi_blast/2.9.0

echo "BLAST searching unmapped contigs"

contigs=/ufrc/lee/braskey/Data/WGS/bwa/contigs/
db=/ufrc/data/reference/blast/202006/
output=/ufrc/lee/braskey/Data/WGS/bwa/BLAST_results/

blastn -query ${contigs}HAV/HAV_unmapped-6.fa -db ${db}nt -out ${output}HAV_unmappedBLAST.out -outfmt 6
