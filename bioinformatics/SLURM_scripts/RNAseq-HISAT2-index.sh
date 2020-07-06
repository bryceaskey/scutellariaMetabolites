#!/bin/bash
#SBATCH --job-name=RNAseq-HISAT2-index         # Job name
#SBATCH --account=lee                          # Account name
#SBATCH --qos=lee                              # QOS name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu            # Where to send mail	
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --mem=16gb                             # Job memory request
#SBATCH --time=120:00:00                       # Time limit hrs:min:sec
#SBATCH --output=RNAseq-STAR-HISAT2_%j.log     # Standard output and error log

pwd; hostname; date

module load hisat2/2.2.0

echo "Generating index for use during HISAT2 mapping"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/GCA005771605.fa
index=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/HISAT2-index/

# Copy reference fasta into index folder, and generate index
cp ${ref} ${index}
hisat2-build ${index}GCA005771605.fa ${index}GCA005771605