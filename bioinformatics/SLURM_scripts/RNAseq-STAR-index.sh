#!/bin/bash
#SBATCH --job-name=RNAseq-STAR-index           # Job name
#SBATCH --account=lee                          # Account name
#SBATCH --qos=lee                              # QOS name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu            # Where to send mail	
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --mem=16gb                             # Job memory request
#SBATCH --time=120:00:00                       # Time limit hrs:min:sec
#SBATCH --output=RNAseq-STAR-index_%j.log      # Standard output and error log

pwd; hostname; date

module load star/2.7.3a

echo "Generating genome index for use during STAR mapping"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_STAR/GCA005771605.fa
index=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_STAR/index/

STAR --runThreadN 1 --runMode genomeGenerate --genomeDir ${index} --genomeFastaFiles ${ref} --genomeSAindexNbases 13