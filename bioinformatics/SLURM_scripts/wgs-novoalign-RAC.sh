#!/bin/bash
#SBATCH --job-name=wgs-novoalign-RAC             # Job name
#SBATCH --account=lee                            # Account name
#SBATCH --qos=lee                                # QOS name
#SBATCH --mail-type=END,FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu              # Where to send mail	
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --mem=4gb                                # Job memory request
#SBATCH --time=120:00:00                         # Time limit hrs:min:sec
#SBATCH --output=wgs-novoalign-RAC_%j.log        # Standard output and error log

pwd; hostname; date

module load novoalign/3.00.02 samtools/1.10

echo "Aligning racemosa WGS data to reference genome"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/
wgs=/ufrc/lee/braskey/Data/WGS/

gunzip -c ${wgs}RAC_1.fq.gz > ${wgs}RAC_1.fq
gunzip -c ${wgs}RAC_2.fq.gz > ${wgs}RAC_2.fq

novoindex ${ref}GCA005771605.nix ${ref}GCA005771605.fa
novoalign -d ${ref}GCA005771605.nix -f  ${wgs}RAC_1.fq ${wgs}RAC_2.fq -o SAM > ${wgs}RAC_novoaln.sam
samtools view -bS ${wgs}RAC_novoaln.sam > ${wgs}RAC_novoaln.bam
samtools sort -O bam -o ${wgs}RAC_novoaln_sorted.bam ${wgs}RAC_novoaln.bam