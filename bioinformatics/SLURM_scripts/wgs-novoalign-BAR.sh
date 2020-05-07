#!/bin/bash
#SBATCH --job-name=wgs-novoalign-BAR             # Job name
#SBATCH --account=lee                            # Account name
#SBATCH --qos=lee                                # QOS name
#SBATCH --mail-type=END,FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu              # Where to send mail	
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --mem=4gb                                # Job memory request
#SBATCH --time=120:00:00                         # Time limit hrs:min:sec
#SBATCH --output=wgs-novoalign-BAR_%j.log        # Standard output and error log

pwd; hostname; date

module load novoalign/3.00.02 samtools/1.10

echo "Aligning barbata WGS data to reference genome"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_BAR/
wgs=/ufrc/lee/braskey/Data/WGS/

#gunzip -c ${wgs}BAR_1.fq.gz > ${wgs}BAR_1.fq
#gunzip -c ${wgs}BAR_2.fq.gz > ${wgs}BAR_2.fq

novoindex ${ref}GCA005771605.nix ${ref}GCA005771605.fa
novoalign -d ${ref}GCA005771605.nix -f  ${wgs}BAR_1.fq ${wgs}BAR_2.fq -o SAM > ${wgs}BAR_novoaln.sam
samtools view -bS ${wgs}BAR_novoaln.sam > ${wgs}BAR_novoaln.bam
samtools sort -O bam -o ${wgs}BAR_novoaln_sorted.bam ${wgs}BAR_novoaln.bam