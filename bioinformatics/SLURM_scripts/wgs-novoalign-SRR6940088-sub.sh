#!/bin/bash
#SBATCH --job-name=wgs-novoalign-SRR6940088-sub         # Job name
#SBATCH --account=lee                                   # Account name
#SBATCH --qos=lee                                       # QOS name
#SBATCH --mail-type=END,FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                     # Where to send mail	
#SBATCH --ntasks=1                                      # Run on a single CPU
#SBATCH --mem=4gb                                       # Job memory request
#SBATCH --time=744:00:00                                # Time limit hrs:min:sec
#SBATCH --output=wgs-novoalign-SRR6940088-sub_%j.log    # Standard output and error log

pwd; hostname; date

module load novoalign/3.00.02 samtools/1.10 seqtk/1.3

echo "Aligning a subset of baicalensis WGS data to reference genome"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_ALT/
wgs=/ufrc/lee/braskey/Data/SRP096180/

seqtk sample -s100 ${wgs}SRR6940088_1.fastq 100000 > ${wgs}SRR6940088_sub1.fastq
seqtk sample -s100 ${wgs}SRR6940088_2.fastq 100000 > ${wgs}SRR6940088_sub2.fastq

novoindex ${ref}GCA005771605.nix ${ref}GCA005771605.fa
novoalign -d ${ref}GCA005771605.nix -f  ${wgs}SRR6940088_sub1.fastq ${wgs}SRR6940088_sub2.fastq -o SAM > ${wgs}SRR6940088_sub_novoaln.sam
samtools view -bS ${wgs}SRR6940088_sub_novoaln.sam > ${wgs}SRR6940088_sub_novoaln.bam
samtools sort -O bam -o ${wgs}SRR6940088_sub_novoaln_sorted.bam ${wgs}SRR6940088_sub_novoaln.bam
samtools index ${wgs}SRR6940088_sub_novoaln_sorted.bam
samtools stats ${wgs}SRR6940088_sub_novoaln_sorted.bam