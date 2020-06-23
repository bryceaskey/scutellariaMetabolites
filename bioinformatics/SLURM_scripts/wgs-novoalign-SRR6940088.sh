#!/bin/bash
#SBATCH --job-name=wgs-novoalign-SRR6940088             # Job name
#SBATCH --account=lee                                   # Account name
#SBATCH --qos=lee                                       # QOS name
#SBATCH --mail-type=END,FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                     # Where to send mail	
#SBATCH --ntasks=1                                      # Run on a single CPU
#SBATCH --mem=4gb                                       # Job memory request
#SBATCH --time=744:00:00                                # Time limit hrs:min:sec
#SBATCH --output=wgs-novoalign-SRR6940088_%j.log        # Standard output and error log

pwd; hostname; date

module load novoalign/3.00.02 samtools/1.10

echo "Aligning baicalensis WGS data to reference genome"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_SRR6940088/
wgs=/ufrc/lee/braskey/Data/SRP096180/
aln=/ufrc/lee/braskey/Data/SRP096180/novoaln/

#novoindex ${ref}GCA005771605.nix ${ref}GCA005771605.fa
novoalign -d ${ref}GCA005771605.nix -f  ${wgs}SRR6940088_1.fastq ${wgs}SRR6940088_2.fastq -o SAM > ${aln}SRR6940088_novoaln.sam
samtools view -bS ${aln}SRR6940088_novoaln.sam > ${aln}SRR6940088_novoaln.bam
samtools sort -O bam -o ${aln}SRR6940088_novoaln_sorted.bam ${aln}SRR6940088_novoaln.bam
samtools index ${aln}SRR6940088_novoaln_sorted.bam
samtools stats ${aln}SRR6940088_novoaln_sorted.bam