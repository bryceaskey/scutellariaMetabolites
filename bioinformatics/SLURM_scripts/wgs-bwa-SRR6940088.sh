#!/bin/bash
#SBATCH --job-name=bwa-SRR6940088           # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=120:00:00                    # Time limit hrs:min:sec
#SBATCH --output=bwa-SRR6940088_%j.log      # Standard output and error log

pwd; hostname; date

module load bwa/0.7.17 samtools/1.10

echo "Aligning SRR6940088 WGS shotgun sequencing data to reference genome"

index=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/
reads=/ufrc/lee/braskey/Data/SRP096180/
aln=/ufrc/lee/braskey/Data/SRP096180/bwa/

#bwa index ${ref}GCA005771605.fa
bwa mem ${index}GCA005771605.fa ${reads}SRR6940088_1.fastq ${reads}SRR6940088_2.fastq > ${aln}SRR6940088.sam
samtools fixmate -O bam ${aln}SRR6940088.sam ${aln}SRR6940088_fixmate.bam
samtools sort -O bam -o ${aln}SRR6940088_sorted.bam ${aln}SRR6940088_fixmate.bam
samtools stats ${aln}SRR6940088_sorted.bam