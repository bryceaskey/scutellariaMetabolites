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

echo "Aligning SRR6940088 transcriptome data to reference genome"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/
main=/ufrc/lee/braskey/Data/SRP096180/

bwa index ${ref}GCA005771605.fa
bwa mem ${ref}GCA005771605.fa ${main}SRR6940088.fastq > ${main}SRR6940088_aln.sam
samtools fixmate -O bam ${main}SRR6940088_aln.sam ${main}SRR6940088_aln_fixmate.bam
samtools sort -O bam -o ${main}SRR6940088_aln_sorted.bam ${main}SRR6940088_aln_fixmate.bam
samtools stats ${main}SRR6940088_aln_sorted.bam