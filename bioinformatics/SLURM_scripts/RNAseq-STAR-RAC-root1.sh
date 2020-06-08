#!/bin/bash
#SBATCH --job-name=RNAseq-STAR-RAC-root1       # Job name
#SBATCH --account=lee                          # Account name
#SBATCH --qos=lee                              # QOS name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu            # Where to send mail	
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --mem=8gb                              # Job memory request
#SBATCH --time=120:00:00                       # Time limit hrs:min:sec
#SBATCH --output=RNAseq-STAR-RAC-root1_%j.log  # Standard output and error log

pwd; hostname; date

module load star/2.7.3a

echo "Mapping racemosa root1 RNAseq data to baicalensis reference genome"

sp=RAC
rep=root1
ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_STAR/index/
reads=/ufrc/lee/braskey/Data/RNAseq/
aln=/ufrc/lee/braskey/Data/RNAseq/alignments/STAR_v2/

#gunzip -c ${reads}${sp}_${rep}_1.fastq.gz > ${reads}${sp}_${rep}_1.fastq
#gunzip -c ${reads}${sp}_${rep}_2.fastq.gz > ${reads}${sp}_${rep}_2.fastq

STAR --runThreadN 1 \
 --outFileNamePrefix ${aln}${sp}_${rep}_ \
 --genomeDir ${ref} \
 --readFilesIn ${reads}${sp}_${rep}_1.fastq ${reads}${sp}_${rep}_2.fastq \
 --outFilterScoreMinOverLread 0 \
 --outFilterMatchNminOverLread 0 \
 --outFilterMatchNmin 0 \
 --outFilterMismatchNmax 2 