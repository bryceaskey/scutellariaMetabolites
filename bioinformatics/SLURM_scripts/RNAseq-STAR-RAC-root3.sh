#!/bin/bash
#SBATCH --job-name=RNAseq-STAR-RAC-root3       # Job name
#SBATCH --account=lee                          # Account name
#SBATCH --qos=lee                              # QOS name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu            # Where to send mail	
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --mem=8gb                              # Job memory request
#SBATCH --time=120:00:00                       # Time limit hrs:min:sec
#SBATCH --output=RNAseq-STAR-RAC-root3_%j.log  # Standard output and error log

pwd; hostname; date

module load star/2.7.3a trimmomatic/0.39

echo "Mapping racemosa root3 RNAseq data to baicalensis reference genome"

sp=RAC
rep=root3
ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_STAR/index/
reads=/ufrc/lee/braskey/Data/RNAseq/
aln=/ufrc/lee/braskey/Data/RNAseq/alignments/STAR_v3/

#gunzip -c ${reads}${sp}_${rep}_1.fastq.gz > ${reads}${sp}_${rep}_1.fastq
#gunzip -c ${reads}${sp}_${rep}_2.fastq.gz > ${reads}${sp}_${rep}_2.fastq

trimmomatic PE -threads 1 \
  ${reads}${sp}_${rep}_1.fastq ${reads}${sp}_${rep}_2.fastq \
  ${reads}${sp}_${rep}_trimmed_1.fastq ${reads}${sp}_${rep}_unpaired_1.fastq \
  ${reads}${sp}_${rep}_trimmed_2.fastq ${reads}${sp}_${rep}_unpaired_2.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

STAR --runThreadN 1 \
  --outFileNamePrefix ${aln}${sp}_${rep}_ \
  --genomeDir ${ref} \
  --readFilesIn ${reads}${sp}_${rep}_trimmed_1.fastq ${reads}${sp}_${rep}_trimmed_2.fastq \