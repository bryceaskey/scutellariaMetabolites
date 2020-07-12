#!/bin/bash
#SBATCH --job-name=RNAseq-HISAT2-SRP179837-leaf1       # Job name
#SBATCH --account=lee                                  # Account name
#SBATCH --qos=lee                                      # QOS name
#SBATCH --mail-type=END,FAIL                           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                    # Where to send mail	
#SBATCH --ntasks=1                                     # Run on a single CPU
#SBATCH --mem=8gb                                      # Job memory request
#SBATCH --time=120:00:00                               # Time limit hrs:min:sec
#SBATCH --output=RNAseq-HISAT2-SRP179837-leaf1_%j.log  # Standard output and error log

pwd; hostname; date

module load trimmomatic/0.39 hisat2/2.2.0

echo "Mapping SRP179837 (baicalensis) leaf1 RNAseq data to baicalensis reference genome"

sp=SRP179837
rep=leaf1
index=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/HISAT2-index/
reads=/ufrc/lee/braskey/Data/SRP179837/
aln=/ufrc/lee/braskey/Data/SRP179837/alignments/HISAT2_v2/

trimmomatic PE -threads 1 \
  ${reads}${sp}_${rep}_1.fastq ${reads}${sp}_${rep}_2.fastq \
  ${reads}${sp}_${rep}_trimmed_1.fastq ${reads}${sp}_${rep}_unpaired_1.fastq \
  ${reads}${sp}_${rep}_trimmed_2.fastq ${reads}${sp}_${rep}_unpaired_2.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:100

hisat2 -x ${index}GCA005771605 \
  -1 ${reads}${sp}_${rep}_trimmed_1.fastq -2 ${reads}${sp}_${rep}_trimmed_2.fastq \
  -S ${aln}${sp}_${rep}.sam \
  --summary-file ${aln}${sp}_${rep}_summary.txt