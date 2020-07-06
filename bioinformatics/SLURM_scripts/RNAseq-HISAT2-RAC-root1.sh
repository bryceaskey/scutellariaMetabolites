#!/bin/bash
#SBATCH --job-name=RNAseq-HISAT2-RAC-root1       # Job name
#SBATCH --account=lee                            # Account name
#SBATCH --qos=lee                                # QOS name
#SBATCH --mail-type=END,FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu              # Where to send mail	
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --mem=8gb                                # Job memory request
#SBATCH --time=120:00:00                         # Time limit hrs:min:sec
#SBATCH --output=RNAseq-HISAT2-RAC-root1_%j.log  # Standard output and error log

pwd; hostname; date

module load adapterremoval/2.2.2 hisat2/2.2.0

echo "Mapping racemosa root1 RNAseq data to baicalensis reference genome"

sp=RAC
rep=root1
index=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/HISAT2-index/
reads=/ufrc/lee/braskey/Data/RNAseq/
aln=/ufrc/lee/braskey/Data/RNAseq/alignments/HISAT2_v1/

#gunzip -c ${reads}${sp}_${rep}_1.fastq.gz > ${reads}${sp}_${rep}_1.fastq
#gunzip -c ${reads}${sp}_${rep}_2.fastq.gz > ${reads}${sp}_${rep}_2.fastq

AdapterRemoval --file1 ${reads}${sp}_${rep}_1.fastq --file2 ${reads}${sp}_${rep}_2.fastq \
 --basename ${reads}${sp}_${rep} --output1 ${reads}${sp}_${rep}_trimmed_1.fastq --output2 ${reads}${sp}_${rep}_trimmed_2.fastq \
 --trimns --trimqualities --minlength 80

hisat2 -x ${index}GCA005771605 \
  -1 ${reads}${sp}_${rep}_trimmed_1.fastq -2 ${reads}${sp}_${rep}_trimmed_2.fastq \
  -S ${aln}${sp}_${rep}.sam \
  --summary-file ${aln}${sp}_${rep}_summary.txt