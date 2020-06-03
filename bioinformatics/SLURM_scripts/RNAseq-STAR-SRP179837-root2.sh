#!/bin/bash
#SBATCH --job-name=RNAseq-STAR-SRP179837-root2       # Job name
#SBATCH --account=lee                                # Account name
#SBATCH --qos=lee                                    # QOS name
#SBATCH --mail-type=END,FAIL                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                  # Where to send mail	
#SBATCH --ntasks=1                                   # Run on a single CPU
#SBATCH --mem=8gb                                    # Job memory request
#SBATCH --time=120:00:00                             # Time limit hrs:min:sec
#SBATCH --output=RNAseq-STAR-SRP179837-root2_%j.log  # Standard output and error log

pwd; hostname; date

module load star/2.7.3a

echo "Mapping SRP179837 (baicalensis) root2 RNAseq data to baicalensis reference genome"

sp=SRP179837
rep=root2
ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_STAR/index/
reads=/ufrc/lee/braskey/Data/SRP179837/
aln=/ufrc/lee/braskey/Data/SRP179837/alignments/

STAR --runThreadN 1 --outFilterMismatchNmax 2 --outFileNamePrefix ${aln}${sp}_${rep}_ --genomeDir ${ref} --readFilesIn ${reads}${sp}_${rep}_1.fastq ${reads}${sp}_${rep}_2.fastq