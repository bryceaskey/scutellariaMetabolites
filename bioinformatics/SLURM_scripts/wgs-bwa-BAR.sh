#!/bin/bash
#SBATCH --job-name=wgs-alignment-BAR        # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=120:00:00                    # Time limit hrs:min:sec
#SBATCH --output=wgs-alignment-BAR_%j.log   # Standard output and error log

pwd; hostname; date

module load bwa/0.7.17 samtools/1.10

echo "Aligning barbata WGS data to reference genome"

sp=BAR
ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/GCA005771605.fa
wgs=/ufrc/lee/braskey/Data/WGS/bwa/

bwa index ${ref}
bwa mem ${ref} ${wgs}${sp}_1.fq.gz ${wgs}${sp}_2.fq.gz > ${wgs}${sp}_aln.sam
samtools fixmate -O bam ${wgs}${sp}_aln.sam ${wgs}${sp}_aln_fixmate.bam
samtools sort -O bam -o ${wgs}${sp}_aln_sorted.bam ${wgs}${sp}_aln_fixmate.bam
samtools stats ${wgs}${sp}_aln_sorted.bam