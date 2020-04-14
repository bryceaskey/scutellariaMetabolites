#!/bin/bash
#SBATCH --job-name=wgs-alignment-ALT        # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=720                          # Time limit hrs:min:sec
#SBATCH --output=wgs-alignment-ALT_%j.log   # Standard output and error log

pwd; hostname; date

module load bwa/0.7.17 samtools/1.10

echo "Aligning altissima WGS data to reference genome"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/GCA005771605.fa
wgs_pe1=/ufrc/lee/braskey/Data/WGS/ALT_1.fq.gz
wgs_pe2=/ufrc/lee/braskey/Data/WGS/ALT_2.fq.gz
output_dir=/ufrc/lee/braskey/Data/WGS/

bwa index ${ref}
bwa mem ${ref} ${wgs_pe1} ${wgs_pe2} > ${output_dir}ALT_aln.sam
samtools fixmate -O bam ${output_dir}ALT_aln.sam ${output_dir}ALT_aln_fixmate.bam
samtools sort -O bam -o ${output_dir}ALT_aln_sorted.bam ${output_dir}ALT_aln_fixmate.bam