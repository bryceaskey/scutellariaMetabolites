#!/bin/bash
#SBATCH --job-name=variant-calling-ALT      # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=720                          # Time limit hrs:min:sec
#SBATCH --output=variant-calling-ALT_%j.log # Standard output and error log

pwd; hostname; date

module bcftools/1.10.2

echo "Calling variants for aligned altissima genomic data"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/GCA005771605.fa
sorted_aln=/ufrc/lee/braskey/Data/WGS/ALT_aln_sorted.bam
output_dir=/ufrc/lee/braskey/Data/WGS/

bcftools mpileup -Ov -f ${ref} ${sorted_aln} | bcftools call -c -V indels -o ${output_dir}ALT_aln.vcf
bcftools filter -Ov -i 'DP>5' ${output_dir}ALT_aln.vcf -o ${output_dir}ALT_aln_filtered.vcf
bcftools view ${output_dir}ALT_aln_filtered.vcf -Oz -o ${output_dir}ALT_aln_filtered.vcf.gz
bcftools index ${output_dir}ALT_aln_filtered.vcf.gz