#!/bin/bash
#SBATCH --job-name=variant-calling-RAC      # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=120:00:00                    # Time limit hrs:min:sec
#SBATCH --output=variant-calling-RAC_%j.log # Standard output and error log

pwd; hostname; date

module load bcftools/1.10.2

echo "Calling variants for aligned racemosa genomic data"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/GCA005771605.fa
sorted_aln=/ufrc/lee/braskey/Data/WGS/RAC_aln_sorted.bam
output_dir=/ufrc/lee/braskey/Data/WGS/

bcftools mpileup -Ov -f ${ref} ${sorted_aln} | bcftools call -c -V indels -o ${output_dir}RAC_aln.vcf
bcftools filter -Ov -i 'DP>7' ${output_dir}RAC_aln.vcf -o ${output_dir}RAC_aln_filtered.vcf
bcftools view ${output_dir}RAC_aln_filtered.vcf -Oz -o ${output_dir}RAC_aln_filtered.vcf.gz
bcftools index ${output_dir}RAC_aln_filtered.vcf.gz