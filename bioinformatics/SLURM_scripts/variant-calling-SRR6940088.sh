#!/bin/bash
#SBATCH --job-name=variant-calling-SRR6940088         # Job name
#SBATCH --account=lee                                 # Account name
#SBATCH --qos=lee                                     # QOS name
#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu                   # Where to send mail	
#SBATCH --ntasks=1                                    # Run on a single CPU
#SBATCH --mem=4gb                                     # Job memory request
#SBATCH --time=120:00:00                              # Time limit hrs:min:sec
#SBATCH --output=variant-calling-SRR6940088_%j.log    # Standard output and error log

pwd; hostname; date

module load bcftools/1.10.2

echo "Calling variants for aligned SRR6940088 (baicalensis) WGS data"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/GCA005771605.fa
sorted_aln=/ufrc/lee/braskey/Data/SRP096180/SRR6940088_aln_sorted.bam
output_dir=/ufrc/lee/braskey/Data/SRP096180/

bcftools mpileup -Ov -Q 20 -f ${ref} ${sorted_aln} | bcftools call -c -o ${output_dir}SRR6940088_aln.vcf
bcftools filter -Ov -i '%QUAL>20' ${output_dir}SRR6940088_aln.vcf -o ${output_dir}SRR6940088_aln_flt.vcf
bcftools index ${output_dir}SRR6940088_aln_flt.vcf
bcftools stats ${output_dir}SRR6940088_aln_flt.vcf
bcftools view ${output_dir}SRR6940088_aln_flt.vcf -Oz -o ${output_dir}SRR6940088_aln_flt.vcf.gz
bcftools index ${output_dir}SRR6940088_aln_flt.vcf.gz