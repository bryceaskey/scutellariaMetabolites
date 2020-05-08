#!/bin/bash
#SBATCH --job-name=bwa-test-HAV-variants     # Job name
#SBATCH --account=lee                         # Account name
#SBATCH --qos=lee                             # QOS name
#SBATCH --mail-type=END,FAIL                  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu           # Where to send mail	
#SBATCH --ntasks=1                            # Run on a single CPU
#SBATCH --mem=4gb                             # Job memory request
#SBATCH --time=120:00:00                      # Time limit hrs:min:sec
#SBATCH --output=bwa-test-HAV-variants_%j.log # Standard output and error log

pwd; hostname; date

module load bcftools/1.10.2

echo "Calling variants for aligned havanesis genomic data"

ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/copy_HAV_bwa/
wgs=/ufrc/lee/braskey/Data/WGS/bwa_test/

bcftools mpileup -Ov -f ${ref}GCA005771605.fa ${wgs}HAV_aln_sorted.bam | bcftools call -c -V indels -o ${wgs}HAV_aln.vcf
bcftools filter -Ov -i 'DP>7' ${wgs}HAV_aln.vcf -o ${wgs}HAV_aln_filtered.vcf
bcftools view ${wgs}HAV_aln_filtered.vcf -Oz -o ${wgs}HAV_aln_filtered.vcf.gz
bcftools index ${wgs}HAV_aln_filtered.vcf.gz