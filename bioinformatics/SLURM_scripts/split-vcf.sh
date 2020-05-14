#!/bin/bash
#SBATCH --job-name=split-vcf                # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=120:00:00                    # Time limit hrs:min:sec
#SBATCH --output=split-vcf_%j.log           # Standard output and error log

pwd; hostname; date

module load bcftools/1.10.2.1

echo "Separating .vcf files by chromosome"

dir=/ufrc/lee/braskey/Data/WGS/bwa/

for name in ALT BAR HAV RAC
do
  for chr in {1..9}
  do
    bcftools view -Oz -r CM01672${chr}.1 ${dir}${name}_aln_flt.vcf.gz > ${dir}${name}_aln_flt_chr${chr}.vcf.gz
  done
done