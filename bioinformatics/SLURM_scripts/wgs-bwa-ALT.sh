#!/bin/bash
#SBATCH --job-name=wgs-bwa-ALT              # Job name
#SBATCH --account=lee                       # Account name
#SBATCH --qos=lee                           # QOS name
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu         # Where to send mail	
#SBATCH --ntasks=1                          # Run on a single CPU
#SBATCH --mem=4gb                           # Job memory request
#SBATCH --time=120:00:00                    # Time limit hrs:min:sec
#SBATCH --output=wgs-bwa-ALT_%j.log         # Standard output and error log

pwd; hostname; date

module load fastx_toolkit/0.0.14 bwa/0.7.17 samtools/1.10

echo "Aligning altissima WGS data to reference genome"

sp=ALT
ref=/ufrc/lee/braskey/Data/ASM577160v1/ncbi_dataset/data/GCA_005771605.1/GCA005771605.fa
wgs=/ufrc/lee/braskey/Data/WGS/bwa/

gunzip -c ${wgs}${sp}_1.fq.gz > ${wgs}${sp}_1.fq
gunzip -c ${wgs}${sp}_2.fq.gz > ${wgs}${sp}_2.fq

fastq_quality_filter -Q 33 -q 36 -i ${wgs}${sp}_1.fq > ${wgs}${sp}_1_flt.fq
fastq_quality_filter -Q 33 -q 36 -i ${wgs}${sp}_2.fq > ${wgs}${sp}_2_flt.fq

bwa index ${ref}
bwa mem ${ref} ${wgs}${sp}_1_flt.fq ${wgs}${sp}_2_flt.fq > ${wgs}${sp}_aln.sam

samtools fixmate -O bam ${wgs}${sp}_aln.sam ${wgs}${sp}_aln_fixmate.bam
samtools sort -O bam -o ${wgs}${sp}_aln_sorted.bam ${wgs}${sp}_aln_fixmate.bam
samtools stats ${wgs}${sp}_aln_sorted.bam