#!/bin/bash
#SBATCH --job-name=chr3_P1groupcall.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prashantbhandari@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=72:00:00
#SBATCH --output=chr3_P1groupcall.%A_%a.out
#SBATCH --array=1-11
#SBATCH --qos=lee-b
pwd; hostname; date;
RUN=${SLURM_ARRAY_TASK_ID}

module load bcftools/1.9

echo '1+2'|bc

INPUT_DIR=/ufrc/lee/share/86538916F2/sorted_bams2020

INPUT_FILE_ONE=$(ls -1 sorted_bams2020/P1*.bwa.sorted.bam | sed -n ${RUN}p)
SAMPLE=$(basename "${INPUT_FILE_ONE}" .bwa.sorted.bam)

OUTPUT_DIR=/ufrc/lee/share/86538916F2/chr_3

bcftools mpileup -Ov -f /ufrc/lee/share/86538916F2/SL4.0/S_lycopersicum_chromosomes.4.00.fa -r SL4.0ch03 ${INPUT_DIR}/${SAMPLE}.bwa.sorted.bam | bcftools call -c -V indels -o ${OUTPUT_DIR}/${SAMPLE}_call_1.vcf   ## individual sorted.bam file
bcftools filter -i 'DP>5' ${OUTPUT_DIR}/${SAMPLE}_call_1.vcf -o ${OUTPUT_DIR}/${SAMPLE}_call_2.vcf
bcftools view ${OUTPUT_DIR}/${SAMPLE}_call_2.vcf -Oz -o ${OUTPUT_DIR}/${SAMPLE}_call_2.vcf.gz
bcftools index ${OUTPUT_DIR}/${SAMPLE}_call_2.vcf.gz

echo '1+2'|bc



