#!/bin/bash
##all shell conditions
#SBATCH --job-name=sample_1_step_1.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prashantbhandari@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6gb
#SBATCH --time=72:00:00
#SBATCH --output=sample_1_step_1.%A_%a.out
#SBATCH --array=1-4
#SBATCH --qos=lee-b
echo '1+2'|bc 
pwd; hostname; date;

RUN=${SLURM_ARRAY_TASK_ID}

module load bwa/0.7.17 samtools/1.9

INPUT_DIR=/ufrc/lee/share/map1/sample_1
BWA_REF_DIR=/ufrc/lee/share/86538916F2/SL4.0
GENOME=S_lycopersicum_chromosomes.4.00.fa

INPUT_FILE_ONE=$(ls -1 sample_1/*_1.fq | sed -n ${RUN}p)
SAMPLE=$(basename "${INPUT_FILE_ONE}" _1.fq)
echo "RUN #${RUN} with sample ${SAMPLE}"

OUT_DIR=/ufrc/lee/share/map1/imfiles
echo "bwa mem ${BWA_REF_DIR}/${GENOME} ${INPUT_DIR}/${SAMPLE}_1.fq ${INPUT_DIR}/${SAMPLE}_2.fq | samtools sort -o ${OUT_DIR}/${SAMPLE}.bwa.sorted.bam"
bwa mem ${BWA_REF_DIR}/${GENOME} ${INPUT_DIR}/${SAMPLE}_1.fq ${INPUT_DIR}/${SAMPLE}_2.fq | samtools sort -o ${OUT_DIR}/${SAMPLE}.bwa.sorted.bam 
samtools index ${OUT_DIR}/${SAMPLE}.bwa.sorted.bam

date
echo '1+2'|bc   

