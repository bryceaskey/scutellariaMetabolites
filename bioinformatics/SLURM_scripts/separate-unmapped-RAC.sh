#!/bin/bash
#SBATCH --job-name=separate-unmapped-RAC       # Job name
#SBATCH --account=lee                          # Account name
#SBATCH --qos=lee                              # QOS name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu            # Where to send mail	
#SBATCH --ntasks=1                             # Run on a single CPU
#SBATCH --mem=16gb                             # Job memory request
#SBATCH --time=120:00:00                       # Time limit hrs:min:sec
#SBATCH --output=separate-unmapped-RAC_%j.log  # Standard output and error log

pwd; hostname; date

module load samtools/1.10 seqtk/1.3

echo "Separating mapped and unmapped reads for racemosa"

sp=RAC
wgs=/ufrc/lee/braskey/Data/WGS/bwa/
aln=/ufrc/lee/braskey/Data/WGS/bwa/alignments/

samtools view -F4 ${aln}${sp}_aln_sorted.bam > ${aln}${sp}_mapped.sam
samtools view -f4 ${aln}${sp}_aln_sorted.bam > ${aln}${sp}_unmapped.sam

cut -f1 ${aln}${sp}_mapped.sam | sort | uniq > ${aln}${sp}_mapped_ids.lst
cut -f1 ${aln}${sp}_unmapped.sam | sort | uniq > ${aln}${sp}_unmapped_ids.lst

seqtk subseq ${wgs}${sp}_1.fq ${aln}${sp}_mapped_ids.lst > ${aln}${sp}_1_mapped.fq
seqtk subseq ${wgs}${sp}_2.fq ${aln}${sp}_mapped_ids.lst > ${aln}${sp}_2_mapped.fq
seqtk subseq ${wgs}${sp}_1.fq ${aln}${sp}_unmapped_ids.lst > ${aln}${sp}_1_unmapped.fq
seqtk subseq ${wgs}${sp}_2.fq ${aln}${sp}_unmapped_ids.lst > ${aln}${sp}_2_unmapped.fq