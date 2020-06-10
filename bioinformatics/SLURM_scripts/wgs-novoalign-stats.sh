#!/bin/bash
#SBATCH --job-name=wgs-novoalign-stats             # Job name
#SBATCH --account=lee                            # Account name
#SBATCH --qos=lee                                # QOS name
#SBATCH --mail-type=END,FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu              # Where to send mail	
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --mem=4gb                                # Job memory request
#SBATCH --time=744:00:00                         # Time limit hrs:min:sec
#SBATCH --output=wgs-novoalign-stats%j.log       # Standard output and error log

pwd; hostname; date

module load samtools/1.10

echo "Calculating stats for alignments generated with novoalign"

aln=/ufrc/lee/braskey/Data/WGS/novoaln/

samtools stats ${aln}ALT_novoaln_sorted.bam
samtools stats ${aln}BAR_novoaln_sorted.bam
samtools stats ${aln}HAV_novoaln_sorted.bam
samtools stats ${aln}RAC_novoaln_sorted.bam