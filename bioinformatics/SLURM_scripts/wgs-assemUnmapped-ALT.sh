#!/bin/bash
#SBATCH --job-name=wgs-assemUnmapped-ALT       # Job name
#SBATCH --account=lee                          # Account name
#SBATCH --qos=lee                              # QOS name
#SBATCH --mail-type=END,FAIL                   # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu            # Where to send mail	
#SBATCH --ntasks=24                            
#SBATCH --cpus-per-task=1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=740:00:00                       # Time limit hrs:min:sec
#SBATCH --output=wgs-assemUnmapped-ALT_%j.log  # Standard output and error log

pwd; hostname; date

module load  intel/2018.1.163 boost/1.71.0 openmpi/3.1.0 sparsehash/2.0.3 abyss/2.1.0

echo "Assembling unmapped altissima reads into contigs"

sp=ALT
reads=/ufrc/lee/braskey/Data/WGS/bwa/unmapped/
contigs=/ufrc/lee/braskey/Data/WGS/bwa/contigs/

abyss-pe name=${contigs}${sp}_unmapped k=64 in="${reads}${sp}_unmapped_1.fq ${reads}${sp}_unmapped_2.fq"