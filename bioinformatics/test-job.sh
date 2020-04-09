#!/bin/sh
#SBATCH --job-name=serial_job_test    # Job name
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=braskey@ufl.edu   # Where to send mail	
#SBATCH --nodes=1                     # Use one node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --mem=600mb                   # Memory limit
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.out   # Standard output and error log

pwd; hostname; date

module load python

python /ufrc/lee/braskey/test-job.py