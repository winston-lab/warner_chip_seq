#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:15                    # Runtime in D-HH:MM format
#SBATCH --mem=2G                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent

module load gcc/9.2.0 java/jdk-1.8u112 qualimap/2.2.1 

for cond in 93_D 93_I 95_D 95_I; do

qualimap bamqc -bam bam/${cond}_input_merged_sorted.bam --outdir QC/ --outfile ${cond}_input_merged_qualimap.pdf --outformat PDF --paint-chromosome-limits

done




