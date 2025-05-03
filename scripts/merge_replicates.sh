#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 2                                 # Requested cores
#SBATCH --time=0-00:10                    # Runtime in D-HH:MM format
#SBATCH --mem=120M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent


for cond in 93_D 93_I 95_D 95_I; do

for IP in V5 8WG16 Flag input; do

cat fastq/${cond}*${IP}*rep*R1*fastq.gz > fastq/${cond}_${IP}_merged_R1_001.fastq.gz
cat fastq/${cond}*${IP}*rep*R2*fastq.gz > fastq/${cond}_${IP}_merged_R2_001.fastq.gz

done

done