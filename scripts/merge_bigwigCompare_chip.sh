#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 2                                 # Requested cores
#SBATCH --time=0-00:15                   # Runtime in D-HH:MM format
#SBATCH --mem=800M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent

module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0


for cond in 93_D 93_I 95_D 95_I; do

for IP in Flag V5 8WG16; do

bigwigCompare -b1 deeptools/si/${cond}_${IP}_merged_si.bw -b2 deeptools/si/${cond}_input_merged_si.bw -o deeptools/ratio/${cond}_${IP}_merged_si_ratio.bw \
	--operation ratio \
	-bs 20 \
	-p max

done

done
