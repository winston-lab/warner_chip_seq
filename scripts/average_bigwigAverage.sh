#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=0-00:10                    # Runtime in D-HH:MM format
#SBATCH --mem=400M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<your_email_here>          # Email to which notifications will be sent


module load gcc/9.2.0 python/3.10.11

source env/deeptools/bin/activate


for cond in 93_D 93_I 95_D 95_I; do

for IP in Flag V5 8WG16 input; do

bigwigAverage -b deeptools/si/${cond}_${IP}_rep*_si.bw  -o deeptools/averaged/${cond}_${IP}_averaged_midpoints_si.bw \
	-bs 10 \
	-p max

done

done

deactivate
