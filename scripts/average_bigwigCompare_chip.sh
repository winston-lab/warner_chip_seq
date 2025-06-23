#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 2                                 # Requested cores
#SBATCH --time=0-00:10                   # Runtime in D-HH:MM format
#SBATCH --mem=800M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent

module load gcc/9.2.0 python/3.10.11

source env/deeptools/bin/activate


for cond in 93_D 93_I 95_D 95_I; do

for IP in Flag V5 8WG16; do

bigwigCompare -b1 deeptools/averaged/${cond}_${IP}_averaged_midpoints_si.bw -b2 deeptools/averaged/${cond}_input_averaged_midpoints_si.bw -o deeptools/averaged/ratio/${cond}_${IP}vinput_averaged_midpoints_si_ratio.bw \
	--operation ratio \
	-bs 10 \
	--pseudocount 0.1 \
	-p max

done

done


for cond in 93_D 93_I 95_D 95_I; do

for IP in Flag V5; do

bigwigCompare -b1 deeptools/averaged/ratio/${cond}_${IP}vinput_averaged_midpoints_si_ratio.bw -b2 deeptools/averaged/ratio/${cond}_8WG16vinput_averaged_midpoints_si_ratio.bw -o deeptools/averaged/perpol/${cond}_${IP}v8WG16_averaged_midpoints_si_ratio.bw \
	--operation ratio \
	-bs 10 \
	--pseudocount 0.1 \
	-p max

done

done

deactivate
