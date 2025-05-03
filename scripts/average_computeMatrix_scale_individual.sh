#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:10                    # Runtime in D-HH:MM format
#SBATCH --mem=200M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent


module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0

for cond in 93_D 93_I 95_D 95_I; do

for IP in V5 Flag 8WG16; do

computeMatrix scale-regions -S deeptools/averaged/${cond}_${IP}_averaged_si_ratio.bw -R genome/annotations/Scer_transcripts_w_verifiedORFs-nonoverlapping_averagedRpb1sorted.bed -m 1000 -o deeptools/averaged/${cond}_${IP}_averaged_si_ratio_scale.gz \
	--outFileNameMatrix deeptools/averaged/tab/${cond}_${IP}_averaged_si_ratio_scale.tab \
	-a 250 \
	-b 250 \
	-bs 10 \
	--averageTypeBins mean \
	--sortRegions keep \

done

done


