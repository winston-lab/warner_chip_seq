#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:12                    # Runtime in D-HH:MM format
#SBATCH --mem=150M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent


module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0


for IP in V5 Flag 8WG16; do

computeMatrix scale-regions -S deeptools/averaged/*${IP}_averaged_si_ratio.bw -R genome/annotations/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed -m 1000 -o deeptools/averaged/${IP}_averaged_si_ratio_scale.gz \
	--outFileNameMatrix deeptools/averaged/tab/${IP}_averaged_si_ratio_scale.tab \
	-a 250 \
	-b 250 \
	-bs 10 \
	--averageTypeBins mean \
	--sortUsing mean \
	--sortRegions descend \
	--sortUsingSamples 2 \

done


computeMatrix scale-regions -S deeptools/averaged/93_I_8WG16_averaged_si_ratio.bw -R genome/annotations/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed -m 1000 -o deeptools/averaged/93_I_8WG16_averaged_si_ratio_scale_sorted.gz \
	--outFileSortedRegions genome/annotations/Scer_transcripts_w_verifiedORFs-nonoverlapping_averagedRpb1sorted.bed \
	-a 250 \
	-b 250 \
	-bs 10 \
	--sortRegions descend \
	--sortUsing mean \
	--averageTypeBins mean \


