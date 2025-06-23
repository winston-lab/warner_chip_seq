#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:15                    # Runtime in D-HH:MM format
#SBATCH --mem=300M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent


module load gcc/9.2.0 python/3.10.11

source env/deeptools/bin/activate


for cond in 93_D 93_I 95_D 95_I; do

for IP in V5 Flag 8WG16; do

computeMatrix reference-point -S deeptools/averaged/ratio/${cond}_${IP}vinput_averaged_midpoints_si_ratio.bw -R genome/annotations/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed -o deeptools/averaged/ratio/gz/${cond}_${IP}vinput_averaged_midpoints_si_ratio_reference.gz \
	--outFileNameMatrix deeptools/averaged/ratio/tab/${cond}_${IP}vinput_averaged_midpoints_si_ratio_reference.tab \
	-a 4500 \
	-b 250 \
	-bs 10 \
	--averageTypeBins mean \
	--nanAfterEnd \
	--sortRegions ascend \
	--sortUsing region_length \

done

done


for cond in 93_D 93_I 95_D 95_I; do

for IP in V5 Flag; do

computeMatrix reference-point -S deeptools/averaged/perpol/${cond}_${IP}v8WG16_averaged_midpoints_si_ratio.bw -R genome/annotations/Scer_transcripts_w_verifiedORFs-nonoverlapping.bed -o deeptools/averaged/perpol/gz/${cond}_${IP}v8WG16_averaged_midpoints_si_ratio_reference.gz \
	--outFileNameMatrix deeptools/averaged/perpol/tab/${cond}_${IP}v8WG16_averaged_midpoints_si_ratio_reference.tab \
	-a 4500 \
	-b 250 \
	-bs 10 \
	--averageTypeBins mean \
	--nanAfterEnd \
	--sortRegions ascend \
	--sortUsing region_length \

done

done

deactivate
