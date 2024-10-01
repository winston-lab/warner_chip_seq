#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:05                    # Runtime in D-HH:MM format
#SBATCH --mem=500M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=james_warner@hms.harvard.edu          # Email to which notifications will be sent

#Use older versions to avoid incompatibility with python 'dpp_inches'
module load gcc/6.2.0 python/2.7.12 deeptools/3.0.2	
	
for IP in V5 8WG16 Flag; do

plotHeatmap -m deeptools/log2/si/${IP}vinput_rep4_si_log2.gz -o deeptools/plots/${IP}vinput_rep4_si_log2_heatmap.png \
	--dpi 300 \
	--startLabel "TSS" \
	--endLabel "TES" \
	-y "normalized counts" \
        --plotTitle "${IP} over input log2" \
        --heatmapWidth 12 \
	--averageTypeSummaryPlot mean \

plotProfile -m deeptools/log2/si/${IP}vinput_rep4_si_log2.gz -o deeptools/plots/${IP}vinput_rep4_si_log2_profile.png \
	--dpi 300 \
	--plotHeight 12 \
	--plotWidth 12 \
	-y "normalized counts" \
	--plotTitle "${IP} over input log2" \
	--startLabel "TSS" \
	--endLabel "TES" \
	--perGroup \
	--averageType mean \

done
