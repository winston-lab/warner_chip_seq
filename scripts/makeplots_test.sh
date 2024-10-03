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

plotHeatmap -m deeptools/log2/si/${IP}vinput_${1%}_si_log2_scale_test.gz -o deeptools/plots/${IP}vinput_${1%}_si_log2_scale_heatmap_test.png \
	--dpi 300 \
	--startLabel "TSS" \
	--endLabel "TES" \
	-y "normalized counts" \
        --plotTitle "${IP} over input log2 ${1%} TEST" \
        --heatmapWidth 8 \
	--averageTypeSummaryPlot mean \

plotHeatmap -m deeptools/log2/si/${IP}vinput_${1%}_si_log2_reference_test.gz -o deeptools/plots/${IP}vinput_${1%}_si_log2_reference_heatmap_test.png \
	--dpi 300 \
	--refPointLabel "TSS" \
	--sortRegions ascend \
	--sortUsing region_length \
	--whatToShow heatmap and colorbar \
	-y "normalized counts" \
        --plotTitle "${IP} over input log2 ${1%} TEST" \
        --heatmapWidth 8 \

plotProfile -m deeptools/log2/si/${IP}vinput_${1%}_si_log2_scale_test.gz -o deeptools/plots/${IP}vinput_${1%}_si_log2_scale_profile_test.png \
	--dpi 300 \
	--plotHeight 12 \
	--plotWidth 12 \
	-y "normalized counts" \
	--plotTitle "${IP} over input log2 ${1%} TEST" \
	--startLabel "TSS" \
	--endLabel "TES" \
	--perGroup \
	--averageType mean \

done
