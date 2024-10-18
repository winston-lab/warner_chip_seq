#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:05                    # Runtime in D-HH:MM format
#SBATCH --mem=750M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent

#Use older versions to avoid incompatibility with python 'dpp_inches'
module load gcc/6.2.0 python/2.7.12 deeptools/3.0.2	
	
for IP in V5 Flag; do

plotHeatmap -m deeptools/ratio/${IP}v8WG16_${1%}_si_ratio_scale.gz -o deeptools/plots/${IP}v8WG16_${1%}_si_ratio_scale_heatmap.png \
	--dpi 300 \
	--startLabel "TSS" \
	--endLabel "TES" \
	-y "normalized counts" \
        --plotTitle "${IP} over 8WG16 ${1%}" \
        --heatmapWidth 8 \
	--averageTypeSummaryPlot mean \

plotHeatmap -m deeptools/ratio/${IP}v8WG16_${1%}_si_ratio_reference.gz -o deeptools/plots/${IP}v8WG16_${1%}_si_ratio_reference_heatmap.png \
	--dpi 300 \
	--refPointLabel "TSS" \
	--sortRegions ascend \
	--whatToShow "heatmap and colorbar" \
	--sortUsing region_length \
	-y "normalized counts" \
        --plotTitle "${IP} over 8WG16 ${1%}" \
        --heatmapWidth 8 \

plotProfile -m deeptools/ratio/${IP}v8WG16_${1%}_si_ratio_scale.gz -o deeptools/plots/${IP}v8WG16_${1%}_si_ratio_scale_profile.png \
	--dpi 300 \
	--plotHeight 12 \
	--plotWidth 12 \
	-y "normalized counts" \
	--plotTitle "${IP} over 8WG16 ${1%}" \
	--startLabel "TSS" \
	--endLabel "TES" \
	--perGroup \
	--averageType mean \

done
