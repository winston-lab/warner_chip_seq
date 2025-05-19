#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 2                                 # Requested cores
#SBATCH --time=0-00:10                   # Runtime in D-HH:MM format
#SBATCH --mem=800M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent

module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0

for strain in 661 666; do

for plasmid in 93 95; do

for IP in 8WG16 myc Flag; do

bigwigCompare -b1 deeptools/si/${strain}_${plasmid}_${IP}_${1%}_midpoints_si.bw -b2 deeptools/si/${strain}_${plasmid}_input_${1%}_midpoints_si.bw -o deeptools/ratio/${strain}_${plasmid}_${IP}vinput_${1%}_midpoints_si_ratio.bw \
	--operation ratio \
	-bs 20 \
	-p max

done

done

done


for strain in 661 666; do

for plasmid in 93 95; do

for IP in myc Flag; do

bigwigCompare -b1 deeptools/ratio/${strain}_${plasmid}_${IP}vinput_${1%}_midpoints_si_ratio.bw -b2 deeptools/ratio/${strain}_${plasmid}_8WG16vinput_${1%}_midpoints_si_ratio.bw -o deeptools/perpol/${strain}_${plasmid}_${IP}v8WG16_${1%}_midpoints_si_ratio.bw \
	--operation ratio \
	-bs 20 \
	-p max

done

done

done
