#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 2                                 # Requested cores
#SBATCH --time=0-00:10                   # Runtime in D-HH:MM format
#SBATCH --mem=800M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=james_warner@hms.harvard.edu          # Email to which notifications will be sent

module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0


bigwigCompare -b1 deeptools/si/93_D_V5_${1%}*_si.bw -b2 deeptools/si/93_D_input_${1%}*_si.bw -o deeptools/log2/si/93_D_V5vinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/93_I_V5_${1%}*_si.bw -b2 deeptools/si/93_I_input_${1%}*_si.bw -o deeptools/log2/si/93_I_V5vinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/95_D_V5_${1%}*_si.bw -b2 deeptools/si/95_D_input_${1%}*_si.bw -o deeptools/log2/si/95_D_V5vinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/95_I_V5_${1%}*_si.bw -b2 deeptools/si/95_I_input_${1%}*_si.bw -o deeptools/log2/si/95_I_V5vinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/93_D_8WG16_${1%}*_si.bw -b2 deeptools/si/93_D_input_${1%}*_si.bw -o deeptools/log2/si/93_D_8WG16vinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/93_I_8WG16_${1%}*_si.bw -b2 deeptools/si/93_I_input_${1%}*_si.bw -o deeptools/log2/si/93_I_8WG16vinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/95_D_8WG16_${1%}*_si.bw -b2 deeptools/si/95_D_input_${1%}*_si.bw -o deeptools/log2/si/95_D_8WG16vinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/95_I_8WG16_${1%}*_si.bw -b2 deeptools/si/95_I_input_${1%}*_si.bw -o deeptools/log2/si/95_I_8WG16vinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/93_D_Flag_${1%}*_si.bw -b2 deeptools/si/93_D_input_${1%}*_si.bw -o deeptools/log2/si/93_D_Flagvinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/93_I_Flag_${1%}*_si.bw -b2 deeptools/si/93_I_input_${1%}*_si.bw -o deeptools/log2/si/93_I_Flagvinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/95_D_Flag_${1%}*_si.bw -b2 deeptools/si/95_D_input_${1%}*_si.bw -o deeptools/log2/si/95_D_Flagvinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max

bigwigCompare -b1 deeptools/si/95_I_Flag_${1%}*_si.bw -b2 deeptools/si/95_I_input_${1%}*_si.bw -o deeptools/log2/si/95_I_Flagvinput_${1%}_si_log2.bw \
        --operation log2 \
        -bs 20 \
        -p max
