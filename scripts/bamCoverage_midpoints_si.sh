#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=0-00:12                    # Runtime in D-HH:MM format
#SBATCH --mem=250M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>         # Email to which notifications will be sent


module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0

base=$(basename ${1%} _sorted.bam)

alpha=$(grep ${base} logs/normalization_table.csv | cut -f2 -d,)

bamCoverage -b ${1%} -o deeptools/si/${base}_midpoints_si.bw \
        -bs 20 \
	--scaleFactor ${alpha} \
        --extendReads \
	--centerReads \
        -p max \
        --smoothLength 60 \
        --ignoreForNormalization chrM \
