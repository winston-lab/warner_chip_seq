#!/bin/bash
# This file must be run in an interactive session with deeptools 3.5.6 installed.

srun --pty -t 0-2:0:0 -p interactive /bin/bash

module load gcc/6.2.0
module load python/2.
module load deeptools/3.5.6


for cond in 93_D 93_I 95_D 95_I; do

for IP in Flag V5 8WG16 input; do

bigwigAverage -b deeptools/si/${cond}_${IP}_rep*_si.bw  -o deeptools/averaged/${cond}_${IP}_averaged_si.bw \
	-bs 10 \
	-p max

done

done
