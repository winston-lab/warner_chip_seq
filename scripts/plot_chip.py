#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 23:33:06 2024

@author: jlwarner
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy


# ### Flag enrichment
# rep2_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep2_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep3_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep3_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep4_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep4_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)

# rep2_FvI = rep2_FvI.transpose()
# rep3_FvI = rep3_FvI.transpose()
# rep4_FvI = rep4_FvI.transpose()

# rep2_FvI_avg = np.array(rep2_FvI.mean(1))
# rep3_FvI_avg = np.array(rep3_FvI.mean(1))
# rep4_FvI_avg = np.array(rep4_FvI.mean(1))

# rep2_FvI_q1 = np.array(rep2_FvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
# rep2_FvI_q3 = np.array(rep2_FvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
# rep3_FvI_q1 = np.array(rep3_FvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
# rep3_FvI_q3 = np.array(rep3_FvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
# rep4_FvI_q1 = np.array(rep4_FvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
# rep4_FvI_q3 = np.array(rep4_FvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])

# rep2_WT_D = np.array(rep2_FvI_avg[0:150])
# rep2_WT_I = np.array(rep2_FvI_avg[150:300])
# rep2_mut_D = np.array(rep2_FvI_avg[300:450])
# rep2_mut_I = np.array(rep2_FvI_avg[450:600])

# rep3_WT_D = np.array(rep3_FvI_avg[0:150])
# rep3_WT_I = np.array(rep3_FvI_avg[150:300])
# rep3_mut_D = np.array(rep3_FvI_avg[300:450])
# rep3_mut_I = np.array(rep3_FvI_avg[450:600])

# rep4_WT_D = np.array(rep4_FvI_avg[0:150])
# rep4_WT_I = np.array(rep4_FvI_avg[150:300])
# rep4_mut_D = np.array(rep4_FvI_avg[300:450])
# rep4_mut_I = np.array(rep4_FvI_avg[450:600])


# positions = np.arange(-200,1300,10)
# ticks = np.array([-200, 0, 1000, 1300])
# labels = ['-0.2 kb', 'TSS', 'TES', '+0.2 kb']
# ybottom = 0
# ytop = 2.25
# yticks = np.array((0, 0.5, 1.0, 1.5, 2.0))


# plt.figure(dpi=300)
# #plt.plot(positions, rep2_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep2_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep2_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep2_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('FLAG ChIP rep 2')
# plt.xlabel('position')
# plt.ylabel('FLAG enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('FLAG_rep2.png')
# plt.show();

# plt.figure(dpi=300)
# #plt.plot(positions, rep3_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep3_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep3_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep3_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('FLAG ChIP rep 3')
# plt.xlabel('position')
# plt.ylabel('FLAG enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('Flag_rep3.png')
# plt.show();

# plt.figure(dpi=300)
# #plt.plot(positions, rep4_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep4_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep4_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep4_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('FLAG ChIP rep 4')
# plt.xlabel('position')
# plt.ylabel('FLAG enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('Flag_rep4.png')
# plt.show();


# ### Rpb1 enrichment
# rep2_PvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep2_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep3_PvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep3_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep4_PvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep4_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)

# rep2_PvI = rep2_PvI.transpose()
# rep3_PvI = rep3_PvI.transpose()
# rep4_PvI = rep4_PvI.transpose()

# rep2_PvI_avg = np.array(rep2_PvI.mean(1))
# rep3_PvI_avg = np.array(rep3_PvI.mean(1))
# rep4_PvI_avg = np.array(rep4_PvI.mean(1))

# rep2_WT_D = np.array(rep2_PvI_avg[0:150])
# rep2_WT_I = np.array(rep2_PvI_avg[150:300])
# rep2_mut_D = np.array(rep2_PvI_avg[300:450])
# rep2_mut_I = np.array(rep2_PvI_avg[450:600])

# rep3_WT_D = np.array(rep3_PvI_avg[0:150])
# rep3_WT_I = np.array(rep3_PvI_avg[150:300])
# rep3_mut_D = np.array(rep3_PvI_avg[300:450])
# rep3_mut_I = np.array(rep3_PvI_avg[450:600])

# rep4_WT_D = np.array(rep4_PvI_avg[0:150])
# rep4_WT_I = np.array(rep4_PvI_avg[150:300])
# rep4_mut_D = np.array(rep4_PvI_avg[300:450])
# rep4_mut_I = np.array(rep4_PvI_avg[450:600])


# ### MEDIAN
# # rep2_PvI_med = np.array(rep2_PvI.median(1))
# # rep3_PvI_med = np.array(rep3_PvI.median(1))
# # rep4_PvI_med = np.array(rep4_PvI.median(1))

# # rep2_WT_D = np.array(rep2_PvI_med[0:150])
# # rep2_WT_I = np.array(rep2_PvI_med[150:300])
# # rep2_mut_D = np.array(rep2_PvI_med[300:450])
# # rep2_mut_I = np.array(rep2_PvI_med[450:600])

# # rep3_WT_D = np.array(rep3_PvI_med[0:150])
# # rep3_WT_I = np.array(rep3_PvI_med[150:300])
# # rep3_mut_D = np.array(rep3_PvI_med[300:450])
# # rep3_mut_I = np.array(rep3_PvI_med[450:600])

# # rep4_WT_D = np.array(rep4_PvI_med[0:150])
# # rep4_WT_I = np.array(rep4_PvI_med[150:300])
# # rep4_mut_D = np.array(rep4_PvI_med[300:450])
# # rep4_mut_I = np.array(rep4_PvI_med[450:600])
# ###


# positions = np.arange(-200,1300,10)
# ticks = np.array([-200, 0, 1000, 1300])
# labels = ['-0.2 kb', 'TSS', 'TES', '+0.2 kb']
# ybottom = 0
# ytop = 2.75
# yticks = np.array((0, 0.5, 1.0, 1.5, 2.0, 2.5))


# plt.figure(dpi=300)
# #plt.plot(positions, rep2_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep2_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep2_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep2_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('Rpb1 ChIP rep 2')
# plt.xlabel('position')
# plt.ylabel('Rpb1 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('Rpb1_rep2.png')
# plt.show();

# plt.figure(dpi=300)
# #plt.plot(positions, rep3_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep3_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep3_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep3_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('Rpb1 ChIP rep 3')
# plt.xlabel('position')
# plt.ylabel('Rpb1 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('Rpb1_rep3.png')
# plt.show();

# plt.figure(dpi=300)
# #plt.plot(positions, rep4_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep4_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep4_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep4_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('Rpb1 ChIP rep 4')
# plt.xlabel('position')
# plt.ylabel('Rpb1 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('Rpb1_rep4.png')
# plt.show();


# ### Flag per Rpb1
# rep2_FvP = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep2_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep3_FvP = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep3_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep4_FvP = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep4_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)

# rep2_FvP = rep2_FvP.transpose()
# rep3_FvP = rep3_FvP.transpose()
# rep4_FvP = rep4_FvP.transpose()

# rep2_FvP_avg = np.array(rep2_FvP.mean(1))
# rep3_FvP_avg = np.array(rep3_FvP.mean(1))
# rep4_FvP_avg = np.array(rep4_FvP.mean(1))



# rep2_WT_D = np.array(rep2_FvP_avg[0:150])
# rep2_WT_I = np.array(rep2_FvP_avg[150:300])
# rep2_mut_D = np.array(rep2_FvP_avg[300:450])
# rep2_mut_I = np.array(rep2_FvP_avg[450:600])

# rep3_WT_D = np.array(rep3_FvP_avg[0:150])
# rep3_WT_I = np.array(rep3_FvP_avg[150:300])
# rep3_mut_D = np.array(rep3_FvP_avg[300:450])
# rep3_mut_I = np.array(rep3_FvP_avg[450:600])

# rep4_WT_D = np.array(rep4_FvP_avg[0:150])
# rep4_WT_I = np.array(rep4_FvP_avg[150:300])
# rep4_mut_D = np.array(rep4_FvP_avg[300:450])
# rep4_mut_I = np.array(rep4_FvP_avg[450:600])


# ### MEDIAN
# # rep2_FvP_med = np.array(rep2_FvP.median(1))
# # rep3_FvP_med = np.array(rep3_FvP.median(1))
# # rep4_FvP_med = np.array(rep4_FvP.median(1))

# # rep2_WT_D = np.array(rep2_FvP_med[0:150])
# # rep2_WT_I = np.array(rep2_FvP_med[150:300])
# # rep2_mut_D = np.array(rep2_FvP_med[300:450])
# # rep2_mut_I = np.array(rep2_FvP_med[450:600])

# # rep3_WT_D = np.array(rep3_FvP_med[0:150])
# # rep3_WT_I = np.array(rep3_FvP_med[150:300])
# # rep3_mut_D = np.array(rep3_FvP_med[300:450])
# # rep3_mut_I = np.array(rep3_FvP_med[450:600])

# # rep4_WT_D = np.array(rep4_FvP_med[0:150])
# # rep4_WT_I = np.array(rep4_FvP_med[150:300])
# # rep4_mut_D = np.array(rep4_FvP_med[300:450])
# # rep4_mut_I = np.array(rep4_FvP_med[450:600])
# ###


# positions = np.arange(-200,1300,10)
# ticks = np.array([-200, 0, 1000, 1300])
# labels = ['-0.2 kb', 'TSS', 'TES', '+0.2 kb']
# ybottom = 0
# ytop = 1.5
# yticks = np.array((0, 0.5, 1.0, 1.5))


# plt.figure(dpi=300)
# #plt.plot(positions, rep2_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep2_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep2_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep2_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('ChIP rep 2')
# plt.xlabel('position')
# plt.ylabel('Flag/Rpb1 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('FlagvRpb1_rep2.png')
# plt.show();

# plt.figure(dpi=300)
# #plt.plot(positions, rep3_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep3_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep3_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep3_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('ChIP rep 3')
# plt.xlabel('position')
# plt.ylabel('Flag/Rpb1 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('FlagvRpb1_rep3.png')
# plt.show();

# plt.figure(dpi=300)
# #plt.plot(positions, rep4_WT_D, label='WT_DMSO', color = 'black')
# plt.plot(positions, rep4_WT_I, label='Spt6', color = 'blue')
# #plt.plot(positions, rep4_mut_D, label='mut_DMSO', color = 'limegreen')
# plt.plot(positions, rep4_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('ChIP rep 4')
# plt.xlabel('position')
# plt.ylabel('Flag/Rpb1 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# plt.savefig('FlagvRpb1_rep4.png')
# plt.show();


# ### V5 enrichment
# rep2_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep2_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep3_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep3_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep4_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep4_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)

# rep2_VvI = rep2_VvI.transpose()
# rep3_VvI = rep3_VvI.transpose()
# rep4_VvI = rep4_VvI.transpose()

# rep2_VvI_avg = np.array(rep2_VvI.mean(1))
# rep3_VvI_avg = np.array(rep3_VvI.mean(1))
# rep4_VvI_avg = np.array(rep4_VvI.mean(1))

# rep2_VvI_q1 = np.array(rep2_VvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
# rep2_VvI_q3 = np.array(rep2_VvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
# rep3_VvI_q1 = np.array(rep3_VvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
# rep3_VvI_q3 = np.array(rep3_VvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
# rep4_VvI_q1 = np.array(rep4_VvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
# rep4_VvI_q3 = np.array(rep4_VvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])

# rep2_WT_D = np.array(rep2_VvI_avg[0:150])
# rep2_WT_I = np.array(rep2_VvI_avg[150:300])
# rep2_mut_D = np.array(rep2_VvI_avg[300:450])
# rep2_mut_I = np.array(rep2_VvI_avg[450:600])

# rep3_WT_D = np.array(rep3_VvI_avg[0:150])
# rep3_WT_I = np.array(rep3_VvI_avg[150:300])
# rep3_mut_D = np.array(rep3_VvI_avg[300:450])
# rep3_mut_I = np.array(rep3_VvI_avg[450:600])

# rep4_WT_D = np.array(rep4_VvI_avg[0:150])
# rep4_WT_I = np.array(rep4_VvI_avg[150:300])
# rep4_mut_D = np.array(rep4_VvI_avg[300:450])
# rep4_mut_I = np.array(rep4_VvI_avg[450:600])


# positions = np.arange(-200,1300,10)
# ticks = np.array([-200, 0, 1000, 1300])
# labels = ['-0.2 kb', 'TSS', 'TES', '+0.2 kb']
# ybottom = 0
# ytop = 1.75
# yticks = np.array((0, 0.5, 1.0, 1.5))


# plt.figure(dpi=300)
# plt.plot(positions, rep2_WT_D, label='Spt6_DMSO', color = 'black')
# plt.plot(positions, rep2_WT_I, label='Spt6', color = 'blue')
# plt.plot(positions, rep2_mut_D, label='Spt6Δ2-238_DMSO', color = 'limegreen')
# plt.plot(positions, rep2_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('V5 ChIP rep 2')
# plt.xlabel('position')
# plt.ylabel('V5 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# #plt.savefig('V5_rep2.png')
# plt.show();

# plt.figure(dpi=300)
# plt.plot(positions, rep3_WT_D, label='Spt6_DMSO', color = 'black')
# plt.plot(positions, rep3_WT_I, label='Spt6', color = 'blue')
# plt.plot(positions, rep3_mut_D, label='Spt6Δ2-238_DMSO', color = 'limegreen')
# plt.plot(positions, rep3_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('V5 ChIP rep 3')
# plt.xlabel('position')
# plt.ylabel('V5 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# #plt.savefig('V5_rep3.png')
# plt.show();

# plt.figure(dpi=300)
# plt.plot(positions, rep4_WT_D, label='Spt6_DMSO', color = 'black')
# plt.plot(positions, rep4_WT_I, label='Spt6', color = 'blue')
# plt.plot(positions, rep4_mut_D, label='Spt6Δ2-238_DMSO', color = 'limegreen')
# plt.plot(positions, rep4_mut_I, label='Spt6Δ2-238', color = 'orange')
# plt.title('V5 ChIP rep 4')
# plt.xlabel('position')
# plt.ylabel('V5 enrichment')
# plt.xlim(left=-200, right=1300)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylim(ybottom, ytop)
# plt.yticks(ticks=yticks)
# plt.legend()
# #plt.savefig('V5_rep4.png')
# plt.show();


### New plotting with interquartile range

# Flag
rep2_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)

rep2_FvI = rep2_FvI.transpose()
rep3_FvI = rep3_FvI.transpose()
rep4_FvI = rep4_FvI.transpose()

rep2_FvI_avg = np.array(rep2_FvI.mean(1))
rep3_FvI_avg = np.array(rep3_FvI.mean(1))
rep4_FvI_avg = np.array(rep4_FvI.mean(1))

rep2_FvI_q1 = np.array(rep2_FvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep2_FvI_q3 = np.array(rep2_FvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
rep3_FvI_q1 = np.array(rep3_FvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep3_FvI_q3 = np.array(rep3_FvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
rep4_FvI_q1 = np.array(rep4_FvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep4_FvI_q3 = np.array(rep4_FvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])


positions = np.arange(-250,1250,10)
ticks = np.array([-250, 0, 1000, 1250])
labels = ['-0.25 kb', 'TSS', 'TES', '+0.25 kb']
ybottom = 0
ytop = 2.5
yticks = np.array((0, 0.5, 1.0, 1.5, 2.0, 2.5))
replicates = [rep2_FvI, rep3_FvI, rep4_FvI]
names = ['rep2', 'rep3', 'rep4']
counter = 0

for df in replicates:
    avg = np.array(df.mean(1))
    q1 = np.array(df.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
    q3 = np.array(df.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
    plt.figure(dpi=300)
    plt.plot(positions, avg[150:300], label='Spt6', color='black')
    plt.fill_between(positions, q1[150:300], q3[150:300], alpha=0.2, color='black')
    plt.plot(positions, avg[450:600], label='Spt6Δ2-238', color='red')
    plt.fill_between(positions, q1[450:600], q3[450:600], alpha=0.2, color='red')
    plt.xlabel('position', fontsize=18)
    plt.ylabel('spike-in normalized\noccupancy', fontsize=18)
    plt.xlim(left=-250, right=1250)
    plt.xticks(ticks=ticks, labels=labels, fontsize=16)
    plt.ylim(ybottom, ytop)
    plt.yticks(ticks=yticks, fontsize=16)
    plt.legend(loc='upper left', fontsize=16, bbox_to_anchor=(1,1))
    plt.title('$\\alpha\mathdefault{-Flag}$', fontsize=20, pad=14)
    plt.savefig('Flag_'+names[counter]+'.svg')
    plt.show();
    counter+=1
    
    
# Rpb1
rep2_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)

rep2_RvI = rep2_RvI.transpose()
rep3_RvI = rep3_RvI.transpose()
rep4_RvI = rep4_RvI.transpose()

rep2_RvI_avg = np.array(rep2_RvI.mean(1))
rep3_RvI_avg = np.array(rep3_RvI.mean(1))
rep4_RvI_avg = np.array(rep4_RvI.mean(1))

rep2_RvI_q1 = np.array(rep2_RvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep2_RvI_q3 = np.array(rep2_RvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
rep3_RvI_q1 = np.array(rep3_RvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep3_RvI_q3 = np.array(rep3_RvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
rep4_RvI_q1 = np.array(rep4_RvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep4_RvI_q3 = np.array(rep4_RvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])


positions = np.arange(-250,1250,10)
ticks = np.array([-250, 0, 1000, 1250])
labels = ['-0.25 kb', 'TSS', 'TES', '+0.25 kb']
ybottom = 0
ytop = 3.0
yticks = np.array((0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))
replicates = [rep2_RvI, rep3_RvI, rep4_RvI]
names = ['rep2', 'rep3', 'rep4']
counter = 0

for df in replicates:
    avg = np.array(df.mean(1))
    q1 = np.array(df.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
    q3 = np.array(df.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
    plt.figure(dpi=300)
    plt.plot(positions, avg[150:300], label='Spt6', color='black')
    plt.fill_between(positions, q1[150:300], q3[150:300], alpha=0.2, color='black')
    plt.plot(positions, avg[450:600], label='Spt6Δ2-238', color='red')
    plt.fill_between(positions, q1[450:600], q3[450:600], alpha=0.2, color='red')
    plt.xlabel('position', fontsize=18)
    plt.ylabel('spike-in normalized\noccupancy', fontsize=18)
    plt.xlim(left=-250, right=1250)
    plt.xticks(ticks=ticks, labels=labels, fontsize=16)
    plt.ylim(ybottom, ytop)
    plt.yticks(ticks=yticks, fontsize=16)
    #plt.legend(loc='upper left', fontsize=16)
    plt.title('$\\alpha\mathdefault{-Rpb1}$', fontsize=20, pad=14)
    plt.savefig('8WG16_'+names[counter]+'.svg')
    plt.show();
    counter+=1

# Flag per Rpb1
rep2_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)

rep2_FvR = rep2_FvR.transpose()
rep3_FvR = rep3_FvR.transpose()
rep4_FvR = rep4_FvR.transpose()

rep2_FvR_avg = np.array(rep2_FvR.mean(1))
rep3_FvR_avg = np.array(rep3_FvR.mean(1))
rep4_FvR_avg = np.array(rep4_FvR.mean(1))

rep2_FvR_q1 = np.array(rep2_FvR.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep2_FvR_q3 = np.array(rep2_FvR.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
rep3_FvR_q1 = np.array(rep3_FvR.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep3_FvR_q3 = np.array(rep3_FvR.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
rep4_FvR_q1 = np.array(rep4_FvR.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep4_FvR_q3 = np.array(rep4_FvR.quantile(q=[0.25,0.75], axis=1).iloc[1,:])


positions = np.arange(-250,1250,10)
ticks = np.array([-250, 0, 1000, 1250])
labels = ['-0.25 kb', 'TSS', 'TES', '+0.25 kb']
ybottom = 0
ytop = 1.5
yticks = np.array((0, 0.5, 1.0, 1.5))
replicates = [rep2_FvR, rep3_FvR, rep4_FvR]
names = ['rep2', 'rep3', 'rep4']
counter = 0

for df in replicates:
    avg = np.array(df.mean(1))
    q1 = np.array(df.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
    q3 = np.array(df.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
    plt.figure(dpi=300)
    plt.plot(positions, avg[150:300], label='Spt6', color='black')
    plt.fill_between(positions, q1[150:300], q3[150:300], alpha=0.2, color='black')
    plt.plot(positions, avg[450:600], label='Spt6Δ2-238', color='red')
    plt.fill_between(positions, q1[450:600], q3[450:600], alpha=0.2, color='red')
    plt.xlabel('position', fontsize=18)
    plt.ylabel('spike-in normalized\noccupancy', fontsize=18)
    plt.xlim(left=-250, right=1250)
    plt.xticks(ticks=ticks, labels=labels, fontsize=16)
    plt.ylim(ybottom, ytop)
    plt.yticks(ticks=yticks, fontsize=16)
    #plt.legend(loc='upper left', fontsize=16)
    plt.title('$\\frac{\\alpha \mathdefault{-Flag}}{\\alpha \mathdefault{-Rpb1}}$', fontsize=26, pad=20)
    plt.savefig('Flag_v_8WG16_'+names[counter]+'.svg')
    plt.show();
    counter+=1
    

# V5
rep2_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)

rep2_VvI = rep2_VvI.transpose()
rep3_VvI = rep3_VvI.transpose()
rep4_VvI = rep4_VvI.transpose()

rep2_VvI_avg = np.array(rep2_VvI.mean(1))
rep3_VvI_avg = np.array(rep3_VvI.mean(1))
rep4_VvI_avg = np.array(rep4_VvI.mean(1))

rep2_VvI_q1 = np.array(rep2_VvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep2_VvI_q3 = np.array(rep2_VvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
rep3_VvI_q1 = np.array(rep3_VvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep3_VvI_q3 = np.array(rep3_VvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
rep4_VvI_q1 = np.array(rep4_VvI.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
rep4_VvI_q3 = np.array(rep4_VvI.quantile(q=[0.25,0.75], axis=1).iloc[1,:])


positions = np.arange(-250,1250,10)
ticks = np.array([-250, 0, 1000, 1250])
labels = ['-0.25 kb', 'TSS', 'TES', '+0.25 kb']
ybottom = -0.1
ytop = 2
yticks = np.array((0, 0.5, 1.0, 1.5))
replicates = [rep2_VvI, rep3_VvI, rep4_VvI]
names = ['rep2', 'rep3', 'rep4']
counter = 0

for df in replicates:
    avg = np.array(df.mean(1))
    q1 = np.array(df.quantile(q=[0.25,0.75], axis=1).iloc[0,:])
    q3 = np.array(df.quantile(q=[0.25,0.75], axis=1).iloc[1,:])
    plt.figure(dpi=300)
    plt.plot(positions, avg[150:300], label='Spt6')
    plt.fill_between(positions, q1[150:300], q3[150:300], alpha=0.2)
    plt.plot(positions, avg[450:600], label='Spt6Δ2-238')
    plt.fill_between(positions, q1[450:600], q3[450:600], alpha=0.2)
    plt.plot(positions, avg[0:150], label='Spt6 DMSO')
    plt.fill_between(positions, q1[0:150], q3[0:150], alpha=0.2)
    plt.plot(positions, avg[300:450], label='Spt6Δ2-238 DMSO')
    plt.fill_between(positions, q1[300:450], q3[300:450], alpha=0.2)
    plt.xlabel('position')
    plt.ylabel('V5 enrichment')
    plt.xlim(left=-250, right=1250)
    plt.xticks(ticks=ticks, labels=labels)
    plt.ylim(ybottom, ytop)
    plt.yticks(ticks=yticks)
    plt.legend(loc='upper left')
    plt.title(names[counter])
    plt.savefig('V5_'+names[counter]+'.png')
    plt.show();
    counter+=1

### Scatterplots

# Flag
rep2_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep2_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep2_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
# V5
rep2_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)


rep2_toplot = pd.DataFrame()
rep2_toplot['Spt6_Flag'] = rep2_FvI.iloc[:,150:300].mean(1)
rep2_toplot['Spt6Δ2-238_Flag'] = rep2_FvI.iloc[:,450:600].mean(1)
rep2_toplot['Spt6_8WG16'] = rep2_RvI.iloc[:,150:300].mean(1)
rep2_toplot['Spt6Δ2-238_8WG16'] = rep2_RvI.iloc[:,450:600].mean(1)
rep2_toplot['Spt6_Flagv8WG16'] = rep2_FvR.iloc[:,150:300].mean(1)
rep2_toplot['Spt6Δ2-238_Flagv8WG16'] = rep2_FvR.iloc[:,450:600].mean(1)
rep2_toplot['Spt6_V5_DMSO'] = rep2_VvI.iloc[:,0:150].mean(1)
rep2_toplot['Spt6_V5'] = rep2_VvI.iloc[:,150:300].mean(1)
rep2_toplot['Spt6Δ2-238_V5_DMSO'] = rep2_VvI.iloc[:,300:450].mean(1)
rep2_toplot['Spt6Δ2-238_V5'] = rep2_VvI.iloc[:,450:600].mean(1)

rep3_toplot = pd.DataFrame()
rep3_toplot['Spt6_Flag'] = rep3_FvI.iloc[:,150:300].mean(1)
rep3_toplot['Spt6Δ2-238_Flag'] = rep3_FvI.iloc[:,450:600].mean(1)
rep3_toplot['Spt6_8WG16'] = rep3_RvI.iloc[:,150:300].mean(1)
rep3_toplot['Spt6Δ2-238_8WG16'] = rep3_RvI.iloc[:,450:600].mean(1)
rep3_toplot['Spt6_Flagv8WG16'] = rep3_FvR.iloc[:,150:300].mean(1)
rep3_toplot['Spt6Δ2-238_Flagv8WG16'] = rep3_FvR.iloc[:,450:600].mean(1)
rep3_toplot['Spt6_V5_DMSO'] = rep3_VvI.iloc[:,0:150].mean(1)
rep3_toplot['Spt6_V5'] = rep3_VvI.iloc[:,150:300].mean(1)
rep3_toplot['Spt6Δ2-238_V5_DMSO'] = rep3_VvI.iloc[:,300:450].mean(1)
rep3_toplot['Spt6Δ2-238_V5'] = rep3_VvI.iloc[:,450:600].mean(1)

rep4_toplot = pd.DataFrame()
rep4_toplot['Spt6_Flag'] = rep4_FvI.iloc[:,150:300].mean(1)
rep4_toplot['Spt6Δ2-238_Flag'] = rep4_FvI.iloc[:,450:600].mean(1)
rep4_toplot['Spt6_8WG16'] = rep4_RvI.iloc[:,150:300].mean(1)
rep4_toplot['Spt6Δ2-238_8WG16'] = rep4_RvI.iloc[:,450:600].mean(1)
rep4_toplot['Spt6_Flagv8WG16'] = rep4_FvR.iloc[:,150:300].mean(1)
rep4_toplot['Spt6Δ2-238_Flagv8WG16'] = rep4_FvR.iloc[:,450:600].mean(1)
rep4_toplot['Spt6_V5_DMSO'] = rep4_VvI.iloc[:,0:150].mean(1)
rep4_toplot['Spt6_V5'] = rep4_VvI.iloc[:,150:300].mean(1)
rep4_toplot['Spt6Δ2-238_V5_DMSO'] = rep4_VvI.iloc[:,300:450].mean(1)
rep4_toplot['Spt6Δ2-238_V5'] = rep4_VvI.iloc[:,450:600].mean(1)


# # Rep 2
# plt.figure(figsize=(4,4), dpi=300)
# plt.plot([-0.5,13],[-0.5,13], color='black', linestyle='--')
# sns.scatterplot(data=rep2_toplot, x='Spt6_Flag', y='Spt6Δ2-238_Flag')
# plt.xlim(-0.5,13)
# plt.xlabel('Spt6')
# plt.ylabel('Spt6Δ2-238')
# plt.ylim(-0.5,13)
# plt.title('Flag rep2')
# plt.show();

# plt.figure(figsize=(4,4), dpi=300)
# plt.plot([-0.5,25],[-0.5,25], color='black', linestyle='--')
# sns.scatterplot(data=rep2_toplot, x='Spt6_8WG16', y='Spt6Δ2-238_8WG16')
# plt.xlim(-0.5,25)
# plt.xlabel('Spt6')
# plt.ylabel('Spt6Δ2-238')
# plt.ylim(-0.5,25)
# plt.title('Rpb1 rep2')
# plt.show();


# # Rep 3
# plt.figure(figsize=(4,4), dpi=300)
# plt.plot([-0.5,13],[-0.5,13], color='black', linestyle='--')
# sns.scatterplot(data=rep3_toplot, x='Spt6_Flag', y='Spt6Δ2-238_Flag')
# #sns.kdeplot(data=rep3_toplot, x='Spt6_Flag', y='Spt6Δ2-238_Flag', color='red')
# plt.xlim(-0.5,13)
# plt.xlabel('Spt6')
# plt.ylabel('Spt6Δ2-238')
# plt.ylim(-0.5,13)
# plt.title('Flag rep3')
# plt.show();

# plt.figure(figsize=(4,4), dpi=300)
# plt.plot([-0.5,25],[-0.5,25], color='black', linestyle='--')
# sns.scatterplot(data=rep3_toplot, x='Spt6_8WG16', y='Spt6Δ2-238_8WG16')
# #sns.kdeplot(data=rep3_toplot, x='Spt6_8WG16', y='Spt6Δ2-238_8WG16', color='red')
# plt.xlim(-0.5,25)
# plt.xlabel('Spt6')
# plt.ylabel('Spt6Δ2-238')
# plt.ylim(-0.5,25)
# plt.title('Rpb1 rep3')
# plt.show();

# # Rep 4
# plt.figure(figsize=(4,4), dpi=300)
# plt.plot([-0.5,13],[-0.5,13], color='black', linestyle='--')
# sns.scatterplot(data=rep4_toplot, x='Spt6_Flag', y='Spt6Δ2-238_Flag')
# #sns.kdeplot(data=rep4_toplot, x='Spt6_Flag', y='Spt6Δ2-238_Flag', color='red')
# plt.xlim(-0.5,13)
# plt.xlabel('Spt6')
# plt.ylabel('Spt6Δ2-238')
# plt.ylim(-0.5,13)
# plt.title('Flag rep4')
# plt.show();

# plt.figure(figsize=(4,4), dpi=300)
# plt.plot([-0.5,25],[-0.5,25], color='black', linestyle='--')
# sns.scatterplot(data=rep4_toplot, x='Spt6_8WG16', y='Spt6Δ2-238_8WG16')
# #sns.kdeplot(data=rep4_toplot, x='Spt6_8WG16', y='Spt6Δ2-238_8WG16', color='red')
# plt.xlim(-0.5,25)
# plt.xlabel('Spt6')
# plt.ylabel('Spt6Δ2-238')
# plt.ylim(-0.5,25)
# plt.title('Rpb1 rep4')
# plt.show();


### Testing showing density of points

#Rep 2
plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep2_toplot['Spt6_Flag'], rep2_toplot['Spt6Δ2-238_Flag']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,13],[-0.5,13], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep2_toplot, x='Spt6_Flag', y='Spt6Δ2-238_Flag', 
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(-0.5,13)
plt.xlabel('Spt6')
plt.ylabel('Spt6Δ2-238')
plt.ylim(-0.5,13)
plt.title('Flag rep2')
plt.savefig('Flag_rep2_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep2_toplot['Spt6_8WG16'], rep2_toplot['Spt6Δ2-238_8WG16']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,25],[-0.5,25], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep2_toplot, x='Spt6_8WG16', y='Spt6Δ2-238_8WG16', 
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(-0.5,25)
plt.xlabel('Spt6')
plt.ylabel('Spt6Δ2-238')
plt.ylim(-0.5,25)
plt.title('Rpb1 rep2')
plt.savefig('8WG16_rep2_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep2_toplot['Spt6_Flagv8WG16'], rep2_toplot['Spt6Δ2-238_Flagv8WG16']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,1.5],[-0.5,1.5], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep2_toplot, x='Spt6_Flagv8WG16', y='Spt6Δ2-238_Flagv8WG16',
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,1.25)
plt.xlabel('Spt6')
plt.ylabel('Spt6Δ2-238')
plt.ylim(0,1.25)
plt.title('Flag/Rpb1 rep2')
plt.savefig('Flag_v_8WG16_rep2_scatter.png')
plt.show();

# V5
plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep2_toplot['Spt6_V5_DMSO'], rep2_toplot['Spt6_V5']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,15],[-0.5,15], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep2_toplot, x='Spt6_V5_DMSO', y='Spt6_V5',
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,12)
plt.xlabel('DMSO')
plt.ylabel('IAA')
plt.ylim(0,12)
plt.title('Spt6 V5 rep2')
plt.savefig('Spt6_V5_rep2_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep2_toplot['Spt6Δ2-238_V5_DMSO'], rep2_toplot['Spt6Δ2-238_V5']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,15],[-0.5,15], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep2_toplot, x='Spt6Δ2-238_V5_DMSO', y='Spt6Δ2-238_V5',
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,12)
plt.xlabel('DMSO')
plt.ylabel('IAA')
plt.ylim(0,12)
plt.title('Spt6Δ2-238 V5 rep2')
plt.savefig('Spt6D2-238_V5_rep2_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep3_toplot['Spt6_V5_DMSO'], rep3_toplot['Spt6_V5']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,15],[-0.5,15], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep3_toplot, x='Spt6_V5_DMSO', y='Spt6_V5',
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,12)
plt.xlabel('DMSO')
plt.ylabel('IAA')
plt.ylim(0,12)
plt.title('Spt6 V5 rep3')
plt.savefig('Spt6_V5_rep3_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep3_toplot['Spt6Δ2-238_V5_DMSO'], rep3_toplot['Spt6Δ2-238_V5']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,15],[-0.5,15], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep3_toplot, x='Spt6Δ2-238_V5_DMSO', y='Spt6Δ2-238_V5',
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,12)
plt.xlabel('DMSO')
plt.ylabel('IAA')
plt.ylim(0,12)
plt.title('Spt6Δ2-238 V5 rep3')
plt.savefig('Spt6D2-238_V5_rep3_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep4_toplot['Spt6_V5_DMSO'], rep4_toplot['Spt6_V5']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,15],[-0.5,15], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep4_toplot, x='Spt6_V5_DMSO', y='Spt6_V5',
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,12)
plt.xlabel('DMSO')
plt.ylabel('IAA')
plt.ylim(0,12)
plt.title('Spt6 V5 rep4')
plt.savefig('Spt6_V5_rep4_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep4_toplot['Spt6Δ2-238_V5_DMSO'], rep4_toplot['Spt6Δ2-238_V5']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,15],[-0.5,15], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep4_toplot, x='Spt6Δ2-238_V5_DMSO', y='Spt6Δ2-238_V5',
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,12)
plt.xlabel('DMSO')
plt.ylabel('IAA')
plt.ylim(0,12)
plt.title('Spt6Δ2-238 V5 rep4')
plt.savefig('Spt6D2-238_V5_rep4_scatter.png')
plt.show();


# Rep 3
plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep3_toplot['Spt6_Flag'], rep3_toplot['Spt6Δ2-238_Flag']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,13],[-0.5,13], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep3_toplot, x='Spt6_Flag', y='Spt6Δ2-238_Flag', 
                s=15, c=kernel, cmap='inferno', zorder=2)
plt.xlim(-0.5,13)
plt.xticks(fontsize=16)
plt.xlabel('Spt6', fontsize=18)
plt.ylabel('Spt6Δ2-238', fontsize=18)
plt.yticks(fontsize=16)
plt.ylim(-0.5,13)
plt.title('$\\alpha\mathdefault{-Flag}$', fontsize=20, pad=10)
plt.savefig('Flag_rep3_scatter.svg')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep3_toplot['Spt6_8WG16'], rep3_toplot['Spt6Δ2-238_8WG16']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,25],[-0.5,25], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep3_toplot, x='Spt6_8WG16', y='Spt6Δ2-238_8WG16', 
                s=15, c=kernel, cmap='inferno', zorder=2)
plt.xlim(-0.5,25)
plt.xticks(fontsize=16)
plt.xlabel('Spt6', fontsize=18)
plt.ylabel('Spt6Δ2-238', fontsize=18)
plt.yticks(fontsize=16)
plt.ylim(-0.5,25)
plt.title('$\\alpha\mathdefault{-Rpb1}$', fontsize=20, pad=10)
plt.savefig('8WG16_rep3_scatter.svg')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep3_toplot['Spt6_Flagv8WG16'], rep3_toplot['Spt6Δ2-238_Flagv8WG16']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,13],[-0.5,13], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep3_toplot, x='Spt6_Flagv8WG16', y='Spt6Δ2-238_Flagv8WG16',
                s=15, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,1.25)
plt.xticks(fontsize=16)
plt.xlabel('Spt6', fontsize=18)
plt.ylabel('Spt6Δ2-238', fontsize=18)
plt.yticks(fontsize=16)
plt.ylim(0,1.25)
plt.title('$\\frac{\\alpha\mathdefault{-Flag}}{\\alpha\mathdefault{-Rpb1}}$', fontsize=26, pad=20)
plt.savefig('Flag_v_8WG16_rep3_scatter.svg')
plt.show();


#Rep 4
plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep4_toplot['Spt6_Flag'], rep4_toplot['Spt6Δ2-238_Flag']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,13],[-0.5,13], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep4_toplot, x='Spt6_Flag', y='Spt6Δ2-238_Flag', 
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(-0.5,13)
plt.xlabel('Spt6')
plt.ylabel('Spt6Δ2-238')
plt.ylim(-0.5,13)
plt.title('Flag rep4')
plt.savefig('Flag_rep4_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep4_toplot['Spt6_8WG16'], rep4_toplot['Spt6Δ2-238_8WG16']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,25],[-0.5,25], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep4_toplot, x='Spt6_8WG16', y='Spt6Δ2-238_8WG16', 
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(-0.5,25)
plt.xlabel('Spt6')
plt.ylabel('Spt6Δ2-238')
plt.ylim(-0.5,25)
plt.title('Rpb1 rep4')
plt.savefig('8WG16_rep4_scatter.png')
plt.show();

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([rep4_toplot['Spt6_Flagv8WG16'], rep4_toplot['Spt6Δ2-238_Flagv8WG16']])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,13],[-0.5,13], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=rep4_toplot, x='Spt6_Flagv8WG16', y='Spt6Δ2-238_Flagv8WG16',
                s=7, c=kernel, cmap='inferno', zorder=2)
plt.xlim(0,1.25)
plt.xlabel('Spt6')
plt.ylabel('Spt6Δ2-238')
plt.ylim(0,1.25)
plt.title('Flag/Rpb1 rep4')
plt.savefig('Flag_v_8WG16_rep4_scatter.png')
plt.show();


### Testing boxplots

rep2_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep2_si_test_ratio_scale.tab', 
                  sep='\t', header=3)
rep3_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep3_si_test_ratio_scale.tab', 
                  sep='\t', header=3)
rep4_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep4_si_test_ratio_scale.tab', 
                  sep='\t', header=3)
rep2_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep2_si_test_ratio_scale.tab', 
                  sep='\t', header=3)
rep3_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep3_si_test_ratio_scale.tab', 
                  sep='\t', header=3)
rep4_RvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/8WG16vinput_rep4_si_test_ratio_scale.tab', 
                  sep='\t', header=3)
# rep2_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep2_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep3_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep3_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# rep4_FvR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagv8WG16_rep4_si_test_ratio_scale.tab', 
#                  sep='\t', header=3)
# V5
rep2_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_VvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/V5vinput_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)

rep2_toplot = pd.DataFrame()
rep2_toplot['Spt6_Flag_DMSO'] = rep2_FvI.iloc[:,0:150].mean(1)
rep2_toplot['Spt6_Flag'] = rep2_FvI.iloc[:,150:300].mean(1)
rep2_toplot['Spt6Δ2-238_Flag_DMSO'] = rep2_FvI.iloc[:,300:450].mean(1)
rep2_toplot['Spt6Δ2-238_Flag'] = rep2_FvI.iloc[:,450:600].mean(1)
rep2_toplot['Spt6_8WG16_DMSO'] = rep2_RvI.iloc[:,0:150].mean(1)
rep2_toplot['Spt6_8WG16'] = rep2_RvI.iloc[:,150:300].mean(1)
rep2_toplot['Spt6Δ2-238_8WG16_DMSO'] = rep2_RvI.iloc[:,300:450].mean(1)
rep2_toplot['Spt6Δ2-238_8WG16'] = rep2_RvI.iloc[:,450:600].mean(1)
# rep2_toplot['Spt6_Flagv8WG16'] = rep2_FvR.iloc[:,150:300].mean(1)
# rep2_toplot['Spt6Δ2-238_Flagv8WG16'] = rep2_FvR.iloc[:,450:600].mean(1)
rep2_toplot['Spt6_V5_DMSO'] = rep2_VvI.iloc[:,0:150].mean(1)
rep2_toplot['Spt6_V5'] = rep2_VvI.iloc[:,150:300].mean(1)
rep2_toplot['Spt6Δ2-238_V5_DMSO'] = rep2_VvI.iloc[:,300:450].mean(1)
rep2_toplot['Spt6Δ2-238_V5'] = rep2_VvI.iloc[:,450:600].mean(1)

rep3_toplot = pd.DataFrame()
rep3_toplot['Spt6_Flag_DMSO'] = rep3_FvI.iloc[:,0:150].mean(1)
rep3_toplot['Spt6_Flag'] = rep3_FvI.iloc[:,150:300].mean(1)
rep3_toplot['Spt6Δ2-238_Flag_DMSO'] = rep3_FvI.iloc[:,300:450].mean(1)
rep3_toplot['Spt6Δ2-238_Flag'] = rep3_FvI.iloc[:,450:600].mean(1)
rep3_toplot['Spt6_8WG16_DMSO'] = rep3_RvI.iloc[:,0:150].mean(1)
rep3_toplot['Spt6_8WG16'] = rep3_RvI.iloc[:,150:300].mean(1)
rep3_toplot['Spt6Δ2-238_8WG16_DMSO'] = rep3_RvI.iloc[:,300:450].mean(1)
rep3_toplot['Spt6Δ2-238_8WG16'] = rep3_RvI.iloc[:,450:600].mean(1)
# rep3_toplot['Spt6_Flagv8WG16'] = rep3_FvR.iloc[:,150:300].mean(1)
# rep3_toplot['Spt6Δ2-238_Flagv8WG16'] = rep3_FvR.iloc[:,450:600].mean(1)
rep3_toplot['Spt6_V5_DMSO'] = rep3_VvI.iloc[:,0:150].mean(1)
rep3_toplot['Spt6_V5'] = rep3_VvI.iloc[:,150:300].mean(1)
rep3_toplot['Spt6Δ2-238_V5_DMSO'] = rep3_VvI.iloc[:,300:450].mean(1)
rep3_toplot['Spt6Δ2-238_V5'] = rep3_VvI.iloc[:,450:600].mean(1)

rep4_toplot = pd.DataFrame()
rep4_toplot['Spt6_Flag_DMSO'] = rep4_FvI.iloc[:,0:150].mean(1)
rep4_toplot['Spt6_Flag'] = rep4_FvI.iloc[:,150:300].mean(1)
rep4_toplot['Spt6Δ2-238_Flag_DMSO'] = rep4_FvI.iloc[:,300:450].mean(1)
rep4_toplot['Spt6Δ2-238_Flag'] = rep4_FvI.iloc[:,450:600].mean(1)
rep4_toplot['Spt6_8WG16_DMSO'] = rep4_RvI.iloc[:,0:150].mean(1)
rep4_toplot['Spt6_8WG16'] = rep4_RvI.iloc[:,150:300].mean(1)
rep4_toplot['Spt6Δ2-238_8WG16_DMSO'] = rep4_RvI.iloc[:,300:450].mean(1)
rep4_toplot['Spt6Δ2-238_8WG16'] = rep4_RvI.iloc[:,450:600].mean(1)
# rep4_toplot['Spt6_Flagv8WG16'] = rep4_FvR.iloc[:,150:300].mean(1)
# rep4_toplot['Spt6Δ2-238_Flagv8WG16'] = rep4_FvR.iloc[:,450:600].mean(1)
rep4_toplot['Spt6_V5_DMSO'] = rep4_VvI.iloc[:,0:150].mean(1)
rep4_toplot['Spt6_V5'] = rep4_VvI.iloc[:,150:300].mean(1)
rep4_toplot['Spt6Δ2-238_V5_DMSO'] = rep4_VvI.iloc[:,300:450].mean(1)
rep4_toplot['Spt6Δ2-238_V5'] = rep4_VvI.iloc[:,450:600].mean(1)


# generate a dataframe that has the sample condition as a column
rep2_toplot.melt(value_vars=['Spt6_V5_DMSO', 'Spt6_V5', 'Spt6Δ2-238_V5_DMSO', 'Spt6Δ2-238_V5'])

sns.set_style('ticks')

# V5
plt.figure(dpi=300)
lala = sns.boxplot(data=rep2_toplot.melt(value_vars=['Spt6_V5_DMSO', 'Spt6_V5', 'Spt6Δ2-238_V5_DMSO', 'Spt6Δ2-238_V5']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='V5 coverage')
lala.set(ylim=(0,2.5))
lala.set(title='rep2')
plt.show();

plt.figure(dpi=300)
lala = sns.boxplot(data=rep3_toplot.melt(value_vars=['Spt6_V5_DMSO', 'Spt6_V5', 'Spt6Δ2-238_V5_DMSO', 'Spt6Δ2-238_V5']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='V5 coverage')
lala.set(ylim=(0,2.5))
lala.set(title='rep3')
plt.show();

plt.figure(dpi=300)
lala = sns.boxplot(data=rep4_toplot.melt(value_vars=['Spt6_V5_DMSO', 'Spt6_V5', 'Spt6Δ2-238_V5_DMSO', 'Spt6Δ2-238_V5']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='V5 coverage')
lala.set(ylim=(0,2.5))
lala.set(title='rep4')
plt.show();

# Flag
plt.figure(dpi=300)
lala = sns.boxplot(data=rep2_toplot.melt(value_vars=['Spt6_Flag_DMSO', 'Spt6_Flag', 'Spt6Δ2-238_Flag_DMSO', 'Spt6Δ2-238_Flag']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='Flag coverage')
lala.set(ylim=(0,3))
lala.set(title='rep2')
plt.show();

plt.figure(dpi=300)
lala = sns.boxplot(data=rep3_toplot.melt(value_vars=['Spt6_Flag_DMSO', 'Spt6_Flag', 'Spt6Δ2-238_Flag_DMSO', 'Spt6Δ2-238_Flag']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='Flag coverage')
lala.set(ylim=(0,3))
lala.set(title='rep3')
plt.show();

plt.figure(dpi=300)
lala = sns.boxplot(data=rep4_toplot.melt(value_vars=['Spt6_Flag_DMSO', 'Spt6_Flag', 'Spt6Δ2-238_Flag_DMSO', 'Spt6Δ2-238_Flag']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='Flag coverage')
lala.set(ylim=(0,3))
lala.set(title='rep4')
plt.show();

# 8WG16
plt.figure(dpi=300)
lala = sns.boxplot(data=rep2_toplot.melt(value_vars=['Spt6_8WG16_DMSO', 'Spt6_8WG16', 'Spt6Δ2-238_8WG16_DMSO', 'Spt6Δ2-238_8WG16']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='8WG16 coverage')
lala.set(ylim=(0,4))
lala.set(title='rep2')
plt.show();

plt.figure(dpi=300)
lala = sns.boxplot(data=rep3_toplot.melt(value_vars=['Spt6_8WG16_DMSO', 'Spt6_8WG16', 'Spt6Δ2-238_8WG16_DMSO', 'Spt6Δ2-238_8WG16']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='8WG16 coverage')
lala.set(ylim=(0,4))
lala.set(title='rep3')
plt.show();

plt.figure(dpi=300)
lala = sns.boxplot(data=rep4_toplot.melt(value_vars=['Spt6_8WG16_DMSO', 'Spt6_8WG16', 'Spt6Δ2-238_8WG16_DMSO', 'Spt6Δ2-238_8WG16']),
            x = 'variable', y='value', fliersize=0, linewidth=2, color='green')
lala.set(xlabel=None)
lala.set(xticklabels=['Spt6\nDMSO','Spt6\nIAA','Spt6Δ2-238\nDMSO','Spt6Δ2-238\nIAA'])
lala.set(ylabel='8WG16 coverage')
lala.set(ylim=(0,4))
lala.set(title='rep4')
plt.show();




# To do:
# Improve graph formatting.  
# Confirm exclusion of outlier points.
# Add Flag DMSO columns and generate a Flag plot.
# Question/Goal: Illustrate the outcompetiton of Flag by V5 in the Spt6 mutant.




    
### Testing heatmaps  
# Flag
rep2_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep2_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep3_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep3_si_test_ratio_scale.tab', 
                 sep='\t', header=3)
rep4_FvI = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/deeptools/tab/Flagvinput_rep4_si_test_ratio_scale.tab', 
                 sep='\t', header=3)



plt.figure()
lala = sns.heatmap(rep2_FvI.iloc[:,150:300])
lala.set(xticklabels=[],yticklabels=[])
plt.show();

#Trying to sort genes in descending order of signal
test = rep2_FvI.iloc[:,150:300]
test['mean'] = test.mean(1)

#sorting, re-indedex
data = test.sort_values(by=['mean'],axis=0,ascending=False,ignore_index=True)

#remove highest 10% of genes
data = data.loc[300:,:]

plt.figure()
lala = sns.heatmap(data)
lala.set(xticklabels=[],yticklabels=[])
plt.show();

    
    