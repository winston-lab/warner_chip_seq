#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 11:16:30 2025

@author: jlwarner
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

sns.set_style('ticks', {'axes.facecolor':'white'})

#Import Data
positions = np.arange(-250,1250,10)

# samples - Rpb1
wt_93_D_Rpb1 = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_D_8WG16vinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_93_I_Rpb1 = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_I_8WG16vinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_D_Rpb1 = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_D_8WG16vinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_I_Rpb1 = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_I_8WG16vinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

# samples - flag
wt_93_D_flag = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_D_Flagvinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_93_I_flag = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_I_Flagvinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_D_flag = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_D_Flagvinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_I_flag = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_I_Flagvinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

# samples - V5
wt_93_D_v5 = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_D_V5vinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_93_I_v5 = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_I_V5vinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_D_v5 = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_D_V5vinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_I_v5 = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_I_V5vinput_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

# samples - flag/Rpb1
wt_93_D_fpR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_D_Flagv8WG16_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_93_I_fpR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_I_Flagv8WG16_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_D_fpR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_D_Flagv8WG16_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_I_fpR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_I_Flagv8WG16_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

# samples - V5/Rpb1
wt_93_D_vpR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_D_V5v8WG16_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_93_I_vpR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/93_I_V5v8WG16_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_D_vpR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_D_V5v8WG16_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)
wt_95_I_vpR = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged_midpoints/95_I_V5v8WG16_averaged_midpoint_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)


### Metagenes
ticks = np.array([-250, 0, 1000, 1250])
labels = ['-250 bp', 'TSS', 'CPS', '+250 bp']
ybottom = -0.1
ytop = 5.0
yticks = np.array((0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))

#Rpb1 metagenes
plt.figure(dpi=300)
plt.plot(positions, wt_93_D_Rpb1.mean(), label = 'WT DMSO', color = 'black')
#plt.fill_between(positions, wt_93_D_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_93_D_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#                 alpha=0.15, color='black')
plt.plot(positions, wt_93_I_Rpb1.mean(), label = 'WT IAA', color = 'red')
plt.fill_between(positions, wt_93_I_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_93_I_Rpb1.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='red')
plt.plot(positions, wt_95_D_Rpb1.mean(), label = 'ΔNTD DMSO', color = 'orange')
#plt.fill_between(positions, wt_95_D_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_95_D_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#                 alpha=0.15, color='orange')
plt.plot(positions, wt_95_I_Rpb1.mean(), label = 'ΔNTD IAA', color = 'blue')
plt.fill_between(positions, wt_95_I_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_95_I_Rpb1.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized Rpb1 occupancy')
#plt.ylim(ybottom, ytop)
#plt.yticks(ticks=yticks)
plt.legend(loc='upper left')
plt.title('Rpb1 ChIP averaged')
plt.tight_layout()
#plt.savefig('Rpb1_mean_metagene_averaged_midpoints.png', dpi=300)
plt.show();

# plt.figure(dpi=300)
# plt.plot(positions, wt_93_D_Rpb1.median(), label = 'WT DMSO', color = 'black')
# plt.fill_between(positions, wt_93_D_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_93_D_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#                  alpha=0.15, color='black')
# plt.plot(positions, wt_93_I_Rpb1.median(), label = 'WT IAA', color = 'red')
# plt.fill_between(positions, wt_93_I_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_93_I_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#                  alpha=0.15, color='red')
# plt.plot(positions, wt_95_D_Rpb1.median(), label = 'ΔNTD DMSO', color = 'orange')
# plt.fill_between(positions, wt_95_D_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_95_D_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#                  alpha=0.15, color='orange')
# plt.plot(positions, wt_95_I_Rpb1.median(), label = 'ΔNTD IAA', color = 'blue')
# plt.fill_between(positions, wt_95_I_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_95_I_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#                  alpha=0.15, color='blue')
# plt.xlabel('position')
# plt.xlim(left=-250, right=1250)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylabel('spike-in normalized Rpb1 occupancy')
# #plt.ylim(ybottom, ytop)
# #plt.yticks(ticks=yticks)
# plt.legend(loc='upper left')
# plt.title('Rpb1 ChIP averaged')
# plt.tight_layout()
# #plt.savefig('Rpb1_median_metagene_averaged_midpoints.png', dpi=300)
# plt.show();


#Flag metagenes
plt.figure(dpi=300)
plt.plot(positions, wt_93_D_flag.mean(), label = 'WT DMSO', color = 'black')
#plt.fill_between(positions, wt_93_D_flag.quantile([0.25,0.75]).iloc[0,:], wt_93_D_flag.quantile([0.25,0.75]).iloc[1,:],
#                 alpha=0.15, color='black')
plt.plot(positions, wt_93_I_flag.mean(), label = 'WT IAA', color = 'red')
plt.fill_between(positions, wt_93_I_flag.quantile([0.25,0.75]).iloc[0,:], wt_93_I_flag.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='red')
plt.plot(positions, wt_95_D_flag.mean(), label = 'ΔNTD DMSO', color = 'orange')
#plt.fill_between(positions, wt_95_D_flag.quantile([0.25,0.75]).iloc[0,:], wt_95_D_flag.quantile([0.25,0.75]).iloc[1,:],
#                 alpha=0.15, color='orange')
plt.plot(positions, wt_95_I_flag.mean(), label = 'ΔNTD IAA', color = 'blue')
plt.fill_between(positions, wt_95_I_flag.quantile([0.25,0.75]).iloc[0,:], wt_95_I_flag.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized Flag occupancy')
#plt.ylim(ybottom, ytop)
#plt.yticks(ticks=yticks)
plt.legend(loc='upper left')
plt.title('Flag ChIP averaged')
plt.tight_layout()
#plt.savefig('flag_mean_metagene_averaged_midpoints.png', dpi=300)
plt.show();

# plt.figure(dpi=300)
# plt.plot(positions, spt6_pob3_30_myc.median(), label = 'SPT6 POB3', color = 'black')
# #plt.fill_between(positions, spt6_pob3_30_myc.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_30_myc.quantile([0.25,0.75]).iloc[1,:],
#  #                alpha=0.2, color='black')
# plt.plot(positions, yw_pob3_30_myc.median(), label = 'spt6-YW POB3', color = 'red')
# #plt.fill_between(positions, yw_pob3_30_myc.quantile([0.25,0.75]).iloc[0,:], yw_pob3_30_myc.quantile([0.25,0.75]).iloc[1,:],
# #                 alpha=0.2, color='red')
# plt.plot(positions, spt6_e154k_30_myc.median(), label = 'SPT6 E154K', color = 'orange')
# #plt.fill_between(positions, spt6_e154k_30_myc.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_30_myc.quantile([0.25,0.75]).iloc[1,:],
# #                 alpha=0.2, color='red')
# plt.plot(positions, yw_e154k_30_myc.median(), label = 'spt6-YW E154K', color = 'blue')
# #plt.fill_between(positions, yw_e154k_30_myc.quantile([0.25,0.75]).iloc[0,:], yw_e154k_30_myc.quantile([0.25,0.75]).iloc[1,:],
# #                 alpha=0.2, color='red')
# plt.xlabel('position')
# plt.xlim(left=-250, right=1250)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylabel('spike-in normalized myc occupancy')
# #plt.ylim(ybottom, ytop)
# #plt.yticks(ticks=yticks)
# plt.legend(loc='upper left')
# plt.title('30°C median midpoints_rep2')
# #plt.savefig('Myc-30_median_metagene_midpoints_rep2.png', dpi=300)
# plt.show();


#V5 metagenes
plt.figure(dpi=300)
plt.plot(positions, wt_93_D_v5.mean(), label = 'WT DMSO', color = 'black')
plt.fill_between(positions, wt_93_D_v5.quantile([0.25,0.75]).iloc[0,:], wt_93_D_v5.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='black')
plt.plot(positions, wt_93_I_v5.mean(), label = 'WT IAA', color = 'red')
plt.fill_between(positions, wt_93_I_v5.quantile([0.25,0.75]).iloc[0,:], wt_93_I_v5.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='red')
plt.plot(positions, wt_95_D_v5.mean(), label = 'ΔNTD DMSO', color = 'orange')
plt.fill_between(positions, wt_95_D_v5.quantile([0.25,0.75]).iloc[0,:], wt_95_D_v5.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='orange')
plt.plot(positions, wt_95_I_v5.mean(), label = 'ΔNTD IAA', color = 'blue')
plt.fill_between(positions, wt_95_I_v5.quantile([0.25,0.75]).iloc[0,:], wt_95_I_v5.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized V5 occupancy')
#plt.ylim(ybottom, ytop)
#plt.yticks(ticks=yticks)
plt.legend(loc='upper left')
plt.title('V5 ChIP averaged')
plt.tight_layout()
#plt.savefig('V5_mean_metagene_averaged_midpoints.png', dpi=300)
plt.show();

#Flag / Rpb1 metagenes
plt.figure(dpi=300)
plt.plot(positions, wt_93_D_fpR.mean(), label = 'WT DMSO', color = 'black')
#plt.fill_between(positions, wt_93_D_fpR.quantile([0.25,0.75]).iloc[0,:], wt_93_D_fpR.quantile([0.25,0.75]).iloc[1,:],
#                 alpha=0.15, color='black')
plt.plot(positions, wt_93_I_fpR.mean(), label = 'WT IAA', color = 'red')
plt.fill_between(positions, wt_93_I_fpR.quantile([0.25,0.75]).iloc[0,:], wt_93_I_fpR.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='red')
plt.plot(positions, wt_95_D_fpR.mean(), label = 'ΔNTD DMSO', color = 'orange')
#plt.fill_between(positions, wt_95_D_fpR.quantile([0.25,0.75]).iloc[0,:], wt_95_D_fpR.quantile([0.25,0.75]).iloc[1,:],
#                 alpha=0.15, color='orange')
plt.plot(positions, wt_95_I_fpR.mean(), label = 'ΔNTD IAA', color = 'blue')
plt.fill_between(positions, wt_95_I_fpR.quantile([0.25,0.75]).iloc[0,:], wt_95_I_fpR.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized $\\frac{Flag}{Rpb1}$ occupancy')
#plt.ylim(ybottom, ytop)
#plt.yticks(ticks=yticks)
plt.legend(loc='upper left')
plt.title('$\\frac{Flag}{Rpb1}$ ChIP averaged')
plt.tight_layout()
#plt.savefig('fpR_mean_metagene_averaged_midpoints.png', dpi=300)
plt.show();

#V5 / Rpb1 metagenes
plt.figure(dpi=300)
plt.plot(positions, wt_93_D_vpR.mean(), label = 'WT DMSO', color = 'black')
plt.fill_between(positions, wt_93_D_vpR.quantile([0.25,0.75]).iloc[0,:], wt_93_D_vpR.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='black')
plt.plot(positions, wt_93_I_vpR.mean(), label = 'WT IAA', color = 'red')
plt.fill_between(positions, wt_93_I_vpR.quantile([0.25,0.75]).iloc[0,:], wt_93_I_vpR.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='red')
plt.plot(positions, wt_95_D_vpR.mean(), label = 'ΔNTD DMSO', color = 'orange')
plt.fill_between(positions, wt_95_D_vpR.quantile([0.25,0.75]).iloc[0,:], wt_95_D_vpR.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='orange')
plt.plot(positions, wt_95_I_vpR.mean(), label = 'ΔNTD IAA', color = 'blue')
plt.fill_between(positions, wt_95_I_vpR.quantile([0.25,0.75]).iloc[0,:], wt_95_I_vpR.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.15, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized $\\frac{V5}{Rpb1}$ occupancy')
#plt.ylim(ybottom, ytop)
#plt.yticks(ticks=yticks)
plt.legend(loc='upper left')
plt.title('$\\frac{V5}{Rpb1}$ ChIP averaged')
plt.tight_layout()
#plt.savefig('vpR_mean_metagene_averaged_midpoints.png', dpi=300)
plt.show();

# plt.figure(dpi=300)
# plt.plot(positions, spt6_pob3_30_mpR.median(), label = 'SPT6 POB3', color = 'black')
# #plt.fill_between(positions, spt6_pob3_30_mpR.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#  #                alpha=0.2, color='black')
# plt.plot(positions, yw_pob3_30_mpR.median(), label = 'spt6-YW POB3', color = 'red')
# #plt.fill_between(positions, yw_pob3_30_mpR.quantile([0.25,0.75]).iloc[0,:], yw_pob3_30_mpR.quantile([0.25,0.75]).iloc[1,:],
# #                 alpha=0.2, color='red')
# plt.plot(positions, spt6_e154k_30_mpR.median(), label = 'SPT6 E154K', color = 'orange')
# #plt.fill_between(positions, spt6_e154k_30_mpR.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_30_mpR.quantile([0.25,0.75]).iloc[1,:],
# #                 alpha=0.2, color='red')
# plt.plot(positions, yw_e154k_30_mpR.median(), label = 'spt6-YW E154K', color = 'blue')
# #plt.fill_between(positions, yw_e154k_30_mpR.quantile([0.25,0.75]).iloc[0,:], yw_e154k_30_mpR.quantile([0.25,0.75]).iloc[1,:],
# #                 alpha=0.2, color='red')
# plt.xlabel('position')
# plt.xlim(left=-250, right=1250)
# plt.xticks(ticks=ticks, labels=labels)
# plt.ylabel('spike-in normalized $\\frac{Myc}{Rpb1}$ occupancy')
# #plt.ylim(ybottom, ytop)
# #plt.yticks(ticks=yticks)
# plt.legend(loc='upper left')
# plt.title('30°C median midpoints_rep2')
# #plt.savefig('MpR-30_median_metagene_midpoints_rep2.png', dpi=300)
# plt.show();




# ## Testing deciles mean vs. median

# # Rpb1
# for (start_tile, end_tile) in [(x/10, (x+1)/10) for x in range(10)]:
#     plt.figure(dpi=300)
#     plt.plot(positions, spt6_pob3_30_Rpb1.iloc[int(spt6_pob3_30_Rpb1.shape[0]*start_tile):int(spt6_pob3_30_Rpb1.shape[0]*end_tile),:].mean(), label = 'SPT6 POB3', color = 'black')
#     #plt.fill_between(positions, spt6_pob3_30_Rpb1.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_30_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#       #                alpha=0.2, color='black')
#     plt.plot(positions, yw_pob3_30_Rpb1.iloc[int(yw_pob3_30_Rpb1.shape[0]*start_tile):int(yw_pob3_30_Rpb1.shape[0]*end_tile),:].mean(), label = 'spt6-YW POB3', color = 'red')
#     #plt.fill_between(positions, yw_pob3_30_Rpb1.quantile([0.25,0.75]).iloc[0,:], yw_pob3_30_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, spt6_e154k_30_Rpb1.iloc[int(spt6_e154k_30_Rpb1.shape[0]*start_tile):int(spt6_e154k_30_Rpb1.shape[0]*end_tile),:].mean(), label = 'SPT6 E154K', color = 'orange')
#     #plt.fill_between(positions, spt6_e154k_30_Rpb1.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_30_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, yw_e154k_30_Rpb1.iloc[int(yw_e154k_30_Rpb1.shape[0]*start_tile):int(yw_e154k_30_Rpb1.shape[0]*end_tile),:].mean(), label = 'spt6-YW E154K', color = 'blue')
#     #plt.fill_between(positions, yw_e154k_30_Rpb1.quantile([0.25,0.75]).iloc[0,:], yw_e154k_30_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.xlabel('position')
#     plt.xlim(left=-250, right=1250)
#     plt.xticks(ticks=ticks, labels=labels)
#     plt.ylabel('spike-in normalized Rpb1 occupancy')
#     #plt.ylim(ybottom, ytop)
#     #plt.yticks(ticks=yticks)
#     plt.legend(loc='upper left')
#     plt.title(str(start_tile)+' to '+str(end_tile)+' (mean) 30°C midpoints_rep2')
#     ##plt.savefig('Rpb1-30_metagene.png', dpi=300)
#     plt.show();
    
#     plt.figure(dpi=300)
#     plt.plot(positions, spt6_pob3_30_Rpb1.iloc[int(spt6_pob3_30_Rpb1.shape[0]*start_tile):int(spt6_pob3_30_Rpb1.shape[0]*end_tile),:].median(), label = 'SPT6 POB3', color = 'black')
#     #plt.fill_between(positions, spt6_pob3_30_Rpb1.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_30_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#       #                alpha=0.2, color='black')
#     plt.plot(positions, yw_pob3_30_Rpb1.iloc[int(yw_pob3_30_Rpb1.shape[0]*start_tile):int(yw_pob3_30_Rpb1.shape[0]*end_tile),:].median(), label = 'spt6-YW POB3', color = 'red')
#     #plt.fill_between(positions, yw_pob3_30_Rpb1.quantile([0.25,0.75]).iloc[0,:], yw_pob3_30_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, spt6_e154k_30_Rpb1.iloc[int(spt6_e154k_30_Rpb1.shape[0]*start_tile):int(spt6_e154k_30_Rpb1.shape[0]*end_tile),:].median(), label = 'SPT6 E154K', color = 'orange')
#     #plt.fill_between(positions, spt6_e154k_30_Rpb1.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_30_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, yw_e154k_30_Rpb1.iloc[int(yw_e154k_30_Rpb1.shape[0]*start_tile):int(yw_e154k_30_Rpb1.shape[0]*end_tile),:].median(), label = 'spt6-YW E154K', color = 'blue')
#     #plt.fill_between(positions, yw_e154k_30_Rpb1.quantile([0.25,0.75]).iloc[0,:], yw_e154k_30_Rpb1.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.xlabel('position')
#     plt.xlim(left=-250, right=1250)
#     plt.xticks(ticks=ticks, labels=labels)
#     plt.ylabel('spike-in normalized Rpb1 occupancy')
#     #plt.ylim(ybottom, ytop)
#     #plt.yticks(ticks=yticks)
#     plt.legend(loc='upper left')
#     plt.title(str(start_tile)+' to '+str(end_tile)+' (median) 30°C midpoints_rep2')
#     ##plt.savefig('Rpb1-30_metagene.png', dpi=300)
#     plt.show();
    
# # Myc 
# for (start_tile, end_tile) in [(x/10, (x+1)/10) for x in range(10)]:
#     plt.figure(dpi=300)
#     plt.plot(positions, spt6_pob3_30_myc.iloc[int(spt6_pob3_30_myc.shape[0]*start_tile):int(spt6_pob3_30_myc.shape[0]*end_tile),:].mean(), label = 'SPT6 POB3', color = 'black')
#     #plt.fill_between(positions, spt6_pob3_30_myc.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_30_myc.quantile([0.25,0.75]).iloc[1,:],
#       #                alpha=0.2, color='black')
#     plt.plot(positions, yw_pob3_30_myc.iloc[int(yw_pob3_30_myc.shape[0]*start_tile):int(yw_pob3_30_myc.shape[0]*end_tile),:].mean(), label = 'spt6-YW POB3', color = 'red')
#     #plt.fill_between(positions, yw_pob3_30_myc.quantile([0.25,0.75]).iloc[0,:], yw_pob3_30_myc.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, spt6_e154k_30_myc.iloc[int(spt6_e154k_30_myc.shape[0]*start_tile):int(spt6_e154k_30_myc.shape[0]*end_tile),:].mean(), label = 'SPT6 E154K', color = 'orange')
#     #plt.fill_between(positions, spt6_e154k_30_myc.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_30_myc.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, yw_e154k_30_myc.iloc[int(yw_e154k_30_myc.shape[0]*start_tile):int(yw_e154k_30_myc.shape[0]*end_tile),:].mean(), label = 'spt6-YW E154K', color = 'blue')
#     #plt.fill_between(positions, yw_e154k_30_myc.quantile([0.25,0.75]).iloc[0,:], yw_e154k_30_myc.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.xlabel('position')
#     plt.xlim(left=-250, right=1250)
#     plt.xticks(ticks=ticks, labels=labels)
#     plt.ylabel('spike-in normalized myc occupancy')
#     #plt.ylim(ybottom, ytop)
#     #plt.yticks(ticks=yticks)
#     plt.legend(loc='upper left')
#     plt.title(str(start_tile)+' to '+str(end_tile)+' (mean) 30°C midpoints_rep2')
#     ##plt.savefig('myc-30_metagene.png', dpi=300)
#     plt.show();
    
#     plt.figure(dpi=300)
#     plt.plot(positions, spt6_pob3_30_myc.iloc[int(spt6_pob3_30_myc.shape[0]*start_tile):int(spt6_pob3_30_myc.shape[0]*end_tile),:].median(), label = 'SPT6 POB3', color = 'black')
#     #plt.fill_between(positions, spt6_pob3_30_myc.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_30_myc.quantile([0.25,0.75]).iloc[1,:],
#       #                alpha=0.2, color='black')
#     plt.plot(positions, yw_pob3_30_myc.iloc[int(yw_pob3_30_myc.shape[0]*start_tile):int(yw_pob3_30_myc.shape[0]*end_tile),:].median(), label = 'spt6-YW POB3', color = 'red')
#     #plt.fill_between(positions, yw_pob3_30_myc.quantile([0.25,0.75]).iloc[0,:], yw_pob3_30_myc.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, spt6_e154k_30_myc.iloc[int(spt6_e154k_30_myc.shape[0]*start_tile):int(spt6_e154k_30_myc.shape[0]*end_tile),:].median(), label = 'SPT6 E154K', color = 'orange')
#     #plt.fill_between(positions, spt6_e154k_30_myc.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_30_myc.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, yw_e154k_30_myc.iloc[int(yw_e154k_30_myc.shape[0]*start_tile):int(yw_e154k_30_myc.shape[0]*end_tile),:].median(), label = 'spt6-YW E154K', color = 'blue')
#     #plt.fill_between(positions, yw_e154k_30_myc.quantile([0.25,0.75]).iloc[0,:], yw_e154k_30_myc.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.xlabel('position')
#     plt.xlim(left=-250, right=1250)
#     plt.xticks(ticks=ticks, labels=labels)
#     plt.ylabel('spike-in normalized myc occupancy')
#     #plt.ylim(ybottom, ytop)
#     #plt.yticks(ticks=yticks)
#     plt.legend(loc='upper left')
#     plt.title(str(start_tile)+' to '+str(end_tile)+' (median) 30°C midpoints_rep2')
#     ##plt.savefig('myc-30_metagene.png', dpi=300)
#     plt.show();


# # Myc/Rpb1
# for (start_tile, end_tile) in [(x/10, (x+1)/10) for x in range(10)]:
#     plt.figure(dpi=300)
#     plt.plot(positions, spt6_pob3_30_mpR.iloc[int(spt6_pob3_30_mpR.shape[0]*start_tile):int(spt6_pob3_30_mpR.shape[0]*end_tile),:].mean(), label = 'SPT6 POB3', color = 'black')
#     #plt.fill_between(positions, spt6_pob3_30_mpR.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#       #                alpha=0.2, color='black')
#     plt.plot(positions, yw_pob3_30_mpR.iloc[int(yw_pob3_30_mpR.shape[0]*start_tile):int(yw_pob3_30_mpR.shape[0]*end_tile),:].mean(), label = 'spt6-YW POB3', color = 'red')
#     #plt.fill_between(positions, yw_pob3_30_mpR.quantile([0.25,0.75]).iloc[0,:], yw_pob3_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, spt6_e154k_30_mpR.iloc[int(spt6_e154k_30_mpR.shape[0]*start_tile):int(spt6_e154k_30_mpR.shape[0]*end_tile),:].mean(), label = 'SPT6 E154K', color = 'orange')
#     #plt.fill_between(positions, spt6_e154k_30_mpR.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, yw_e154k_30_mpR.iloc[int(yw_e154k_30_mpR.shape[0]*start_tile):int(yw_e154k_30_mpR.shape[0]*end_tile),:].mean(), label = 'spt6-YW E154K', color = 'blue')
#     #plt.fill_between(positions, yw_e154k_30_mpR.quantile([0.25,0.75]).iloc[0,:], yw_e154k_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.xlabel('position')
#     plt.xlim(left=-250, right=1250)
#     plt.xticks(ticks=ticks, labels=labels)
#     plt.ylabel('spike-in normalized $\\frac{myc}{Rpb1}$ occupancy')
#     #plt.ylim(ybottom, ytop)
#     #plt.yticks(ticks=yticks)
#     plt.legend(loc='upper left')
#     plt.title(str(start_tile)+' to '+str(end_tile)+' (mean) 30°C midpoints_rep2')
#     ##plt.savefig('mpR-30_metagene.png', dpi=300)
#     plt.show();
    
#     plt.figure(dpi=300)
#     plt.plot(positions, spt6_pob3_30_mpR.iloc[int(spt6_pob3_30_mpR.shape[0]*start_tile):int(spt6_pob3_30_mpR.shape[0]*end_tile),:].median(), label = 'SPT6 POB3', color = 'black')
#     #plt.fill_between(positions, spt6_pob3_30_mpR.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#       #                alpha=0.2, color='black')
#     plt.plot(positions, yw_pob3_30_mpR.iloc[int(yw_pob3_30_mpR.shape[0]*start_tile):int(yw_pob3_30_mpR.shape[0]*end_tile),:].median(), label = 'spt6-YW POB3', color = 'red')
#     #plt.fill_between(positions, yw_pob3_30_mpR.quantile([0.25,0.75]).iloc[0,:], yw_pob3_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, spt6_e154k_30_mpR.iloc[int(spt6_e154k_30_mpR.shape[0]*start_tile):int(spt6_e154k_30_mpR.shape[0]*end_tile),:].median(), label = 'SPT6 E154K', color = 'orange')
#     #plt.fill_between(positions, spt6_e154k_30_mpR.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.plot(positions, yw_e154k_30_mpR.iloc[int(yw_e154k_30_mpR.shape[0]*start_tile):int(yw_e154k_30_mpR.shape[0]*end_tile),:].median(), label = 'spt6-YW E154K', color = 'blue')
#     #plt.fill_between(positions, yw_e154k_30_mpR.quantile([0.25,0.75]).iloc[0,:], yw_e154k_30_mpR.quantile([0.25,0.75]).iloc[1,:],
#     #                 alpha=0.2, color='red')
#     plt.xlabel('position')
#     plt.xlim(left=-250, right=1250)
#     plt.xticks(ticks=ticks, labels=labels)
#     plt.ylabel('spike-in normalized $\\frac{myc}{Rpb1}$ occupancy')
#     #plt.ylim(ybottom, ytop)
#     #plt.yticks(ticks=yticks)
#     plt.legend(loc='upper left')
#     plt.title(str(start_tile)+' to '+str(end_tile)+' (median) 30°C midpoints_rep2')
#     ##plt.savefig('mpR-30_metagene.png', dpi=300)
#     plt.show();

## Violin Plots
# Rpb1
wt_D = np.log2(wt_93_D_Rpb1.iloc[:,25:125].mean(1))
wt_I = np.log2((wt_93_I_Rpb1.iloc[:,25:125].mean(1)))
mut_D = np.log2((wt_95_D_Rpb1.iloc[:,25:125].mean(1)))
mut_I = np.log2((wt_95_I_Rpb1.iloc[:,25:125].mean(1)))

df = pd.DataFrame({
#    'WT DMSO':wt_D, 
    'WT IAA':wt_I, 
#    'ΔNTD DMSO':mut_D, 
    'ΔNTD IAA':mut_I
    })
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 Rpb1 occupancy', fontsize=14)
plt.title('Rpb1 ChIP per gene', fontsize=16, pad=10)
plt.tight_layout()
#plt.savefig('Rpb1_violin_30_midpoints_rep2.png', dpi=300)
plt.show();

# Flag
wt_D = np.log2(wt_93_D_flag.iloc[:,25:125].mean(1))
wt_I = np.log2((wt_93_I_flag.iloc[:,25:125].mean(1)))
mut_D = np.log2((wt_95_D_flag.iloc[:,25:125].mean(1)))
mut_I = np.log2((wt_95_I_flag.iloc[:,25:125].mean(1)))

df = pd.DataFrame({
#    'WT DMSO':wt_D, 
    'WT IAA':wt_I, 
#    'ΔNTD DMSO':mut_D, 
    'ΔNTD IAA':mut_I
    })
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 Flag occupancy', fontsize=14)
plt.title('Flag ChIP per gene', fontsize=16, pad=10)
plt.tight_layout()
#plt.savefig('flag_violin_30_midpoints_rep2.png', dpi=300)
plt.show();

# V5
wt_D = np.log2(wt_93_D_v5.iloc[:,25:125].mean(1))
wt_I = np.log2((wt_93_I_v5.iloc[:,25:125].mean(1)))
mut_D = np.log2((wt_95_D_v5.iloc[:,25:125].mean(1)))
mut_I = np.log2((wt_95_I_v5.iloc[:,25:125].mean(1)))

df = pd.DataFrame({
    'WT DMSO':wt_D, 
    'WT IAA':wt_I, 
    'ΔNTD DMSO':mut_D, 
    'ΔNTD IAA':mut_I
    })
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 V5 occupancy', fontsize=14)
plt.title('V5 ChIP per gene', fontsize=16, pad=10)
plt.tight_layout()
#plt.savefig('v5_violin_30_midpoints_rep2.png', dpi=300)
plt.show();

# Flag / Rpb1
wt_D = np.log2(wt_93_D_fpR.iloc[:,25:125].mean(1))
wt_I = np.log2((wt_93_I_fpR.iloc[:,25:125].mean(1)))
mut_D = np.log2((wt_95_D_fpR.iloc[:,25:125].mean(1)))
mut_I = np.log2((wt_95_I_fpR.iloc[:,25:125].mean(1)))

df = pd.DataFrame({
#    'WT DMSO':wt_D, 
    'WT IAA':wt_I, 
#    'ΔNTD DMSO':mut_D, 
    'ΔNTD IAA':mut_I
    })
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 $\\frac{Flag}{Rpb1}$ occupancy', fontsize=14)
plt.title('$\\frac{Flag}{Rpb1}$ ChIP per gene', fontsize=16, pad=10)
plt.tight_layout()
#plt.savefig('fpR_violin_30_midpoints_rep2.png', dpi=300)
plt.show();

# V5 / Rpb1
wt_D = np.log2(wt_93_D_vpR.iloc[:,25:125].mean(1))
wt_I = np.log2((wt_93_I_vpR.iloc[:,25:125].mean(1)))
mut_D = np.log2((wt_95_D_vpR.iloc[:,25:125].mean(1)))
mut_I = np.log2((wt_95_I_vpR.iloc[:,25:125].mean(1)))

df = pd.DataFrame({
    'WT DMSO':wt_D, 
    'WT IAA':wt_I, 
    'ΔNTD DMSO':mut_D, 
    'ΔNTD IAA':mut_I
    })
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 $\\frac{V5}{Rpb1}$ occupancy', fontsize=14)
plt.title('$\\frac{V5}{Rpb1}$ ChIP per gene', fontsize=16, pad=10)
plt.tight_layout()
#plt.savefig('vpR_violin_30_midpoints_rep2.png', dpi=300)
plt.show();




#Scatterplots
# x = spt6_pob3_30_Rpb1.iloc[:,25:125].mean(1)
# y = yw_pob3_30_Rpb1.iloc[:,25:125].mean(1)
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-0.5,60.5],[-0.5,60.5], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-0.5,60.5)
# plt.xticks(np.arange(0,60.5,10), fontsize=16)
# plt.xlabel('Spt6', fontsize=18)
# plt.ylabel('Spt6-YW', fontsize=18)
# plt.yticks(np.arange(0,60.5,10), fontsize=16)
# plt.ylim(-0.5,60.5)
# plt.title('Rpb1 midpoints_rep2', fontsize=20, pad=10)
# plt.tight_layout()
# ##plt.savefig('Rpb1_midpoints_rep2_scatter.png')
# plt.show();


# x = spt6_pob3_30_myc.iloc[:,25:125].mean(1)
# y = yw_pob3_30_myc.iloc[:,25:125].mean(1)
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-0.5,70.5],[-0.5,70.5], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-0.5,70.5)
# plt.xticks(np.arange(0,70.5,10), fontsize=16)
# plt.xlabel('Spt6', fontsize=18)
# plt.ylabel('Spt6-YW', fontsize=18)
# plt.yticks(np.arange(0,70.5,10), fontsize=16)
# plt.ylim(-0.5,70.5)
# plt.title('Myc midpoints_rep2', fontsize=20, pad=10)
# plt.tight_layout()
# ##plt.savefig('Myc_midpoints_rep2_scatter.png')
# plt.show();


# x = spt6_pob3_30_mpR.iloc[:,25:125].mean(1)
# y = yw_pob3_30_mpR.iloc[:,25:125].mean(1)
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-0.5,4.5],[-0.5,4.5], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-0.5,4.5)
# plt.xticks(np.arange(0,4.5,1), fontsize=16)
# plt.xlabel('$SPT6 \: POB3$', fontsize=18)
# plt.ylabel('$YW \: POB3$', fontsize=18)
# plt.yticks(np.arange(0,4.5,1), fontsize=16)
# plt.ylim(-0.5,4.5)
# plt.title('$\\frac{Myc}{Rpb1}$ midpoints_rep2', fontsize=20, pad=10)
# plt.tight_layout()
# ##plt.savefig('Myc_midpoints_rep2_scatter.png')
# plt.show();

# x = spt6_pob3_30_mpR.iloc[:,25:125].mean(1)
# y = spt6_e154k_30_mpR.iloc[:,25:125].mean(1)
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-0.5,4.5],[-0.5,4.5], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-0.5,4.5)
# plt.xticks(np.arange(0,4.5,1), fontsize=16)
# plt.xlabel('$SPT6 \: POB3$', fontsize=18)
# plt.ylabel('$SPT6 \: E154K$', fontsize=18)
# plt.yticks(np.arange(0,4.5,1), fontsize=16)
# plt.ylim(-0.5,4.5)
# plt.title('$\\frac{Myc}{Rpb1}$ midpoints_rep2', fontsize=20, pad=10)
# plt.tight_layout()
# ##plt.savefig('Myc_midpoints_rep2_scatter.png')
# plt.show();
#
# x = spt6_pob3_30_mpR.iloc[:,25:125].mean(1)
# y = yw_e154k_30_mpR.iloc[:,25:125].mean(1)
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-0.5,4.5],[-0.5,4.5], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-0.5,4.5)
# plt.xticks(np.arange(0,4.5,1), fontsize=16)
# plt.xlabel('$SPT6 \: POB3$', fontsize=18)
# plt.ylabel('$YW \: E154K$', fontsize=18)
# plt.yticks(np.arange(0,4.5,1), fontsize=16)
# plt.ylim(-0.5,4.5)
# plt.title('$\\frac{Myc}{Rpb1}$ midpoints_rep2', fontsize=20, pad=10)
# plt.tight_layout()
# ##plt.savefig('Myc_midpoints_rep2_scatter.png')
# plt.show();


x = np.log2(wt_93_I_Rpb1.iloc[:,25:125].mean(1))
y = np.log2((wt_95_I_flag.iloc[:,25:125].mean(1)/wt_93_I_flag.iloc[:,25:125].mean(1)))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-1.5,6.5],[0,0], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-2.5,5)
plt.xticks(np.arange(-2,6,1), fontsize=12)
plt.xlabel('WT IAA:\n Rpb1 enrichment', fontsize=14)
plt.ylabel('log2 $\\frac{ΔNTD}{WT}$\nFlag enrichment', fontsize=14)
plt.text(3.9, 1.8, str((y>0).sum()))
plt.text(3.9, -1.9, str((y<0).sum()))
#plt.yticks(np.arange(0,4.5,1), fontsize=16)
plt.ylim(-2,2)
#plt.title('ΔNTD vs WT', fontsize=16, pad=10)
plt.tight_layout()
#plt.savefig('log2fc_yw-wt_scatter_30_midpoints_rep2.png', dpi=300)
plt.show();


# x = np.log2(spt6_pob3_30_Rpb1.iloc[:,25:125].mean(1))
# y = np.log2((spt6_e154k_30_myc.iloc[:,25:125].mean(1)/spt6_pob3_30_myc.iloc[:,25:125].mean(1)))
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-1.5,6.5],[0,0], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-1.5,6)
# plt.xticks(np.arange(-1.5,6.5,1.5), fontsize=12)
# plt.xlabel('SPT6 POB3:\n Rpb1 enrichment', fontsize=14)
# plt.ylabel('log2 $\\frac{SPT6 \: pob3-E154K}{wild-type}$\nMyc enrichment', fontsize=14)
# plt.text(4.9, 4.1, str((y>0).sum()))
# plt.text(4.9, -1.8, str((y<0).sum()))
# #plt.yticks(np.arange(0,4.5,1), fontsize=16)
# plt.ylim(-2,4.5)
# plt.title('30°C midpoints_rep2', fontsize=16, pad=10)
# plt.tight_layout()
# #plt.savefig('log2fc_e154k-wt_scatter_30_midpoints_rep2.png', dpi=300)
# plt.show();

# x = np.log2(spt6_pob3_30_Rpb1.iloc[:,25:125].mean(1))
# y = np.log2((yw_e154k_30_myc.iloc[:,25:125].mean(1)/spt6_pob3_30_myc.iloc[:,25:125].mean(1)))
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-1.5,6.5],[0,0], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-1.5,6)
# plt.xticks(np.arange(-1.5,6.5,1.5), fontsize=12)
# plt.xlabel('SPT6 POB3:\n Rpb1 enrichment', fontsize=14)
# plt.ylabel('log2 $\\frac{spt6-YW \: pob3-E154K}{wild-type}$\nMyc enrichment', fontsize=14)
# plt.text(4.9, 4.1, str((y>0).sum()))
# plt.text(4.9, -1.8, str((y<0).sum()))
# #plt.yticks(np.arange(0,4.5,1), fontsize=16)
# plt.ylim(-2,4.5)
# plt.title('30°C midpoints_rep2', fontsize=16, pad=10)
# plt.tight_layout()
# #plt.savefig('log2fc_ywe154k-wt_scatter_30_midpoints_rep2.png', dpi=300)
# plt.show();


# # Violin plots of Rpb1 log2 fc
# yw = np.log2((yw_pob3_30_Rpb1.iloc[:,25:125].mean(1)/spt6_pob3_30_Rpb1.iloc[:,25:125].mean(1)))
# e154k = np.log2((spt6_e154k_30_Rpb1.iloc[:,25:125].mean(1)/spt6_pob3_30_Rpb1.iloc[:,25:125].mean(1)))
# yw_e154k = np.log2((yw_e154k_30_Rpb1.iloc[:,25:125].mean(1)/spt6_pob3_30_Rpb1.iloc[:,25:125].mean(1)))

# df = pd.DataFrame({'spt6-YW':yw, 'pob3-E154K':e154k, 'spt6-YW\npob3-E154K':yw_e154k})
# plt.figure(figsize=(4,3), dpi=300)
# sns.violinplot(data=df,
#                linecolor='black')
# plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
# plt.ylabel('log2 $\\frac{mutant}{wild-type}$ Rpb1', fontsize=14)
# plt.title('30°C midpoints_rep2', fontsize=16, pad=10)
# plt.tight_layout()
# #plt.savefig('log2fc_Rpb1_violin_30_midpoints_rep2.png', dpi=300)
# plt.show();

# # Violin plots of myc log2 fc
# yw = np.log2((yw_pob3_30_myc.iloc[:,25:125].mean(1)/spt6_pob3_30_myc.iloc[:,25:125].mean(1)))
# e154k = np.log2((spt6_e154k_30_myc.iloc[:,25:125].mean(1)/spt6_pob3_30_myc.iloc[:,25:125].mean(1)))
# yw_e154k = np.log2((yw_e154k_30_myc.iloc[:,25:125].mean(1)/spt6_pob3_30_myc.iloc[:,25:125].mean(1)))

# df = pd.DataFrame({'spt6-YW':yw, 'pob3-E154K':e154k, 'spt6-YW\npob3-E154K':yw_e154k})
# plt.figure(figsize=(4,3), dpi=300)
# sns.violinplot(data=df,
#                linecolor='black')
# plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
# plt.ylabel('log2 $\\frac{mutant}{wild-type}$ Myc', fontsize=14)
# plt.title('30°C midpoints_rep2', fontsize=16, pad=10)
# plt.tight_layout()
# #plt.savefig('log2fc_Myc_violin_30_midpoints_rep2.png', dpi=300)
# plt.show();

# # Violin plots of myc/Rpb1 log2 fc
# yw = np.log2((yw_pob3_30_mpR.iloc[:,25:125].mean(1)/spt6_pob3_30_mpR.iloc[:,25:125].mean(1)))
# e154k = np.log2((spt6_e154k_30_mpR.iloc[:,25:125].mean(1)/spt6_pob3_30_mpR.iloc[:,25:125].mean(1)))
# yw_e154k = np.log2((yw_e154k_30_mpR.iloc[:,25:125].mean(1)/spt6_pob3_30_mpR.iloc[:,25:125].mean(1)))

# df = pd.DataFrame({'spt6-YW':yw, 'pob3-E154K':e154k, 'spt6-YW\npob3-E154K':yw_e154k})
# plt.figure(figsize=(4,3), dpi=300)
# sns.violinplot(data=df,
#                linecolor='black')
# plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
# plt.ylabel('log2 $\\frac{mutant}{wild-type} \: \\frac{Myc}{Rpb1}$ ', fontsize=14)
# plt.title('30°C midpoints_rep2', fontsize=16, pad=10)
# plt.tight_layout()
# #plt.savefig('log2fc_MpR_violin_30_midpoints_rep2.png', dpi=300)
# plt.show();





#Heatmaps
ticks = [0,25,125,149]
labels = ['-250 bp', 'TSS', 'CPS', '+250 bp']

sns.set_style('ticks', {'axes.facecolor':'white'})

# plt.figure(figsize=(3,8), dpi=300)
# ax = sns.heatmap(data=wt_I_Rpb1, cmap=sns.color_palette("Greens", as_cmap=True), 
#             vmax=wt_I_Rpb1.quantile([0.90]).max(1), 
#             yticklabels=False)
# for spine in ax.spines.values():
#     spine.set(visible=True, lw=2, color='black')

# plt.ylabel('genes')
# plt.xticks(ticks=ticks, labels=labels)
# plt.title('Spt6')
# #plt.savefig("Rpb1_ChIP_WT.png")
# plt.show();

# Rpb1 Heatmaps
minimum = np.min([wt_93_I_Rpb1.quantile([0.10]).min(1),wt_95_I_Rpb1.quantile([0.10]).min(1),
                  ])
maximum = np.max([wt_93_I_Rpb1.quantile([0.90]).max(1),wt_95_I_Rpb1.quantile([0.90]).max(1),
                  ])
f = plt.figure(figsize=(6, 16), dpi=300)
gs = f.add_gridspec(1,2)
f.suptitle('Rpb1 ChIP occupancy\nmidpoints averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=wt_93_I_Rpb1, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('WT', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=wt_95_I_Rpb1, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('ΔNTD', fontsize=22)  

f.tight_layout()
#plt.savefig("Rpb1_heatmap_30_midpoints_rep2.png", dpi=300)
plt.show();

#colorbar
# a = np.array([[0,4]])
# plt.figure(figsize=(1,12))
# img = plt.imshow(a, cmap="Greens")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# ##plt.savefig('Rpb1_colorbar.svg', dpi=300)
# plt.show();

# Flag heatmaps
minimum = np.min([wt_93_I_flag.quantile([0.10]).min(1),wt_95_I_flag.quantile([0.10]).min(1),
                  ])
maximum = np.max([wt_93_I_flag.quantile([0.90]).max(1),wt_95_I_flag.quantile([0.90]).max(1),
                  ])
f = plt.figure(figsize=(6, 16), dpi=300)
gs = f.add_gridspec(1,2)
f.suptitle('Flag ChIP occupancy\nmidpoints averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=wt_93_I_flag, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('WT', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=wt_95_I_flag, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('ΔNTD', fontsize=22)  

# ax = f.add_subplot(gs[0,2])
# sns.heatmap(data=spt6_e154k_30_flag, cmap=sns.color_palette("Greens", as_cmap=True), 
#             vmax=maximum,
#             vmin=minimum, 
#             yticklabels=False, cbar=False)

# for spine in ax.spines.values():
#     spine.set(visible=True, lw=2, color='black')
# #plt.ylabel('genes', fontsize=20)
# plt.xticks(ticks=ticks, labels=labels, fontsize=20)
# plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
# plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
# plt.title('SPT6\npob3-E154K', fontsize=22)    

# ax = f.add_subplot(gs[0,3])
# sns.heatmap(data=yw_e154k_30_flag, cmap=sns.color_palette("Greens", as_cmap=True), 
#             vmax=maximum,
#             vmin=minimum, 
#             yticklabels=False, cbar=True)

# for spine in ax.spines.values():
#     spine.set(visible=True, lw=2, color='black')
# #plt.ylabel('genes', fontsize=20)
# plt.xticks(ticks=ticks, labels=labels, fontsize=20)
# plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
# plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
# plt.title('spt6-YW\npob3-E154K', fontsize=22)    

f.tight_layout()
#plt.savefig("flag_heatmap_30_midpoints_rep2.png", dpi=300)
plt.show();


#Flag / Rpb1 heatmap
minimum = np.min([wt_93_I_fpR.quantile([0.10]).min(1),wt_95_I_fpR.quantile([0.10]).min(1),
                  ])
maximum = np.max([wt_93_I_fpR.quantile([0.90]).max(1),wt_95_I_fpR.quantile([0.90]).max(1),
                  ])
f = plt.figure(figsize=(6, 16), dpi=300)
gs = f.add_gridspec(1,2)
f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy\nmidpoints averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=wt_93_I_fpR, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('WT', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=wt_95_I_fpR, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('ΔNTD', fontsize=22)  
 
f.tight_layout()
#plt.savefig("fpR_heatmap_30_midpoints_rep2.png", dpi=300)
plt.show();


#fold-change
fc_Rpb1 = np.log2(wt_95_I_Rpb1/wt_93_I_Rpb1)
fc_flag = np.log2(wt_95_I_flag/wt_93_I_flag)
fc_fpR = np.log2(wt_95_I_fpR/wt_93_I_fpR)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{ΔNTD}{WT}$ ChIP occupancy\nmidpoints averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_Rpb1, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('Rpb1', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=fc_flag, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('Flag', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=fc_fpR, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('$\\frac{Flag}{Rpb1}$', fontsize=22)   

f.tight_layout()
#plt.savefig("Rpb1_log2fc_heatmap_30_midpoints_rep2.png", dpi=300)
plt.show()




# fold-change


#colorbar
# a = np.array([[0,3.5]])
# plt.figure(figsize=(1,12))
# img = plt.imshow(a, cmap="Oranges")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical', ticks=[0,1,2,3,3.5])
# ##plt.savefig('Flag_colorbar.svg', dpi=300)
# plt.show();

#fold-change
# fc = np.log2(yw_pob3_30_myc/spt6_pob3_30_myc)

# plt.figure(figsize=(3,8), dpi=300)
# ax = sns.heatmap(data=fc, cmap=sns.color_palette("vlag", as_cmap=True),
#                  vmin=0-fc.quantile([0.95]).max(1), vmax=fc.quantile([0.95]).max(1),
#                  yticklabels=False)
# for spine in ax.spines.values():
#     spine.set(visible=True, lw=2, color='black')

# #plt.ylabel('genes')
# plt.xticks(ticks=ticks, labels=labels, fontsize=14)
# plt.title('Myc ChIP occupancy\nlog2$\\left(\\frac{spt6-YW}{SPT6}\\right)$', fontsize=16)
# plt.tight_layout()
# ##plt.savefig("Myc_heatmap_log2fc.png")
# plt.show();


# plt.figure(figsize=(2,4))
# lala = sns.boxplot(data=fc.iloc[:,25:125].mean(1), 
#             color='Orange', fliersize=0, linewidth=1, linecolor='black')
# lala.set(xlabel=None)
# lala.set(ylabel='log2$\\left(\\frac{spt6-YW}{SPT6}\\right)$ ChIP occupancy')
# lala.set(ylim=(-1,1))
# lala.set(title='Myc')
# plt.tight_layout()
# ##plt.savefig('Myc_boxplot_log2fc', dpi=300)
# plt.show();

# plt.figure(figsize=(2,4))
# lala = sns.violinplot(data=fc.iloc[:,25:125].mean(1), 
#             linewidth=1, inner= 'box',
#             linecolor='black', color = 'Orange',
#             inner_kws={'box_width':6, 'whis_width':2})
# lala.set(xlabel=None)
# lala.set(ylabel='log2$\\left(\\frac{spt6-YW}{SPT6}\\right)$ ChIP occupancy')
# #lala.set(ylim=(-1,1))
# lala.set(title='Flag')
# plt.tight_layout()
# ##plt.savefig('Myc_violin_log2fc', dpi=300)
# plt.show();


# #Violinplots
# wt = spt6_pob3_30_Rpb1.iloc[:,25:125].mean(1)
# mut = yw_pob3_30_Rpb1.iloc[:,25:125].mean(1)
# scipy.stats.wilcoxon(wt,mut)

# df_Rpb1 = pd.DataFrame({'Spt6':wt, 'Spt6-YW':mut})
# plt.figure(figsize=(4,3))
# sns.violinplot(data=df_Rpb1,
#                linecolor='black', color='Green')
# plt.ylabel('ChIP occupancy')
# plt.title('Rpb1')
# ##plt.savefig('Rpb1_violin.png', dpi=300)
# plt.show();


# plt.figure(figsize=(3,4.5))
# sns.violinplot(data=df_Rpb1, log_scale=2,
#                linecolor='black', color='Green')
# plt.ylabel('log2(ChIP occupancy)')
# plt.title('Rpb1')
# ##plt.savefig('Rpb1_log2_violin.png', dpi=300)
# plt.show();

# # df.to_csv('df.csv',index=False,header=True)


# wt = spt6_pob3_30_myc.iloc[:,25:125].mean(1)
# mut = yw_pob3_30_myc.iloc[:,25:125].mean(1)
# scipy.stats.wilcoxon(wt,mut)

# df_Flag = pd.DataFrame({'Spt6':wt, 'Spt6-YW':mut})
# plt.figure(figsize=(3,4.5))
# sns.violinplot(data=df_Flag, log_scale=2,
#                linecolor='black', color='Orange')
# plt.ylabel('log2(ChIP occupancy)')
# plt.title('Myc')
# ##plt.savefig('Flag_violin_log2.png', dpi=300)
# plt.show();

# # df.to_csv('df.csv',index=False,header=True)


# wt = (spt6_pob3_30_myc/spt6_pob3_30_Rpb1).mean(1)
# mut = (yw_pob3_30_myc/yw_pob3_30_Rpb1).mean(1)

# scipy.stats.wilcoxon(wt,mut)

# df_FlagvRpb1 = pd.DataFrame({'Spt6':wt, 'Spt6-YW':mut})
# plt.figure(figsize=(3,4.5))
# sns.violinplot(data=df_FlagvRpb1, 
#                log_scale=2,
#                linecolor='black')
# plt.ylabel('ChIP occupancy')
# plt.title('Myc/Rpb1')
# ##plt.savefig('MycvRpb1_violin.png', dpi=300)
# plt.show();

# #df_FlagvRpb1.to_csv('df.csv',index=False,header=True)
