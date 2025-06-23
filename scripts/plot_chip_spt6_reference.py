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
positions = np.arange(-250,4500,10)

#averaged samples
spt6_dmso_Rpb1_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/93_D_8WG16vinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_iaa_Rpb1_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/93_I_8WG16vinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_dmso_Rpb1_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/95_D_8WG16vinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_iaa_Rpb1_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/95_I_8WG16vinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_dmso_Flag_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/93_D_Flagvinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_iaa_Flag_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/93_I_Flagvinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_dmso_Flag_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/95_D_Flagvinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_iaa_Flag_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/95_I_Flagvinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_dmso_V5_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/93_D_V5vinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_iaa_V5_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/93_I_V5vinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_dmso_V5_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/95_D_V5vinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_iaa_V5_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/ratio/tab/95_I_V5vinput_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_dmso_FpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/perpol/tab/93_D_Flagv8WG16_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_iaa_FpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/perpol/tab/93_I_Flagv8WG16_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_dmso_FpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/perpol/tab/95_D_Flagv8WG16_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_iaa_FpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/perpol/tab/95_I_Flagv8WG16_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_dmso_VpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/perpol/tab/93_D_V5v8WG16_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_iaa_VpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/perpol/tab/93_I_V5v8WG16_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_dmso_VpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/perpol/tab/95_D_V5v8WG16_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)

dntd_iaa_VpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/averaged/perpol/tab/95_I_V5v8WG16_averaged_midpoints_si_ratio_reference.tab', 
            sep='\t', header=3, names=positions)



#Heatmaps
ticks = [25,125,225,325,425]
labels = ['TSS', '+1 kb', '+2 kb', '+3 kb', '+4 kb']

sns.set_style('ticks', {'axes.facecolor':'#DCDCDC'})

# Rpb1 Heatmaps
minimum = np.min([spt6_dmso_Rpb1_averaged.quantile([0.10]).min(1),dntd_dmso_Rpb1_averaged.quantile([0.10]).min(1),
                  spt6_iaa_Rpb1_averaged.quantile([0.10]).min(1),dntd_iaa_Rpb1_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_dmso_Rpb1_averaged.quantile([0.90]).max(1),dntd_dmso_Rpb1_averaged.quantile([0.90]).max(1),
                  spt6_iaa_Rpb1_averaged.quantile([0.90]).max(1),dntd_iaa_Rpb1_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_dmso_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nDMSO', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_dmso_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nDMSO', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_iaa_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nIAA', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_iaa_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nIAA', fontsize=22)    

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/Rpb1_heatmap_reference.png", dpi=300)
plt.show();

#colorbar
# a = np.array([[0,4]])
# plt.figure(figsize=(1,12))
# img = plt.imshow(a, cmap="Greens")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# ##plt.savefig('Rpb1_colorbar.svg', dpi=300)
# plt.show();


#fold-change
fc_0 = np.log2(dntd_dmso_Rpb1_averaged/spt6_dmso_Rpb1_averaged)
fc_1 = np.log2(spt6_iaa_Rpb1_averaged/spt6_dmso_Rpb1_averaged)
fc_2 = np.log2(dntd_iaa_Rpb1_averaged/spt6_dmso_Rpb1_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_0, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=fc_1, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=fc_2, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/Rpb1_log2fc_heatmap_reference.png", dpi=300)
plt.show()

fc_3 = np.log2(dntd_iaa_Rpb1_averaged/spt6_iaa_Rpb1_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/Rpb1_log2fc_dntd-e154kvsdntd_heatmap_reference.png", dpi=300)
plt.show();




# V5 Heatmaps
minimum = np.min([spt6_dmso_V5_averaged.quantile([0.10]).min(1),dntd_dmso_V5_averaged.quantile([0.10]).min(1),
                  spt6_iaa_V5_averaged.quantile([0.10]).min(1),dntd_iaa_V5_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_dmso_V5_averaged.quantile([0.90]).max(1),dntd_dmso_V5_averaged.quantile([0.90]).max(1),
                  spt6_iaa_V5_averaged.quantile([0.90]).max(1),dntd_iaa_V5_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('V5 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_dmso_V5_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_dmso_V5_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_iaa_V5_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_iaa_V5_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/V5_heatmap_reference.png", dpi=300)
plt.show();

#fold-change
fc_0 = np.log2(dntd_dmso_V5_averaged/spt6_dmso_V5_averaged)
fc_1 = np.log2(spt6_iaa_V5_averaged/spt6_dmso_V5_averaged)
fc_2 = np.log2(dntd_iaa_V5_averaged/spt6_dmso_V5_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -5
maximum = 5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ V5 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_0, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=fc_1, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=fc_2, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

# ax = f.add_subplot(gs[0,3])
# a = np.array([[minimum, maximum]])
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/V5_log2fc_heatmap_reference.png", dpi=300)
plt.show();

# a = np.array([[minimum, maximum]])
# plt.figure(figsize=(2,16), dpi=300)
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# plt.show();

fc_3 = np.log2(dntd_iaa_V5_averaged/dntd_dmso_V5_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('log2 $\\frac{spt6ΔNTD\:pob3-E154K}{spt6ΔNTD\:POB3}$ V5 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/V5_log2fc_dntd-e154kvsdntd_heatmap_reference.png", dpi=300)
plt.show();


# Flag Heatmaps
minimum = np.min([spt6_dmso_Flag_averaged.quantile([0.10]).min(1),dntd_dmso_Flag_averaged.quantile([0.10]).min(1),
                  spt6_iaa_Flag_averaged.quantile([0.10]).min(1),dntd_iaa_Flag_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_dmso_Flag_averaged.quantile([0.90]).max(1),dntd_dmso_Flag_averaged.quantile([0.90]).max(1),
                  spt6_iaa_Flag_averaged.quantile([0.90]).max(1),dntd_iaa_Flag_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('Flag ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_dmso_Flag_averaged, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_dmso_Flag_averaged, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_iaa_Flag_averaged, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_iaa_Flag_averaged, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/Flag_heatmap_reference.png", dpi=300)
plt.show();

#fold-change
fc_0 = np.log2(dntd_dmso_Flag_averaged/spt6_dmso_Flag_averaged)
fc_1 = np.log2(spt6_iaa_Flag_averaged/spt6_dmso_Flag_averaged)
fc_2 = np.log2(dntd_iaa_Flag_averaged/spt6_dmso_Flag_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ Flag ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_0, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=fc_1, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22, style='italic')  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=fc_2, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22, style='italic')    

# ax = f.add_subplot(gs[0,3])
# a = np.array([[minimum, maximum]])
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/Flag_log2fc_heatmap_reference.png", dpi=300)
plt.show();

fc_3 = np.log2(dntd_iaa_Flag_averaged/spt6_iaa_Flag_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('log2 $\\frac{spt6ΔNTD\:IAA}{SPT6\:IAA}$ Flag ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/Flag_log2fc_dntde154kvsdntd_heatmap_reference.png", dpi=300)
plt.show();



# VpR Heatmaps
minimum = np.min([spt6_dmso_VpR_averaged.quantile([0.10]).min(1),dntd_dmso_VpR_averaged.quantile([0.10]).min(1),
                  spt6_iaa_VpR_averaged.quantile([0.10]).min(1),dntd_iaa_VpR_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_dmso_VpR_averaged.quantile([0.90]).max(1),dntd_dmso_VpR_averaged.quantile([0.90]).max(1),
                  spt6_iaa_VpR_averaged.quantile([0.90]).max(1),dntd_iaa_VpR_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('$\\frac{V5}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_dmso_VpR_averaged, cmap=sns.color_palette("Reds", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_dmso_VpR_averaged, cmap=sns.color_palette("Reds", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_iaa_VpR_averaged, cmap=sns.color_palette("Reds", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_iaa_VpR_averaged, cmap=sns.color_palette("Reds", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/VpR_heatmap_reference.png", dpi=300)
plt.show();

#fold-change
fc_0 = np.log2(dntd_dmso_VpR_averaged/spt6_dmso_VpR_averaged)
fc_1 = np.log2(spt6_iaa_VpR_averaged/spt6_dmso_VpR_averaged)
fc_2 = np.log2(dntd_iaa_VpR_averaged/spt6_dmso_VpR_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ VpR ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_0, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=fc_1, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=fc_2, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

# ax = f.add_subplot(gs[0,3])
# a = np.array([[minimum, maximum]])
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/VpR_log2fc_heatmap_reference.png", dpi=300)
plt.show();

fc_3 = np.log2(dntd_iaa_VpR_averaged/dntd_dmso_VpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('log2 $\\frac{spt6ΔNTD\:IAA}{spt6ΔNTD\:DMSO}$ VpR ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/VpR_log2fc_dntd-e154kvsdntd_heatmap_reference.png", dpi=300)
plt.show();


# FpR Heatmaps
minimum = np.min([spt6_dmso_FpR_averaged.quantile([0.10]).min(1),dntd_dmso_FpR_averaged.quantile([0.10]).min(1),
                  spt6_iaa_FpR_averaged.quantile([0.10]).min(1),dntd_iaa_FpR_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_dmso_FpR_averaged.quantile([0.90]).max(1),dntd_dmso_FpR_averaged.quantile([0.90]).max(1),
                  spt6_iaa_FpR_averaged.quantile([0.90]).max(1),dntd_iaa_FpR_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_dmso_FpR_averaged, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_dmso_FpR_averaged, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_iaa_FpR_averaged, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_iaa_FpR_averaged, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/FpR_heatmap_reference.png", dpi=300)
plt.show();

#fold-change
fc_0 = np.log2(dntd_dmso_FpR_averaged/spt6_dmso_FpR_averaged)
fc_1 = np.log2(spt6_iaa_FpR_averaged/spt6_dmso_FpR_averaged)
fc_2 = np.log2(dntd_iaa_FpR_averaged/spt6_dmso_FpR_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ FpR ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_0, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=fc_1, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=fc_2, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

# ax = f.add_subplot(gs[0,3])
# a = np.array([[minimum, maximum]])
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/FpR_log2fc_heatmap_reference.png", dpi=300)
plt.show();

fc_3 = np.log2(dntd_iaa_FpR_averaged/spt6_iaa_FpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ FpR ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/FpR_log2fc_dntde154kvsdntd_heatmap_reference.png", dpi=300)
plt.show();


#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,12))
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,-1,-0.5,0,0.5,1,1.5])
plt.tight_layout()
##plt.savefig('plots/heatmaps/fc_colorbar.png', dpi=300)
plt.show();



### Fig 1 heatmaps
ticks = [25,125,225,325,425]
labels = ['TSS', '+1 kb', '+2 kb', '+3 kb', '+4 kb']

sns.set_style('ticks', {'axes.facecolor':'#DCDCDC'})
## Rpb1
# ΔNTD IAA vs SPT6 IAA
minimum = np.min([spt6_iaa_Rpb1_averaged.quantile([0.05]).min(1),dntd_iaa_Rpb1_averaged.quantile([0.05]).min(1)])
maximum = np.max([spt6_iaa_Rpb1_averaged.quantile([0.95]).max(1),dntd_iaa_Rpb1_averaged.quantile([0.95]).max(1)])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_iaa_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/reference/Rpb1_heatmap_reference_wt_iaa.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_iaa_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('ΔNTD\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/reference/Rpb1_heatmap_reference_dntd_iaa.png", dpi=600)
plt.show();

#colorbar
a = np.array([[minimum,maximum]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="Greens")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.savefig('plots/heatmaps/Rpb1_heatmaps_colorbar_reference.png', dpi=600)
plt.show();

fc_3 = np.log2(dntd_iaa_Rpb1_averaged/spt6_iaa_Rpb1_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/reference/Rpb1_heatmap_reference_fc_dntd_wt_iaa.png", dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/Rpb1_fc_colorbar_reference.png', dpi=600)
plt.show();

#Rpb1 IAA vs DMSO
fc_3 = np.log2(spt6_iaa_Rpb1_averaged/spt6_dmso_Rpb1_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/reference/Rpb1_heatmap_reference_fc_wt_iaa_dmso.png", dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/Rpb1_fc_colorbar_reference.png', dpi=600)
plt.show();


## V5
# SPT6 IAA vs SPT6 DMSO
minimum = np.min([spt6_dmso_V5_averaged.quantile([0.05]).min(1),
                  spt6_iaa_V5_averaged.quantile([0.05]).min(1),])
maximum = np.max([spt6_dmso_V5_averaged.quantile([0.95]).max(1),
                  spt6_iaa_V5_averaged.quantile([0.95]).max(1),])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('V5 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_dmso_V5_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nDMSO', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/reference/V5_heatmap_reference_wt_dmso.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('V5 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_iaa_V5_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/reference/V5_heatmap_reference_wt_iaa.png", dpi=600)
plt.show();

# #colorbar
# a = np.array([[minimum,maximum]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="Oranges")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# plt.tight_layout()
# plt.savefig('plots/heatmaps/V5_heatmaps_colorbar_reference.png', dpi=600)
# plt.show();

fc_3 = np.log2(spt6_iaa_V5_averaged/spt6_dmso_V5_averaged)
minimum = -5
maximum = 5

f = plt.figure(figsize=(1.5, 2.5), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{SPT6\:IAA}{SPT6\:DMSO}$ V5 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=7)
plt.axvline(25, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/reference/V5_heatmap_reference_fc_wt_iaa_dmso.png", dpi=600)
plt.show();

# #colorbar
# a = np.array([[-5,5]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="bwr")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical', ticks=[-5,-2.5,0,2.5,5])
# plt.tight_layout()
# plt.savefig('plots/heatmaps/V5_fc_colorbar_reference.png', dpi=600)
# plt.show();


## Flag per Rpb1
# ΔNTD IAA vs SPT6 IAA
minimum = np.min([spt6_iaa_FpR_averaged.quantile([0.05]).min(1),dntd_iaa_FpR_averaged.quantile([0.05]).min(1)])
maximum = np.max([spt6_iaa_FpR_averaged.quantile([0.95]).max(1),dntd_iaa_FpR_averaged.quantile([0.95]).max(1)])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_iaa_FpR_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/reference/FpR_heatmap_reference_wt_iaa.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_iaa_FpR_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('ΔNTD\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/reference/FpR_heatmap_reference_dntd_iaa.png", dpi=600)
plt.show();

#colorbar
a = np.array([[minimum,maximum]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="Greens")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.savefig('plots/heatmaps/FpR_heatmaps_colorbar_reference.png', dpi=600)
plt.show();

fc_3 = np.log2(dntd_iaa_FpR_averaged/spt6_iaa_FpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ FpR ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/reference/FpR_heatmap_reference_fc_dntd_wt_iaa.png", dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/FpR_fc_colorbar_reference.png', dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(4,1), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
cbar = plt.colorbar(ticks=[-1.5,0,1.5], orientation='horizontal')
cbar.ax.set_xticklabels(['$\\leq$$-1.5$', '0', '$\\geq$$1.5$'])
plt.tight_layout()
plt.savefig('plots/heatmaps/fc_colorbar_horizontal.png', dpi=600)
plt.show();

#colorbar
a = np.array([[-5,5]])
plt.figure(figsize=(4,1), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
cbar = plt.colorbar(ticks=[-5,0,5], orientation='horizontal')
cbar.ax.set_xticklabels(['$\\leq$$-5$', '0', '$\\geq$$5$'])
plt.tight_layout()
plt.savefig('plots/heatmaps/V5_fc_colorbar_horizontal.png', dpi=600)
plt.show();







