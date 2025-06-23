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

#averaged samples
spt6_pob3_Rpb1_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/661_93_8WG16vinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_e154k_Rpb1_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/666_93_8WG16vinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_pob3_Rpb1_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/661_95_8WG16vinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_e154k_Rpb1_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/666_95_8WG16vinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_pob3_Myc_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/661_93_mycvinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_e154k_Myc_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/666_93_mycvinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_pob3_Myc_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/661_95_mycvinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_e154k_Myc_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/666_95_mycvinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_pob3_Flag_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/661_93_Flagvinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_e154k_Flag_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/666_93_Flagvinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_pob3_Flag_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/661_95_Flagvinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_e154k_Flag_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/ratio/tab/666_95_Flagvinput_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_pob3_MpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/perpol/tab/661_93_mycv8WG16_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_e154k_MpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/perpol/tab/666_93_mycv8WG16_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_pob3_MpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/perpol/tab/661_95_mycv8WG16_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_e154k_MpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/perpol/tab/666_95_mycv8WG16_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_pob3_FpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/perpol/tab/661_93_Flagv8WG16_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

spt6_e154k_FpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/perpol/tab/666_93_Flagv8WG16_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_pob3_FpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/perpol/tab/661_95_Flagv8WG16_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

dntd_e154k_FpR_averaged = pd.read_csv('/Users/jlwarner/Desktop/spt16_chip_seq/transfer/averaged/perpol/tab/666_95_Flagv8WG16_averaged_midpoints_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)


### Fig 6 heatmaps
ticks = [0,25,125,149]
labels = ['-250 bp', 'TSS', 'CPS', '+250 bp']

sns.set_style('ticks', {'axes.facecolor':'white'})

#fold-change colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(0.5,2), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
cbar = plt.colorbar(ticks=[-1.5,0,1.5], orientation='vertical')
cbar.ax.set_yticklabels(['$\\leq$$-1.5$', '0', '$\\geq$$1.5$'], fontsize=8)
plt.tight_layout()
#plt.savefig('plots/heatmaps/scale/fc_colorbar.png', dpi=600)
plt.show();

##Rpb1 log2 fold-change
#ΔNTD POB3 vs SPT6 POB3
fc_3 = np.log2(dntd_pob3_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#SPT6 E154K vs SPT6 POB3
fc_3 = np.log2(spt6_e154k_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_fc_spt6-e154k_spt6-pob3.png", dpi=600)
plt.show();

#ΔNTD E154K vs SPT6 POB3
fc_3 = np.log2(dntd_e154k_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_fc_dntd-e154k_spt6-pob3.png", dpi=600)
plt.show();

#ΔNTD E154K vs ΔNTD POB3
fc_3 = np.log2(dntd_e154k_Rpb1_averaged/dntd_pob3_Rpb1_averaged)
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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();



##Myc log2 fold-change
#ΔNTD POB3 vs SPT6 POB3
fc_3 = np.log2(dntd_pob3_Myc_averaged/spt6_pob3_Myc_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Myc ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#SPT6 E154K vs SPT6 POB3
fc_3 = np.log2(spt6_e154k_Myc_averaged/spt6_pob3_Myc_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Myc ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_fc_spt6-e154k_spt6-pob3.png", dpi=600)
plt.show();

#ΔNTD E154K vs SPT6 POB3
fc_3 = np.log2(dntd_e154k_Myc_averaged/spt6_pob3_Myc_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Myc ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_fc_dntd-e154k_spt6-pob3.png", dpi=600)
plt.show();

#ΔNTD E154K vs ΔNTD POB3
fc_3 = np.log2(dntd_e154k_Myc_averaged/dntd_pob3_Myc_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Myc ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();



##MpR log2 fold-change
#ΔNTD POB3 vs SPT6 POB3
fc_3 = np.log2(dntd_pob3_MpR_averaged/spt6_pob3_MpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ MpR ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#SPT6 E154K vs SPT6 POB3
fc_3 = np.log2(spt6_e154k_MpR_averaged/spt6_pob3_MpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ MpR ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_fc_spt6-e154k_spt6-pob3.png", dpi=600)
plt.show();

#ΔNTD E154K vs SPT6 POB3
fc_3 = np.log2(dntd_e154k_MpR_averaged/spt6_pob3_MpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ MpR ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_fc_dntd-e154k_spt6-pob3.png", dpi=600)
plt.show();

#ΔNTD E154K vs ΔNTD POB3
fc_3 = np.log2(dntd_e154k_MpR_averaged/dntd_pob3_MpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ MpR ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();



##FpR log2 fold-change
#ΔNTD POB3 vs SPT6 POB3
fc_3 = np.log2(dntd_pob3_FpR_averaged/spt6_pob3_FpR_averaged)
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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#SPT6 E154K vs SPT6 POB3
fc_3 = np.log2(spt6_e154k_FpR_averaged/spt6_pob3_FpR_averaged)
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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_fc_spt6-e154k_spt6-pob3.png", dpi=600)
plt.show();

#ΔNTD E154K vs SPT6 POB3
fc_3 = np.log2(dntd_e154k_FpR_averaged/spt6_pob3_FpR_averaged)
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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_fc_dntd-e154k_spt6-pob3.png", dpi=600)
plt.show();

#ΔNTD E154K vs ΔNTD POB3
fc_3 = np.log2(dntd_e154k_FpR_averaged/dntd_pob3_FpR_averaged)
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
plt.axvline(125, color='black', linestyle='--', linewidth=0.5, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();














### Metagenes
ticks = np.array([-250, 0, 1000, 1250])
labels = ['-250 bp', 'TSS', 'CPS', '+250 bp']
#ybottom = -0.1
#ytop = 4.0
#yticks = np.array((0, 1.0, 2.0, 3.0, 4.0))

#Rpb1
plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_Rpb1_averaged.mean(), label = 'SPT6 POB3', color = 'black')
plt.fill_between(positions, spt6_pob3_Rpb1_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_Rpb1_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='black')
plt.plot(positions, dntd_pob3_Rpb1_averaged.mean(), label = 'ΔNTD POB3', color = 'red')
plt.fill_between(positions, dntd_pob3_Rpb1_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_pob3_Rpb1_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='red')
plt.plot(positions, spt6_e154k_Rpb1_averaged.mean(), label = 'SPT6 E154K', color = 'orange')
plt.fill_between(positions, spt6_e154k_Rpb1_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_Rpb1_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='orange')
plt.plot(positions, dntd_e154k_Rpb1_averaged.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.fill_between(positions, dntd_e154k_Rpb1_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_e154k_Rpb1_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized Rpb1 occupancy')
plt.ylim(0, 3.25)
plt.yticks(ticks=np.array((0, 1.0, 2.0, 3.0, 3.25)))
plt.legend(loc='upper left')
plt.title('Rpb1 mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/Rpb1_mean_metagene_midpoints_averaged.png', dpi=300)
plt.show();

fc_0 = np.log2(dntd_pob3_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
fc_1 = np.log2(spt6_e154k_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
fc_2 = np.log2(dntd_e154k_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
fc_3 = np.log2(dntd_e154k_Rpb1_averaged/dntd_pob3_Rpb1_averaged)
plt.figure(dpi=300)
plt.plot(positions, fc_0.mean(), label = 'ΔNTD POB3', color = 'red')
plt.plot(positions, fc_1.mean(), label = 'SPT6 E154K', color = 'orange')
plt.plot(positions, fc_2.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.plot(positions, fc_3.mean(), label = 'ΔNTD E154K / ΔNTD', color = 'black')
plt.axhline(0, linestyle='--', color='black')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
#plt.ylabel('spike-in normalized Rpb1 occupancy')
plt.ylim(-0.35,0.5)
#plt.yticks(ticks=np.array((0, 1.0, 2.0, 3.0, 4.0)))
plt.legend(loc='upper left')
plt.title('log2 fold-change Rpb1 mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/log2fc_Rpb1_metagene_averaged.png', dpi=300)
plt.show();


#Myc
plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_Myc_averaged.mean(), label = 'SPT6 POB3', color = 'black')
plt.fill_between(positions, spt6_pob3_Myc_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_Myc_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='black')
plt.plot(positions, dntd_pob3_Myc_averaged.mean(), label = 'ΔNTD POB3', color = 'red')
plt.fill_between(positions, dntd_pob3_Myc_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_pob3_Myc_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='red')
plt.plot(positions, spt6_e154k_Myc_averaged.mean(), label = 'SPT6 E154K', color = 'orange')
plt.fill_between(positions, spt6_e154k_Myc_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_Myc_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='orange')
plt.plot(positions, dntd_e154k_Myc_averaged.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.fill_between(positions, dntd_e154k_Myc_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_e154k_Myc_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized Myc occupancy')
plt.ylim(0, 4.0)
plt.yticks(ticks=np.array((0, 1.0, 2.0, 3.0, 4.0)))
plt.legend(loc='upper left')
plt.title('Myc mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/Myc_mean_metagene_midpoints_averaged.png', dpi=300)
plt.show();


fc_0 = np.log2(dntd_pob3_Myc_averaged/spt6_pob3_Myc_averaged)
fc_1 = np.log2(spt6_e154k_Myc_averaged/spt6_pob3_Myc_averaged)
fc_2 = np.log2(dntd_e154k_Myc_averaged/spt6_pob3_Myc_averaged)
fc_3 = np.log2(dntd_e154k_Myc_averaged/dntd_pob3_Myc_averaged)
plt.figure(dpi=300)
plt.plot(positions, fc_0.mean(), label = 'ΔNTD POB3', color = 'red')
plt.plot(positions, fc_1.mean(), label = 'SPT6 E154K', color = 'orange')
plt.plot(positions, fc_2.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.plot(positions, fc_3.mean(), label = 'ΔNTD E154K / ΔNTD', color = 'black')
plt.axhline(0, linestyle='--', color='black')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
#plt.ylabel('spike-in normalized Myc occupancy')
plt.ylim(-0.35,0.5)
#plt.yticks(ticks=np.array((0, 1.0, 2.0, 3.0, 4.0)))
plt.legend(loc='upper left')
plt.title('log2 fold-change Myc mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/log2fc_Myc_metagene_averaged.png', dpi=300)
plt.show();



#Flag
plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_Flag_averaged.mean(), label = 'SPT6 POB3', color = 'black')
plt.fill_between(positions, spt6_pob3_Flag_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_Flag_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='black')
plt.plot(positions, dntd_pob3_Flag_averaged.mean(), label = 'ΔNTD POB3', color = 'red')
plt.fill_between(positions, dntd_pob3_Flag_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_pob3_Flag_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='red')
plt.plot(positions, spt6_e154k_Flag_averaged.mean(), label = 'SPT6 E154K', color = 'orange')
plt.fill_between(positions, spt6_e154k_Flag_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_Flag_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='orange')
plt.plot(positions, dntd_e154k_Flag_averaged.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.fill_between(positions, dntd_e154k_Flag_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_e154k_Flag_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized Flag occupancy')
plt.ylim(0, 2.5)
plt.yticks(ticks=np.array((0, 1.0, 2.0, 2.5)))
plt.legend(loc='upper left')
plt.title('Flag mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/Flag_mean_metagene_midpoints_averaged.png', dpi=300)
plt.show();

fc_0 = np.log2(dntd_pob3_Flag_averaged/spt6_pob3_Flag_averaged)
fc_1 = np.log2(spt6_e154k_Flag_averaged/spt6_pob3_Flag_averaged)
fc_2 = np.log2(dntd_e154k_Flag_averaged/spt6_pob3_Flag_averaged)
fc_3 = np.log2(dntd_e154k_Flag_averaged/dntd_pob3_Flag_averaged)
plt.figure(dpi=300)
plt.plot(positions, fc_0.mean(), label = 'ΔNTD POB3', color = 'red')
plt.plot(positions, fc_1.mean(), label = 'SPT6 E154K', color = 'orange')
plt.plot(positions, fc_2.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.plot(positions, fc_3.mean(), label = 'ΔNTD E154K / ΔNTD', color = 'black')
plt.axhline(0, linestyle='--', color='black')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
#plt.ylabel('spike-in normalized Flag occupancy')
plt.ylim(-0.35,0.5)
#plt.yticks(ticks=np.array((0, 1.0, 2.0, 3.0, 4.0)))
plt.legend(loc='upper left')
plt.title('log2 fold-change Flag mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/log2fc_Flag_metagene_averaged.png', dpi=300)
plt.show();

#Myc / Rpb1
plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_MpR_averaged.mean(), label = 'SPT6 POB3', color = 'black')
plt.fill_between(positions, spt6_pob3_MpR_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_MpR_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='black')
plt.plot(positions, dntd_pob3_MpR_averaged.mean(), label = 'ΔNTD POB3', color = 'red')
plt.fill_between(positions, dntd_pob3_MpR_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_pob3_MpR_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='red')
plt.plot(positions, spt6_e154k_MpR_averaged.mean(), label = 'SPT6 E154K', color = 'orange')
plt.fill_between(positions, spt6_e154k_MpR_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_MpR_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='orange')
plt.plot(positions, dntd_e154k_MpR_averaged.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.fill_between(positions, dntd_e154k_MpR_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_e154k_MpR_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized $\\frac{Myc}{Rpb1}$ occupancy')
plt.ylim(0, 2)
plt.yticks(ticks=np.array((0, 1.0, 2.0)))
plt.legend(loc='upper left')
plt.title('MpR mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/MpR_mean_metagene_midpoints_averaged.png', dpi=300)
plt.show();

fc_0 = np.log2(dntd_pob3_MpR_averaged/spt6_pob3_MpR_averaged)
fc_1 = np.log2(spt6_e154k_MpR_averaged/spt6_pob3_MpR_averaged)
fc_2 = np.log2(dntd_e154k_MpR_averaged/spt6_pob3_MpR_averaged)
fc_3 = np.log2(dntd_e154k_MpR_averaged/dntd_pob3_MpR_averaged)
plt.figure(dpi=300)
plt.plot(positions, fc_0.mean(), label = 'ΔNTD POB3', color = 'red')
plt.plot(positions, fc_1.mean(), label = 'SPT6 E154K', color = 'orange')
plt.plot(positions, fc_2.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.plot(positions, fc_3.mean(), label = 'ΔNTD E154K / ΔNTD', color = 'black')
plt.axhline(0, linestyle='--', color='black')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
#plt.ylabel('spike-in normalized MpR occupancy')
plt.ylim(-0.35,0.5)
#plt.yticks(ticks=np.array((0, 1.0, 2.0, 3.0, 4.0)))
plt.legend(loc='upper left')
plt.title('log2 fold-change MpR mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/log2fc_MpR_metagene_averaged.png', dpi=300)
plt.show();

#Flag / Rpb1
plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_FpR_averaged.mean(), label = 'SPT6 POB3', color = 'black')
plt.fill_between(positions, spt6_pob3_FpR_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_pob3_FpR_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='black')
plt.plot(positions, dntd_pob3_FpR_averaged.mean(), label = 'ΔNTD POB3', color = 'red')
plt.fill_between(positions, dntd_pob3_FpR_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_pob3_FpR_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='red')
plt.plot(positions, spt6_e154k_FpR_averaged.mean(), label = 'SPT6 E154K', color = 'orange')
plt.fill_between(positions, spt6_e154k_FpR_averaged.quantile([0.25,0.75]).iloc[0,:], spt6_e154k_FpR_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='orange')
plt.plot(positions, dntd_e154k_FpR_averaged.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.fill_between(positions, dntd_e154k_FpR_averaged.quantile([0.25,0.75]).iloc[0,:], dntd_e154k_FpR_averaged.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.1, color='blue')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized $\\frac{Flag}{Rpb1}$ occupancy')
plt.ylim(0, 2)
plt.yticks(ticks=np.array((0, 1.0, 2.0)))
plt.legend(loc='upper left')
plt.title('FpR mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/FpR_mean_metagene_midpoints_averaged.png', dpi=300)
plt.show();

fc_0 = np.log2(dntd_pob3_FpR_averaged/spt6_pob3_FpR_averaged)
fc_1 = np.log2(spt6_e154k_FpR_averaged/spt6_pob3_FpR_averaged)
fc_2 = np.log2(dntd_e154k_FpR_averaged/spt6_pob3_FpR_averaged)
fc_3 = np.log2(dntd_e154k_FpR_averaged/dntd_pob3_FpR_averaged)
plt.figure(dpi=300)
plt.plot(positions, fc_0.mean(), label = 'ΔNTD POB3', color = 'red')
plt.plot(positions, fc_1.mean(), label = 'SPT6 E154K', color = 'orange')
plt.plot(positions, fc_2.mean(), label = 'ΔNTD E154K', color = 'blue')
plt.plot(positions, fc_3.mean(), label = 'ΔNTD E154K / ΔNTD', color = 'black')
plt.axhline(0, linestyle='--', color='black')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
#plt.ylabel('spike-in normalized FpR occupancy')
plt.ylim(-0.35,0.5)
#plt.yticks(ticks=np.array((0, 1.0, 2.0, 3.0, 4.0)))
plt.legend(loc='upper left')
plt.title('log2 fold-change FpR mean midpoints_averaged')
plt.tight_layout()
#plt.savefig('plots/metagenes/log2fc_FpR_metagene_averaged.png', dpi=300)
plt.show();





## Violin Plots
# Rpb1
wt = np.log2(spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
dntd = np.log2((dntd_pob3_Rpb1_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_Rpb1_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_Rpb1_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'WT':wt, 'spt6ΔNTD':dntd, 'pob3-E154K':e154k, 'spt6ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 Rpb1', fontsize=14)
plt.title('Rpb1 midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/Rpb1_violin_midpoints_averaged.png', dpi=300)
plt.show();

# Myc
wt = np.log2(spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1))
dntd = np.log2((dntd_pob3_Myc_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_Myc_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_Myc_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'WT':wt, 'spt6ΔNTD':dntd, 'pob3-E154K':e154k, 'spt6ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 Myc', fontsize=14)
plt.title('Myc midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/Myc_violin_midpoints_averaged.png', dpi=300)
plt.show();

#Flag
wt = np.log2(spt6_pob3_Flag_averaged.iloc[:,25:125].mean(1))
dntd = np.log2((dntd_pob3_Flag_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_Flag_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_Flag_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'WT':wt, 'spt6ΔNTD':dntd, 'pob3-E154K':e154k, 'spt6ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 Flag', fontsize=14)
plt.title('Flag midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/Flag_violin_midpoints_averaged.png', dpi=300)
plt.show();

# Myc / Rpb1
wt = np.log2(spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1))
dntd = np.log2((dntd_pob3_MpR_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_MpR_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_MpR_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'WT':wt, 'spt6ΔNTD':dntd, 'pob3-E154K':e154k, 'spt6ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2$\\frac{Myc}{Rpb1}$', fontsize=14)
plt.title('MpR midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/MpR_violin_midpoints_averaged.png', dpi=300)
plt.show();

# Flag / Rpb1
wt = np.log2(spt6_pob3_FpR_averaged.iloc[:,25:125].mean(1))
dntd = np.log2((dntd_pob3_FpR_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_FpR_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_FpR_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'WT':wt, 'spt6ΔNTD':dntd, 'pob3-E154K':e154k, 'spt6ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2$\\frac{Flag}{Rpb1}$', fontsize=14)
plt.title('FpR midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/FpR_violin_midpoints_averaged.png', dpi=300)
plt.show();


## Scatterplots

#Rpb1 Scatterplots
x = np.log2(spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
y = np.log2(dntd_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-5.25,5.25)
plt.xticks(np.arange(-5,5.5,2.5), fontsize=14)
plt.ylim(-5.25,5.25)
plt.yticks(np.arange(-5,5.5,2.5), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD POB3', fontsize=18, style='italic')
plt.title('midpoints_averaged\nlog2(Rpb1 enrichment)', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntdvswt_Rpb1_midpoints_averaged_scatter.png')
plt.show();

x = np.log2(spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
y = np.log2(spt6_e154k_Rpb1_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-5.25,5.25)
plt.xticks(np.arange(-5,5.5,2.5), fontsize=14)
plt.ylim(-5.25,5.25)
plt.yticks(np.arange(-5,5.5,2.5), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('SPT6 pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\nlog2(Rpb1 enrichment)', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/e154kvswt_Rpb1_midpoints_averaged_scatter.png')
plt.show();

x = np.log2(spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
y = np.log2(dntd_e154k_Rpb1_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-5.25,5.25)
plt.xticks(np.arange(-5,5.5,2.5), fontsize=14)
plt.ylim(-5.25,5.25)
plt.yticks(np.arange(-5,5.5,2.5), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\nlog2(Rpb1 enrichment)', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntde154kvswt_Rpb1_midpoints_averaged_scatter.png')
plt.show();

x = np.log2(dntd_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
y = np.log2(dntd_e154k_Rpb1_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-5.25,5.25)
plt.xticks(np.arange(-5,5.5,2.5), fontsize=14)
plt.ylim(-5.25,5.25)
plt.yticks(np.arange(-5,5.5,2.5), fontsize=14)
plt.xlabel('spt6ΔNTD POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\nlog2(Rpb1 enrichment)', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntde154kvsdntd_Rpb1_midpoints_averaged_scatter.png')
plt.show();

#Myc Scatterplots
x = np.log2(spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1))
y = np.log2(dntd_pob3_Myc_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-2,5.25)
plt.xticks(np.arange(-2,5.5,2), fontsize=14)
plt.ylim(-2,5.25)
plt.yticks(np.arange(-2,5.5,2), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD POB3', fontsize=18, style='italic')
plt.title('midpoints_averaged\nlog2(Myc enrichment)', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntdvswt_Myc_midpoints_averaged_scatter.png')
plt.show();

x = np.log2(spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1))
y = np.log2(spt6_e154k_Myc_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-2,5.25)
plt.xticks(np.arange(-2,5.5,2), fontsize=14)
plt.ylim(-2,5.25)
plt.yticks(np.arange(-2,5.5,2), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('SPT6 pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\nlog2(Myc enrichment)', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/e154kvswt_Myc_midpoints_averaged_scatter.png')
plt.show();

x = np.log2(spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1))
y = np.log2(dntd_e154k_Myc_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-2,5.25)
plt.xticks(np.arange(-2,5.5,2), fontsize=14)
plt.ylim(-2,5.25)
plt.yticks(np.arange(-2,5.5,2), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\nlog2(Myc enrichment)', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntde154kvswt_Myc_midpoints_averaged_scatter.png')
plt.show();

x = np.log2(dntd_pob3_Myc_averaged.iloc[:,25:125].mean(1))
y = np.log2(dntd_e154k_Myc_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-2,5.25)
plt.xticks(np.arange(-2,5.5,2), fontsize=14)
plt.ylim(-2,5.25)
plt.yticks(np.arange(-2,5.5,2), fontsize=14)
plt.xlabel('spt6ΔNTD POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\nlog2(Myc enrichment)', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntde154kvsdntd_Myc_midpoints_averaged_scatter.png')
plt.show();


#MpR Scatterplots
x = (spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1))
y = (dntd_pob3_MpR_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=8, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(0.75,2)
plt.xticks(np.arange(0.75,2.1,0.25), fontsize=14)
plt.ylim(0.75,2)
plt.yticks(np.arange(0.75,2.1,0.25), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD POB3', fontsize=18, style='italic')
plt.title('midpoints_averaged\n$\\frac{Myc}{Rpb1}$ enrichment', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntdvswt_MpR_midpoints_averaged_scatter.png')
plt.show();

x = (spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1))
y = (spt6_e154k_MpR_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=8, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(0.75,2)
plt.xticks(np.arange(0.75,2.1,0.25), fontsize=14)
plt.ylim(0.75,2)
plt.yticks(np.arange(0.75,2.1,0.25), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('SPT6 pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\n$\\frac{Myc}{Rpb1}$ enrichment', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/e154kvswt_MpR_midpoints_averaged_scatter.png')
plt.show();

x = (spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1))
y = (dntd_e154k_MpR_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=8, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(0.75,2)
plt.xticks(np.arange(0.75,2.1,0.25), fontsize=14)
plt.ylim(0.75,2)
plt.yticks(np.arange(0.75,2.1,0.25), fontsize=14)
plt.xlabel('SPT6 POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\n$\\frac{Myc}{Rpb1}$ enrichment', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntde154kvswt_MpR_midpoints_averaged_scatter.png')
plt.show();

x = (dntd_pob3_MpR_averaged.iloc[:,25:125].mean(1))
y = (dntd_e154k_MpR_averaged.iloc[:,25:125].mean(1))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=8, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(0.75,2)
plt.xticks(np.arange(0.75,2.1,0.25), fontsize=14)
plt.ylim(0.75,2)
plt.yticks(np.arange(0.75,2.1,0.25), fontsize=14)
plt.xlabel('spt6ΔNTD POB3', fontsize=18, style='italic')
plt.ylabel('spt6ΔNTD pob3-E154K', fontsize=18, style='italic')
plt.title('midpoints_averaged\n$\\frac{Myc}{Rpb1}$ enrichment', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/dntde154kvsdntd_MpR_midpoints_averaged_scatter.png')
plt.show();


#Old scatters
x = spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1)
y = dntd_pob3_Myc_averaged.iloc[:,25:125].mean(1)
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,35.5],[-0.5,35.5], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-0.5,35.5)
plt.xticks(np.arange(0,35.5,10), fontsize=16)
plt.xlabel('Spt6', fontsize=18)
plt.ylabel('spt6ΔNTD', fontsize=18)
plt.yticks(np.arange(0,35.5,10), fontsize=16)
plt.ylim(-0.5,35.5)
plt.title('Myc midpoints_averaged', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/Myc_midpoints_averaged_scatter.png')
plt.show();


x = spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1)
y = dntd_pob3_MpR_averaged.iloc[:,25:125].mean(1)
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,4.5],[-0.5,4.5], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-0.5,4.5)
plt.xticks(np.arange(0,4.5,1), fontsize=16)
plt.xlabel('$SPT6 \: POB3$', fontsize=18)
plt.ylabel('$ΔNTD \: POB3$', fontsize=18)
plt.yticks(np.arange(0,4.5,1), fontsize=16)
plt.ylim(-0.5,4.5)
plt.title('$\\frac{Myc}{Rpb1}$ midpoints_averaged', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/MpR_midpoints_averaged_scatter.png')
plt.show();

x = spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1)
y = spt6_e154k_MpR_averaged.iloc[:,25:125].mean(1)
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,4.5],[-0.5,4.5], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-0.5,4.5)
plt.xticks(np.arange(0,4.5,1), fontsize=16)
plt.xlabel('$SPT6 \: POB3$', fontsize=18)
plt.ylabel('$SPT6 \: E154K$', fontsize=18)
plt.yticks(np.arange(0,4.5,1), fontsize=16)
plt.ylim(-0.5,4.5)
plt.title('$\\frac{Myc}{Rpb1}$ midpoints_averaged', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/MpR_midpoints_averaged_scatter.png')
plt.show();

x = spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1)
y = dntd_e154k_MpR_averaged.iloc[:,25:125].mean(1)
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,4.5],[-0.5,4.5], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(0.5,2)
#plt.xticks(np.arange(0,4.5,1), fontsize=16)
plt.xlabel('$SPT6 \: POB3$', fontsize=18)
plt.ylabel('$ΔNTD \: E154K$', fontsize=18)
#plt.yticks(np.arange(0,4.5,1), fontsize=16)
plt.ylim(0.5,2)
plt.title('$\\frac{Myc}{Rpb1}$ midpoints_averaged', fontsize=20, pad=10)
plt.tight_layout()
#plt.savefig('plots/scatter/MpR_midpoints_averaged_scatter.png')
plt.show();


#log2fc scatters
x = np.log2(spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
y = np.log2((dntd_pob3_Myc_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1)))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-1.5,6.5],[0,0], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-1.5,6)
plt.xticks(np.arange(-1.5,6.5,1.5), fontsize=12)
plt.xlabel('SPT6 POB3:\n Rpb1 enrichment', fontsize=14)
plt.ylabel('log2 $\\frac{ΔNTD \: POB3}{wild-type}$\nMyc enrichment', fontsize=14)
plt.text(4.9, 1.75, str((y>0).sum()))
plt.text(4.9, -1.9, str((y<0).sum()))
#plt.yticks(np.arange(0,4.5,1), fontsize=16)
plt.ylim(-2,2)
plt.title('midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/scatter/log2fc_dntd-wt_scatter_midpoints_averaged.png', dpi=300)
plt.show();


x = np.log2(spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
y = np.log2((spt6_e154k_Myc_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1)))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-1.5,6.5],[0,0], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-1.5,6)
plt.xticks(np.arange(-1.5,6.5,1.5), fontsize=12)
plt.xlabel('SPT6 POB3:\n Rpb1 enrichment', fontsize=14)
plt.ylabel('log2 $\\frac{SPT6 \: E154K}{wild-type}$\nMyc enrichment', fontsize=14)
plt.text(4.9, 1.75, str((y>0).sum()))
plt.text(4.9, -1.9, str((y<0).sum()))
#plt.yticks(np.arange(0,4.5,1), fontsize=16)
plt.ylim(-2,2)
plt.title('midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/scatter/log2fc_e154k-wt_scatter_midpoints_averaged.png', dpi=300)
plt.show();

x = np.log2(spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1))
y = np.log2((dntd_e154k_Myc_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1)))
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-1.5,6.5],[0,0], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-1.5,6)
plt.xticks(np.arange(-1.5,6.5,1.5), fontsize=12)
plt.xlabel('SPT6 POB3:\n Rpb1 enrichment', fontsize=14)
plt.ylabel('log2 $\\frac{ΔNTD \: E154K}{wild-type}$\nMyc enrichment', fontsize=14)
plt.text(4.9, 1.75, str((y>0).sum()))
plt.text(4.9, -1.9, str((y<0).sum()))
#plt.yticks(np.arange(0,4.5,1), fontsize=16)
plt.ylim(-2,2)
plt.title('midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/scatter/log2fc_dntd_e154k-wt_scatter_midpoints_averaged.png', dpi=300)
plt.show();



# #NEW Scatterplots
# x = np.log2(spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1))
# y = np.log2(spt6_pob3_FpR_averaged.iloc[:,25:125].mean(1))
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-1,2)
# #plt.xticks(np.arange(-5,5.5,2.5), fontsize=14)
# plt.ylim(-1,2)
# #plt.yticks(np.arange(-5,5.5,2.5), fontsize=14)
# plt.xlabel('FACT/Rpb1', fontsize=18, style='italic')
# plt.ylabel('Spt6/Rpb1', fontsize=18, style='italic')
# plt.title('WT', fontsize=20, pad=10)
# plt.tight_layout()
# #plt.savefig('plots/scatter/dntdvswt_Rpb1_midpoints_averaged_scatter.png')
# plt.show();

# x = np.log2(dntd_pob3_MpR_averaged.iloc[:,25:125].mean(1))
# y = np.log2(dntd_pob3_FpR_averaged.iloc[:,25:125].mean(1))
# df = pd.DataFrame({'x':x, 'y':y})

# plt.figure(figsize=(4,4), dpi=300)
# values = np.vstack([x, y])
# kernel = scipy.stats.gaussian_kde(values)(values)
# plt.plot([-5.5,8],[-5.5,8], color='black', linestyle='--', linewidth=1, zorder=1)
# sns.scatterplot(data=df, x='x', y='y', 
#                 s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
# plt.xlim(-1,2)
# #plt.xticks(np.arange(-5,5.5,2.5), fontsize=14)
# plt.ylim(-1,2)
# #plt.yticks(np.arange(-5,5.5,2.5), fontsize=14)
# plt.xlabel('FACT/Rpb1', fontsize=18, style='italic')
# plt.ylabel('Spt6/Rpb1', fontsize=18, style='italic')
# plt.title('DeltaNTD', fontsize=20, pad=10)
# plt.tight_layout()
# #plt.savefig('plots/scatter/dntdvswt_Rpb1_midpoints_averaged_scatter.png')
# plt.show();




# Violin plots of Rpb1 log2 fc
dntd = np.log2((dntd_pob3_Rpb1_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_Rpb1_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_Rpb1_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Rpb1_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'ΔNTD':dntd, 'pob3-E154K':e154k, 'ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 $\\frac{mutant}{wild-type}$ Rpb1', fontsize=14)
plt.title('Rpb1 midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/log2fc_Rpb1_violin_midpoints_averaged.png', dpi=300)
plt.show();

# Violin plots of myc log2 fc
dntd = np.log2((dntd_pob3_Myc_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_Myc_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_Myc_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Myc_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'ΔNTD':dntd, 'pob3-E154K':e154k, 'ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 $\\frac{mutant}{wild-type}$ Myc', fontsize=14)
plt.title('Myc midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/log2fc_Myc_violin_midpoints_averaged.png', dpi=300)
plt.show();

# Violin plots of Flag log2 fc
dntd = np.log2((dntd_pob3_Flag_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Flag_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_Flag_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Flag_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_Flag_averaged.iloc[:,25:125].mean(1)/spt6_pob3_Flag_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'ΔNTD':dntd, 'pob3-E154K':e154k, 'ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 $\\frac{mutant}{wild-type}$ Flag', fontsize=14)
plt.title('Flag midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/log2fc_Flag_violin_midpoints_averaged.png', dpi=300)
plt.show();

# Violin plots of myc/Rpb1 log2 fc
dntd = np.log2((dntd_pob3_MpR_averaged.iloc[:,25:125].mean(1)/spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_MpR_averaged.iloc[:,25:125].mean(1)/spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_MpR_averaged.iloc[:,25:125].mean(1)/spt6_pob3_MpR_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'ΔNTD':dntd, 'pob3-E154K':e154k, 'ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 $\\frac{mutant}{wild-type} \: \\frac{Myc}{Rpb1}$ ', fontsize=14)
plt.title('MpR midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/log2fc_MpR_violin_midpoints_averaged.png', dpi=300)
plt.show();

# Violin plots of Flag/Rpb1 log2 fc
dntd = np.log2((dntd_pob3_FpR_averaged.iloc[:,25:125].mean(1)/spt6_pob3_FpR_averaged.iloc[:,25:125].mean(1)))
e154k = np.log2((spt6_e154k_FpR_averaged.iloc[:,25:125].mean(1)/spt6_pob3_FpR_averaged.iloc[:,25:125].mean(1)))
dntd_e154k = np.log2((dntd_e154k_FpR_averaged.iloc[:,25:125].mean(1)/spt6_pob3_FpR_averaged.iloc[:,25:125].mean(1)))

df = pd.DataFrame({'ΔNTD':dntd, 'pob3-E154K':e154k, 'ΔNTD\npob3-E154K':dntd_e154k})
plt.figure(figsize=(4,3), dpi=300)
sns.violinplot(data=df,
               linecolor='black')
plt.axhline(0, color = 'black', linestyle='--', linewidth = 1, zorder=1)
plt.ylabel('log2 $\\frac{mutant}{wild-type} \: \\frac{Flag}{Rpb1}$ ', fontsize=14)
plt.title('FpR midpoints_averaged', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('plots/violin/log2fc_FpR_violin_midpoints_averaged.png', dpi=300)
plt.show();





#Heatmaps
ticks = [0,25,125,149]
labels = ['-250 bp', 'TSS', 'CPS', '+250 bp']

sns.set_style('ticks', {'axes.facecolor':'white'})

# Rpb1 Heatmaps
minimum = np.min([spt6_pob3_Rpb1_averaged.quantile([0.10]).min(1),dntd_pob3_Rpb1_averaged.quantile([0.10]).min(1),
                  spt6_e154k_Rpb1_averaged.quantile([0.10]).min(1),dntd_e154k_Rpb1_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_pob3_Rpb1_averaged.quantile([0.90]).max(1),dntd_pob3_Rpb1_averaged.quantile([0.90]).max(1),
                  spt6_e154k_Rpb1_averaged.quantile([0.90]).max(1),dntd_e154k_Rpb1_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('Rpb1 ChIP occupancy midpoints_averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_e154k_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_e154k_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/Rpb1_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

#colorbar
# a = np.array([[0,4]])
# plt.figure(figsize=(1,12))
# img = plt.imshow(a, cmap="Greens")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# #plt.savefig('Rpb1_colorbar.svg', dpi=300)
# plt.show();


#fold-change
fc_0 = np.log2(dntd_pob3_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
fc_1 = np.log2(spt6_e154k_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
fc_2 = np.log2(dntd_e154k_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ Rpb1 ChIP occupancy midpoints_averaged', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/Rpb1_log2fc_heatmap_midpoints_averaged.png", dpi=300)
plt.show()

fc_3 = np.log2(dntd_e154k_Rpb1_averaged/dntd_pob3_Rpb1_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('midpoints_averaged\nlog2 $\\frac{spt6ΔNTD\:pob3-E154K}{spt6ΔNTD\:POB3}$ Rpb1 ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/Rpb1_log2fc_dntd-e154kvsdntd_heatmap_midpoints_averaged.png", dpi=300)
plt.show();


# Myc Heatmaps
minimum = np.min([spt6_pob3_Myc_averaged.quantile([0.10]).min(1),dntd_pob3_Myc_averaged.quantile([0.10]).min(1),
                  spt6_e154k_Myc_averaged.quantile([0.10]).min(1),dntd_e154k_Myc_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_pob3_Myc_averaged.quantile([0.90]).max(1),dntd_pob3_Myc_averaged.quantile([0.90]).max(1),
                  spt6_e154k_Myc_averaged.quantile([0.90]).max(1),dntd_e154k_Myc_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('Myc ChIP occupancy midpoints_averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_Myc_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_Myc_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_e154k_Myc_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_e154k_Myc_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/Myc_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

#fold-change
fc_0 = np.log2(dntd_pob3_Myc_averaged/spt6_pob3_Myc_averaged)
fc_1 = np.log2(spt6_e154k_Myc_averaged/spt6_pob3_Myc_averaged)
fc_2 = np.log2(dntd_e154k_Myc_averaged/spt6_pob3_Myc_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ Myc ChIP occupancy midpoints_averaged', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

# ax = f.add_subplot(gs[0,3])
# a = np.array([[minimum, maximum]])
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')

f.tight_layout()
plt.savefig("plots/heatmaps/Myc_log2fc_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

# a = np.array([[minimum, maximum]])
# plt.figure(figsize=(2,16), dpi=300)
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# plt.show();

fc_3 = np.log2(dntd_e154k_Myc_averaged/dntd_pob3_Myc_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('midpoints_averaged\nlog2 $\\frac{spt6ΔNTD\:pob3-E154K}{spt6ΔNTD\:POB3}$ Myc ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/Myc_log2fc_dntd-e154kvsdntd_heatmap_midpoints_averaged.png", dpi=300)
plt.show();


# Flag Heatmaps
minimum = np.min([spt6_pob3_Flag_averaged.quantile([0.10]).min(1),dntd_pob3_Flag_averaged.quantile([0.10]).min(1),
                  spt6_e154k_Flag_averaged.quantile([0.10]).min(1),dntd_e154k_Flag_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_pob3_Flag_averaged.quantile([0.90]).max(1),dntd_pob3_Flag_averaged.quantile([0.90]).max(1),
                  spt6_e154k_Flag_averaged.quantile([0.90]).max(1),dntd_e154k_Flag_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('Flag ChIP occupancy midpoints_averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_Flag_averaged, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_Flag_averaged, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_e154k_Flag_averaged, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_e154k_Flag_averaged, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/Flag_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

#fold-change
fc_0 = np.log2(dntd_pob3_Flag_averaged/spt6_pob3_Flag_averaged)
fc_1 = np.log2(spt6_e154k_Flag_averaged/spt6_pob3_Flag_averaged)
fc_2 = np.log2(dntd_e154k_Flag_averaged/spt6_pob3_Flag_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ Flag ChIP occupancy midpoints_averaged', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22, style='italic')    

# ax = f.add_subplot(gs[0,3])
# a = np.array([[minimum, maximum]])
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')

f.tight_layout()
#plt.savefig("plots/heatmaps/Flag_log2fc_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

fc_3 = np.log2(dntd_e154k_Flag_averaged/dntd_pob3_Flag_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('midpoints_averaged\nlog2 $\\frac{spt6ΔNTD\:pob3-E154K}{spt6ΔNTD\:POB3}$ Flag ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
#plt.savefig("plots/heatmaps/Flag_log2fc_dntde154kvsdntd_heatmap_midpoints_averaged.png", dpi=300)
plt.show();



# MpR Heatmaps
minimum = np.min([spt6_pob3_MpR_averaged.quantile([0.10]).min(1),dntd_pob3_MpR_averaged.quantile([0.10]).min(1),
                  spt6_e154k_MpR_averaged.quantile([0.10]).min(1),dntd_e154k_MpR_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_pob3_MpR_averaged.quantile([0.90]).max(1),dntd_pob3_MpR_averaged.quantile([0.90]).max(1),
                  spt6_e154k_MpR_averaged.quantile([0.90]).max(1),dntd_e154k_MpR_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('$\\frac{Myc}{Rpb1}$ ChIP occupancy midpoints_averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_MpR_averaged, cmap=sns.color_palette("Reds", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_MpR_averaged, cmap=sns.color_palette("Reds", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_e154k_MpR_averaged, cmap=sns.color_palette("Reds", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_e154k_MpR_averaged, cmap=sns.color_palette("Reds", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/MpR_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

#fold-change
fc_0 = np.log2(dntd_pob3_MpR_averaged/spt6_pob3_MpR_averaged)
fc_1 = np.log2(spt6_e154k_MpR_averaged/spt6_pob3_MpR_averaged)
fc_2 = np.log2(dntd_e154k_MpR_averaged/spt6_pob3_MpR_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ MpR ChIP occupancy midpoints_averaged', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

# ax = f.add_subplot(gs[0,3])
# a = np.array([[minimum, maximum]])
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')

f.tight_layout()
plt.savefig("plots/heatmaps/MpR_log2fc_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

fc_3 = np.log2(dntd_e154k_MpR_averaged/dntd_pob3_MpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('midpoints_averaged\nlog2 $\\frac{spt6ΔNTD\:pob3-E154K}{spt6ΔNTD\:POB3}$ MpR ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/MpR_log2fc_dntd-e154kvsdntd_heatmap_midpoints_averaged.png", dpi=300)
plt.show();


# FpR Heatmaps
minimum = np.min([spt6_pob3_FpR_averaged.quantile([0.10]).min(1),dntd_pob3_FpR_averaged.quantile([0.10]).min(1),
                  spt6_e154k_FpR_averaged.quantile([0.10]).min(1),dntd_e154k_FpR_averaged.quantile([0.10]).min(1)])
maximum = np.max([spt6_pob3_FpR_averaged.quantile([0.90]).max(1),dntd_pob3_FpR_averaged.quantile([0.90]).max(1),
                  spt6_e154k_FpR_averaged.quantile([0.90]).max(1),dntd_e154k_FpR_averaged.quantile([0.90]).max(1)])
f = plt.figure(figsize=(16, 16), dpi=300)
gs = f.add_gridspec(1,4)
f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy midpoints_averaged', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_FpR_averaged, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\nPOB3', fontsize=22)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_FpR_averaged, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\nPOB3', fontsize=22)  

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=spt6_e154k_FpR_averaged, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('SPT6\npob3-E154K', fontsize=22)    

ax = f.add_subplot(gs[0,3])
sns.heatmap(data=dntd_e154k_FpR_averaged, cmap=sns.color_palette("Purples", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/FpR_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

#fold-change
fc_0 = np.log2(dntd_pob3_FpR_averaged/spt6_pob3_FpR_averaged)
fc_1 = np.log2(spt6_e154k_FpR_averaged/spt6_pob3_FpR_averaged)
fc_2 = np.log2(dntd_e154k_FpR_averaged/spt6_pob3_FpR_averaged)
# Set min and max artifically so that the colorbar is symmetrical and 0 = white
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(12, 16), dpi=300)
gs = f.add_gridspec(1,3)
f.suptitle('log2 $\\frac{mutant}{wild-type}$ FpR ChIP occupancy midpoints_averaged', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
plt.title('spt6ΔNTD\npob3-E154K', fontsize=22)    

# ax = f.add_subplot(gs[0,3])
# a = np.array([[minimum, maximum]])
# img = plt.imshow(a, cmap='bwr')
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')

f.tight_layout()
plt.savefig("plots/heatmaps/FpR_log2fc_heatmap_midpoints_averaged.png", dpi=300)
plt.show();

fc_3 = np.log2(dntd_e154k_FpR_averaged/dntd_pob3_FpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(4, 16), dpi=300)
gs = f.add_gridspec(1,1)
f.suptitle('midpoints_averaged\nlog2 $\\frac{spt6ΔNTD\:pob3-E154K}{spt6ΔNTD\:POB3}$ FpR ChIP occupancy', fontsize = 30)

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
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
#plt.savefig("plots/heatmaps/FpR_log2fc_dntde154kvsdntd_heatmap_midpoints_averaged.png", dpi=300)
plt.show();


#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,12))
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,-1,-0.5,0,0.5,1,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/fc_colorbar.png', dpi=300)
plt.show();



### Fig 6 heatmaps
ticks = [0,25,125,149]
labels = ['-250 bp', 'TSS', 'CPS', '+250 bp']

sns.set_style('ticks', {'axes.facecolor':'white'})
## Rpb1
# ΔNTD POB3 vs SPT6 POB3
minimum = np.min([spt6_pob3_Rpb1_averaged.quantile([0.05]).min(1),dntd_pob3_Rpb1_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([spt6_pob3_Rpb1_averaged.quantile([0.95]).max(1),dntd_pob3_Rpb1_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_spt6_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[minimum,maximum]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="Greens")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/Rpb1_heatmaps_colorbar_scale.png', dpi=600)
plt.show();

fc_3 = np.log2(dntd_pob3_Rpb1_averaged/spt6_pob3_Rpb1_averaged)
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
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/Rpb1_fc_colorbar_scale.png', dpi=600)
plt.show();




## Flag
# ΔNTD POB3 vs SPT6 POB3
minimum = np.min([spt6_pob3_Flag_averaged.quantile([0.05]).min(1),dntd_pob3_Flag_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([spt6_pob3_Flag_averaged.quantile([0.95]).max(1),dntd_pob3_Flag_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Flag}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_Flag_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Flag_heatmap_scale_spt6_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_Flag_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Flag_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[minimum,maximum]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="Greens")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/Flag_heatmaps_colorbar_scale.png', dpi=600)
plt.show();

fc_3 = np.log2(dntd_pob3_Flag_averaged/spt6_pob3_Flag_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Flag_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/Flag_fc_colorbar_scale.png', dpi=600)
plt.show();


## FpR
# ΔNTD POB3 vs SPT6 POB3
minimum = np.min([spt6_pob3_FpR_averaged.quantile([0.05]).min(1),dntd_pob3_FpR_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([spt6_pob3_FpR_averaged.quantile([0.95]).max(1),dntd_pob3_FpR_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{FpR}{FpR}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_FpR_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_spt6_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_FpR_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[minimum,maximum]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="Greens")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/FpR_heatmaps_colorbar_scale.png', dpi=600)
plt.show();

fc_3 = np.log2(dntd_pob3_FpR_averaged/spt6_pob3_FpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/FpR_fc_colorbar_scale.png', dpi=600)
plt.show();



## Myc
# ΔNTD POB3 vs SPT6 POB3
minimum = np.min([spt6_pob3_Myc_averaged.quantile([0.05]).min(1),dntd_pob3_Myc_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([spt6_pob3_Myc_averaged.quantile([0.95]).max(1),dntd_pob3_Myc_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Myc}{Myc}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_Myc_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_spt6_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_Myc_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[minimum,maximum]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="Oranges")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/Myc_heatmaps_colorbar_scale.png', dpi=600)
plt.show();

fc_3 = np.log2(dntd_pob3_Myc_averaged/spt6_pob3_Myc_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/Myc_fc_colorbar_scale.png', dpi=600)
plt.show();


## MpR
# ΔNTD POB3 vs SPT6 POB3
minimum = np.min([spt6_pob3_MpR_averaged.quantile([0.05]).min(1),dntd_pob3_MpR_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([spt6_pob3_MpR_averaged.quantile([0.95]).max(1),dntd_pob3_MpR_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{MpR}{MpR}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_MpR_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_spt6_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_MpR_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[minimum,maximum]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="Oranges")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/MpR_heatmaps_colorbar_scale.png', dpi=600)
plt.show();

fc_3 = np.log2(dntd_pob3_MpR_averaged/spt6_pob3_MpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_fc_dntd-pob3_spt6-pob3.png", dpi=600)
plt.show();

#colorbar
a = np.array([[-1.5,1.5]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="bwr")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
plt.tight_layout()
plt.savefig('plots/heatmaps/scale/MpR_fc_colorbar_scale.png', dpi=600)
plt.show();




## Rpb1
# ΔNTD E154K vs ΔNTD POB3
minimum = np.min([dntd_pob3_Rpb1_averaged.quantile([0.05]).min(1),dntd_e154k_Rpb1_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([dntd_pob3_Rpb1_averaged.quantile([0.95]).max(1),dntd_e154k_Rpb1_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Rpb1}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_e154k_Rpb1_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_dntd_e154k.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[minimum,maximum]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="Greens")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/Rpb1_heatmaps_colorbar_scale.png', dpi=600)
# plt.show();

fc_3 = np.log2(dntd_e154k_Rpb1_averaged/dntd_pob3_Rpb1_averaged)
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
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Rpb1_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[-1.5,1.5]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="bwr")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/Rpb1_fc_colorbar_scale.png', dpi=600)
# plt.show();




## Flag
# ΔNTD E154K vs ΔNTD POB3
minimum = np.min([dntd_pob3_Flag_averaged.quantile([0.05]).min(1),dntd_e154k_Flag_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([dntd_pob3_Flag_averaged.quantile([0.95]).max(1),dntd_e154k_Flag_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Flag}{Flag}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_Flag_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Flag_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_e154k_Flag_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Flag_heatmap_scale_dntd_e154k.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[minimum,maximum]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="Greens")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/Flag_heatmaps_colorbar_scale.png', dpi=600)
# plt.show();

fc_3 = np.log2(dntd_e154k_Flag_averaged/dntd_pob3_Flag_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Flag_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[-1.5,1.5]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="bwr")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/Flag_fc_colorbar_scale.png', dpi=600)
# plt.show();


## FpR
# ΔNTD E154K vs ΔNTD POB3
minimum = np.min([dntd_pob3_FpR_averaged.quantile([0.05]).min(1),dntd_e154k_FpR_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([dntd_pob3_FpR_averaged.quantile([0.95]).max(1),dntd_e154k_FpR_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{FpR}{FpR}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_FpR_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_e154k_FpR_averaged, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_dntd_e154k.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[minimum,maximum]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="Greens")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/FpR_heatmaps_colorbar_scale.png', dpi=600)
# plt.show();

fc_3 = np.log2(dntd_e154k_FpR_averaged/dntd_pob3_FpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/FpR_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[-1.5,1.5]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="bwr")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/FpR_fc_colorbar_scale.png', dpi=600)
# plt.show();



## Myc
# ΔNTD E154K vs ΔNTD POB3
minimum = np.min([dntd_pob3_Myc_averaged.quantile([0.05]).min(1),dntd_e154k_Myc_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([dntd_pob3_Myc_averaged.quantile([0.95]).max(1),dntd_e154k_Myc_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{Myc}{Myc}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_Myc_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_e154k_Myc_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_dntd_e154k.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[minimum,maximum]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="Oranges")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/Myc_heatmaps_colorbar_scale.png', dpi=600)
# plt.show();

fc_3 = np.log2(dntd_e154k_Myc_averaged/dntd_pob3_Myc_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/Myc_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[-1.5,1.5]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="bwr")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/Myc_fc_colorbar_scale.png', dpi=600)
# plt.show();


## MpR
# ΔNTD E154K vs ΔNTD POB3
minimum = np.min([dntd_pob3_MpR_averaged.quantile([0.05]).min(1),dntd_e154k_MpR_averaged.quantile([0.05]).min(1)
                  ])
maximum = np.max([dntd_pob3_MpR_averaged.quantile([0.95]).max(1),dntd_e154k_MpR_averaged.quantile([0.95]).max(1)
                  ])

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('$\\frac{MpR}{MpR}$ ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_MpR_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_dntd_pob3.png", dpi=600)
plt.show();

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_e154k_MpR_averaged, cmap=sns.color_palette("Oranges", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('SPT6\nIAA', fontsize=22)    

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_dntd_e154k.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[minimum,maximum]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="Oranges")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical')
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/MpR_heatmaps_colorbar_scale.png', dpi=600)
# plt.show();

fc_3 = np.log2(dntd_e154k_MpR_averaged/dntd_pob3_MpR_averaged)
minimum = -1.5
maximum = 1.5

f = plt.figure(figsize=(1.75, 3.4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc_3, cmap=sns.color_palette("bwr", as_cmap=True), 
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.axvline(25, color='black', linestyle='--', linewidth=1, zorder=1)
plt.axvline(125, color='black', linestyle='--', linewidth=1, zorder=1)
#plt.title('spt6ΔNTD\nPOB3', fontsize=22, style='italic')

f.tight_layout()
plt.savefig("plots/heatmaps/scale/MpR_heatmap_scale_fc_dntd-e154k_dntd-pob3.png", dpi=600)
plt.show();

#colorbar
# a = np.array([[-1.5,1.5]])
# plt.figure(figsize=(1,4), dpi=600)
# img = plt.imshow(a, cmap="bwr")
# plt.gca().set_visible(False)
# plt.colorbar(orientation='vertical', ticks=[-1.5,0,1.5])
# plt.tight_layout()
# plt.savefig('plots/heatmaps/scale/MpR_fc_colorbar_scale.png', dpi=600)
# plt.show();




