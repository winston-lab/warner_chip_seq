#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 30 16:23:24 2025

@author: jlwarner
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

renamer = {
    "'95_I_8WG16_rep4_S24_midpoints_si'": 'ΔNTD_IAA_8WG16_rep4',
    "'95_I_8WG16_rep2_S16_midpoints_si'": 'ΔNTD_IAA_8WG16_rep2',
    "'95_I_8WG16_rep3_S20_midpoints_si'": 'ΔNTD_IAA_8WG16_rep3',
    "'93_D_8WG16_rep3_S17_midpoints_si'": 'SPT6_DMSO_8WG16_rep3',
    "'93_D_8WG16_rep4_S21_midpoints_si'": 'SPT6_DMSO_8WG16_rep4',
    "'95_D_8WG16_rep3_S19_midpoints_si'": 'ΔNTD_DMSO_8WG16_rep3',
    "'95_D_8WG16_rep4_S23_midpoints_si'": 'ΔNTD_DMSO_8WG16_rep4',
    "'93_D_8WG16_rep2_S13_midpoints_si'": 'SPT6_DMSO_8WG16_rep2',
    "'95_D_8WG16_rep2_S15_midpoints_si'": 'ΔNTD_DMSO_8WG16_rep2',
    "'93_I_8WG16_rep2_S14_midpoints_si'": 'SPT6_IAA_8WG16_rep2',
    "'93_I_8WG16_rep3_S18_midpoints_si'": 'SPT6_IAA_8WG16_rep3',
    "'93_I_8WG16_rep4_S22_midpoints_si'": 'SPT6_IAA_8WG16_rep4',
    "'95_I_Flag_rep3_S32_midpoints_si'": 'ΔNTD_IAA_Flag_rep3',
    "'95_I_Flag_rep2_S28_midpoints_si'": 'ΔNTD_IAA_Flag_rep2',
    "'95_I_Flag_rep4_S36_midpoints_si'": 'ΔNTD_IAA_Flag_rep4',
    "'95_D_Flag_rep4_S35_midpoints_si'": 'ΔNTD_DMSO_Flag_rep4',
    "'93_D_Flag_rep4_S33_midpoints_si'": 'SPT6_DMSO_Flag_rep4',
    "'93_I_Flag_rep4_S34_midpoints_si'": 'SPT6_IAA_Flag_rep4',
    "'95_D_Flag_rep3_S31_midpoints_si'": 'ΔNTD_DMSO_Flag_rep3',
    "'93_D_Flag_rep3_S29_midpoints_si'": 'SPT6_DMSO_Flag_rep3',
    "'93_I_Flag_rep3_S30_midpoints_si'": 'SPT6_IAA_Flag_rep3',
    "'95_D_Flag_rep2_S27_midpoints_si'": 'ΔNTD_DMSO_Flag_rep2',
    "'93_D_Flag_rep2_S25_midpoints_si'": 'SPT6_DMSO_Flag_rep2',
    "'93_I_Flag_rep2_S26_midpoints_si'": 'SPT6_IAA_Flag_rep2',
    "'95_I_V5_rep2_S40_midpoints_si'": 'ΔNTD_IAA_V5_rep2',
    "'93_D_V5_rep4_S45_midpoints_si'": 'SPT6_DMSO_V5_rep4',
    "'93_D_V5_rep3_S41_midpoints_si'": 'SPT6_DMSO_V5_rep3',
    "'95_D_V5_rep3_S43_midpoints_si'": 'ΔNTD_DMSO_V5_rep3',
    "'95_D_V5_rep4_S47_midpoints_si'": 'ΔNTD_DMSO_V5_rep4',
    "'93_D_V5_rep2_S37_midpoints_si'": 'SPT6_DMSO_V5_rep2',
    "'95_D_V5_rep2_S39_midpoints_si'": 'ΔNTD_DMSO_V5_rep2',
    "'95_I_V5_rep4_S48_midpoints_si'": 'ΔNTD_IAA_V5_rep4',
    "'93_I_V5_rep2_S38_midpoints_si'": 'SPT6_IAA_V5_rep2',
    "'95_I_V5_rep3_S44_midpoints_si'": 'ΔNTD_IAA_V5_rep3',
    "'93_I_V5_rep3_S42_midpoints_si'": 'SPT6_IAA_V5_rep3',
    "'93_I_V5_rep4_S46_midpoints_si'": 'SPT6_IAA_V5_rep4',
    }

# plotting correlation
pearson_corr = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/correlation/tab/8WG16_spearman.tab', 
                 sep='\t', header=1)
pearson_corr.rename(columns=renamer,inplace=True)
pearson_corr.rename(index=renamer,inplace=True)
sns.set_style('ticks')
plt.figure(dpi=300)
g = sns.clustermap(pearson_corr, 
               xticklabels=True, yticklabels=True, figsize=(8,8),
               annot=True, annot_kws={'fontweight':'bold'},
               )
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontweight='bold', fontsize=12)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontweight='bold', fontsize=12)
plt.savefig('plots/correlation/8WG16_spearman_corr.png', dpi=300)
plt.show();

pearson_corr = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/correlation/tab/Flag_spearman.tab', 
                 sep='\t', header=1)
pearson_corr.rename(columns=renamer,inplace=True)
pearson_corr.rename(index=renamer,inplace=True)
sns.set_style('ticks')
plt.figure(dpi=300)
g = sns.clustermap(pearson_corr, 
               xticklabels=True, yticklabels=True, figsize=(8,8),
               annot=True, annot_kws={'fontweight':'bold'},
               )
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontweight='bold', fontsize=12)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontweight='bold', fontsize=12)
plt.savefig('plots/correlation/Flag_spearman_corr.png', dpi=300)
plt.show();

pearson_corr = pd.read_csv('/Users/jlwarner/Desktop/chip_seq/transfer/correlation/tab/V5_spearman.tab', 
                 sep='\t', header=1)
pearson_corr.rename(columns=renamer,inplace=True)
pearson_corr.rename(index=renamer,inplace=True)
sns.set_style('ticks')
plt.figure(dpi=300)
g = sns.clustermap(pearson_corr, 
               xticklabels=True, yticklabels=True, figsize=(8,8),
               annot=True, annot_kws={'fontweight':'bold'},
               )
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontweight='bold', fontsize=12)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontweight='bold', fontsize=12)
plt.savefig('plots/correlation/V5_spearman_corr.png', dpi=300)
plt.show();


