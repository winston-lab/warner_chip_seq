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



#Import Data
positions = np.arange(-250,1250,10)

wt_I_Rpb1 = pd.read_csv('/Users/jlwarner/Desktop/transfer/data to plot/tab/93_I_8WG16_averaged_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

mut_I_Rpb1 = pd.read_csv('/Users/jlwarner/Desktop/transfer/data to plot/tab/95_I_8WG16_averaged_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

wt_I_Flag = pd.read_csv('/Users/jlwarner/Desktop/transfer/data to plot/tab/93_I_Flag_averaged_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

mut_I_Flag = pd.read_csv('/Users/jlwarner/Desktop/transfer/data to plot/tab/95_I_Flag_averaged_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

wt_I_V5 = pd.read_csv('/Users/jlwarner/Desktop/transfer/data to plot/tab/93_I_V5_averaged_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

wt_D_V5 = pd.read_csv('/Users/jlwarner/Desktop/transfer/data to plot/tab/93_D_V5_averaged_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

mut_I_V5 = pd.read_csv('/Users/jlwarner/Desktop/transfer/data to plot/tab/95_I_V5_averaged_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)

mut_D_V5 = pd.read_csv('/Users/jlwarner/Desktop/transfer/data to plot/tab/95_D_V5_averaged_si_ratio_scale.tab', 
            sep='\t', header=3, names=positions)



#Metagenes
ticks = np.array([-250, 0, 1000, 1250])
labels = ['-250 bp', 'TSS', 'TES', '+250 bp']
ybottom = -0.1
ytop = 3.0
yticks = np.array((0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0))

plt.figure(dpi=300)
plt.plot(positions, wt_I_Rpb1.mean(), label = 'Spt6', color = 'black')
plt.fill_between(positions, wt_I_Rpb1.quantile([0.25,0.75]).iloc[0,:], wt_I_Rpb1.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.2, color='black')
plt.plot(positions, mut_I_Rpb1.mean(), label = 'Spt6Δ2-238', color = 'red')
plt.fill_between(positions, mut_I_Rpb1.quantile([0.25,0.75]).iloc[0,:], mut_I_Rpb1.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.2, color='red')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized Rpb1 occupancy')
plt.ylim(ybottom, ytop)
plt.yticks(ticks=yticks)
plt.legend(loc='upper left')
plt.savefig('Rpb1_metagene.png', dpi=300)
plt.show();

plt.figure(dpi=300)
plt.plot(positions, wt_I_Flag.mean(), label = 'Spt6', color = 'black')
plt.fill_between(positions, wt_I_Flag.quantile([0.25,0.75]).iloc[0,:], wt_I_Flag.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.2, color='black')
plt.plot(positions, mut_I_Flag.mean(), label = 'Spt6Δ2-238', color = 'red')
plt.fill_between(positions, mut_I_Flag.quantile([0.25,0.75]).iloc[0,:], mut_I_Flag.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.2, color='red')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized Flag occupancy')
plt.ylim(ybottom, ytop)
plt.yticks(ticks=yticks)
plt.legend(loc='upper left')
plt.savefig('Flag_metagene.png', dpi=300)
plt.show();

plt.figure(dpi=300)
plt.plot(positions, wt_D_V5.mean(), label = 'DMSO', color = 'blue')
plt.fill_between(positions, wt_D_V5.quantile([0.25,0.75]).iloc[0,:], wt_D_V5.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.2, color='blue')
plt.plot(positions, wt_I_V5.mean(), label = 'IAA', color = 'orange')
plt.fill_between(positions, wt_I_V5.quantile([0.25,0.75]).iloc[0,:], wt_I_V5.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.2, color='orange')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized V5 occupancy')
plt.ylim(ybottom, 1.0)
plt.yticks(ticks=np.array((0,0.5,1.0)))
plt.legend(loc='upper left')
plt.title(label='Spt6-V5 depletion in Spt6-Flag')
plt.show();

plt.figure(dpi=300)
plt.plot(positions, mut_D_V5.mean(), label = 'DMSO', color = 'blue')
plt.fill_between(positions, mut_D_V5.quantile([0.25,0.75]).iloc[0,:], mut_D_V5.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.2, color='blue')
plt.plot(positions, mut_I_V5.mean(), label = 'IAA', color = 'orange')
plt.fill_between(positions, mut_I_V5.quantile([0.25,0.75]).iloc[0,:], mut_I_V5.quantile([0.25,0.75]).iloc[1,:],
                 alpha=0.2, color='orange')
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized V5 occupancy')
plt.ylim(ybottom, 2.0)
plt.yticks(ticks=np.array((0,0.5,1.0,1.5,2.0)))
plt.legend(loc='upper left')
plt.title(label='Spt6-V5 depletion in Spt6Δ2-238-Flag')
plt.show();



#Scatterplots
x = wt_I_Rpb1.iloc[:,25:125].mean(1)
y = mut_I_Rpb1.iloc[:,25:125].mean(1)
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,25],[-0.5,25], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-0.5,25)
plt.xticks(np.arange(0,30,5), fontsize=16)
plt.xlabel('Spt6', fontsize=18)
plt.ylabel('Spt6Δ2-238', fontsize=18)
plt.yticks(np.arange(0,30,5), fontsize=16)
plt.ylim(-0.5,25)
plt.title('$\\alpha\mathdefault{-Rpb1}$', fontsize=20, pad=10)
plt.tight_layout()
plt.savefig('Rpb1_averaged_scatter.png')
plt.show();


x = wt_I_Flag.iloc[:,25:125].mean(1)
y = mut_I_Flag.iloc[:,25:125].mean(1)
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,15],[-0.5,15], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-0.5,15)
plt.xticks(np.arange(0,20,5), fontsize=16)
plt.xlabel('Spt6', fontsize=18)
plt.ylabel('Spt6Δ2-238', fontsize=18)
plt.yticks(np.arange(0,20,5), fontsize=16)
plt.ylim(-0.5,15)
plt.title('$\\alpha\mathdefault{-Flag}$', fontsize=20, pad=10)
plt.tight_layout()
plt.savefig('Flag_averaged_scatter.png')
plt.show();


x = mut_D_V5.iloc[:,25:125].mean(1)
y = mut_I_V5.iloc[:,25:125].mean(1)
df = pd.DataFrame({'x':x, 'y':y})

plt.figure(figsize=(4,4), dpi=300)
values = np.vstack([x, y])
kernel = scipy.stats.gaussian_kde(values)(values)
plt.plot([-0.5,15],[-0.5,15], color='black', linestyle='--', linewidth=1, zorder=1)
sns.scatterplot(data=df, x='x', y='y', 
                s=15, c=kernel, cmap='inferno', zorder=2, linewidth=0)
plt.xlim(-0.5,15)
plt.xticks(np.arange(0,20,5), fontsize=16)
plt.xlabel('DMSO', fontsize=18)
plt.ylabel('IAA', fontsize=18)
plt.yticks(np.arange(0,20,5), fontsize=16)
plt.ylim(-0.5,15)
plt.title('$\\alpha\mathdefault{-V5}$', fontsize=20, pad=10)
plt.tight_layout()
plt.savefig('V5_averaged_scatter.png')
plt.show();




#Heatmaps
ticks = [0,25,125,149]
labels = ['-250 bp', 'TSS', 'TES', '+250 bp']

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
# plt.savefig("Rpb1_ChIP_WT.png")
# plt.show();


f = plt.figure(figsize=(6,16))
gs = f.add_gridspec(1,2)
f.suptitle('Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=wt_I_Rpb1, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=wt_I_Rpb1.quantile([0.90]).max(1), 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.title('Spt6', fontsize=24)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=mut_I_Rpb1, cmap=sns.color_palette("Greens", as_cmap=True), 
            vmax=wt_I_Rpb1.quantile([0.90]).max(1), 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.title('Spt6Δ2-238', fontsize=24)  
f.tight_layout()
plt.savefig("Rpb1_heatmap.png", dpi=300)

#colorbar
a = np.array([[0,4]])
plt.figure(figsize=(1,12))
img = plt.imshow(a, cmap="Greens")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical')
plt.savefig('Rpb1_colorbar.svg', dpi=300)
plt.show();


#fold-change
fc = np.log2(mut_I_Rpb1/wt_I_Rpb1)

plt.figure(figsize=(3,8), dpi=300)
ax = sns.heatmap(data=fc, cmap=sns.color_palette("vlag", as_cmap=True),
                 vmin=0-fc.quantile([0.95]).max(1), vmax=fc.quantile([0.95]).max(1),
                 yticklabels=False)
for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')

#plt.ylabel('genes')
plt.xticks(ticks=ticks, labels=labels, fontsize=14)
plt.title('Rpb1 ChIP occupancy\nlog2$\\left(\\frac{mutant}{WT}\\right)$', fontsize=16)
plt.tight_layout()
plt.savefig("Rpb1_heatmap_log2fc.png")
plt.show();

plt.figure(figsize=(2,4))
lala = sns.boxplot(data=fc.iloc[:,25:125].mean(1), 
            fliersize=1, linewidth=1, 
            linecolor='black', color = 'Green')
lala.set(xlabel=None)
lala.set(ylabel='log2$\\left(\\frac{mutant}{WT}\\right)$ ChIP occupancy')
#lala.set(ylim=(-1,1))
lala.set(title='Rpb1')
plt.tight_layout()
#plt.savefig('Rpb1_boxplot_log2fc', dpi=300)
plt.show();

plt.figure(figsize=(2,4))
lala = sns.violinplot(data=fc.iloc[:,25:125].mean(1), 
            linewidth=1, inner= 'box',
            linecolor='black', color = 'Green',
            inner_kws={'box_width':6, 'whis_width':2})
lala.set(xlabel=None)
lala.set(ylabel='log2$\\left(\\frac{mutant}{WT}\\right)$ ChIP occupancy')
#lala.set(ylim=(-1,1))
lala.set(title='Rpb1')
plt.tight_layout()
plt.savefig('Rpb1_violin_log2fc', dpi=300)
plt.show();



f = plt.figure(figsize=(6,16))
gs = f.add_gridspec(1,2)
f.suptitle('Flag ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=wt_I_Flag, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=wt_I_Flag.quantile([0.90]).max(1), 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.title('Spt6', fontsize=24)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=mut_I_Flag, cmap=sns.color_palette("Blues", as_cmap=True), 
            vmax=wt_I_Flag.quantile([0.90]).max(1), 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')   
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=20)
plt.title('Spt6Δ2-238', fontsize=24)  
f.tight_layout()
plt.savefig("Flag_heatmap.png", dpi=300)

#colorbar
a = np.array([[0,3.5]])
plt.figure(figsize=(1,12))
img = plt.imshow(a, cmap="Blues")
plt.gca().set_visible(False)
plt.colorbar(orientation='vertical', ticks=[0,1,2,3,3.5])
plt.savefig('Flag_colorbar.svg', dpi=300)
plt.show();

#fold-change
fc = np.log2(mut_I_Flag/wt_I_Flag)

plt.figure(figsize=(3,8), dpi=300)
ax = sns.heatmap(data=fc, cmap=sns.color_palette("vlag", as_cmap=True),
                 vmin=0-fc.quantile([0.95]).max(1), vmax=fc.quantile([0.95]).max(1),
                 yticklabels=False)
for spine in ax.spines.values():
    spine.set(visible=True, lw=2, color='black')

#plt.ylabel('genes')
plt.xticks(ticks=ticks, labels=labels, fontsize=14)
plt.title('Flag ChIP occupancy\nlog2$\\left(\\frac{mutant}{WT}\\right)$', fontsize=16)
plt.tight_layout()
plt.savefig("Flag_heatmap_log2fc.png")
plt.show();


plt.figure(figsize=(2,4))
lala = sns.boxplot(data=fc.iloc[:,25:125].mean(1), 
            color='Blue', fliersize=0, linewidth=1, linecolor='black')
lala.set(xlabel=None)
lala.set(ylabel='log2$\\left(\\frac{mutant}{WT}\\right)$ ChIP occupancy')
lala.set(ylim=(-1,1))
lala.set(title='Flag')
plt.tight_layout()
plt.savefig('Flag_boxplot_log2fc', dpi=300)
plt.show();

plt.figure(figsize=(2,4))
lala = sns.violinplot(data=fc.iloc[:,25:125].mean(1), 
            linewidth=1, inner= 'box',
            linecolor='black', color = 'Blue',
            inner_kws={'box_width':6, 'whis_width':2})
lala.set(xlabel=None)
lala.set(ylabel='log2$\\left(\\frac{mutant}{WT}\\right)$ ChIP occupancy')
#lala.set(ylim=(-1,1))
lala.set(title='Flag')
plt.tight_layout()
plt.savefig('Flag_violin_log2fc', dpi=300)
plt.show();


#Violinplots
wt = wt_I_Rpb1.iloc[:,25:125].mean(1)
mut = mut_I_Rpb1.iloc[:,25:125].mean(1)
scipy.stats.wilcoxon(wt,mut)

df_Rpb1 = pd.DataFrame({'Spt6':wt, 'Spt6Δ2-238':mut})
plt.figure(figsize=(6,3))
sns.violinplot(data=df_Rpb1,
               linecolor='black', color='Green')
plt.ylabel('ChIP occupancy')
plt.title('Rpb1')
plt.savefig('Rpb1_violin.png', dpi=300)
plt.show();


plt.figure(figsize=(6,3))
sns.violinplot(data=df_Rpb1, log_scale=2,
               linecolor='black', color='Green')
plt.ylabel('log2(ChIP occupancy)')
plt.title('Rpb1')
plt.savefig('Rpb1_log2_violin.png', dpi=300)
plt.show();

# df.to_csv('df.csv',index=False,header=True)


wt = wt_I_Flag.iloc[:,25:125].mean(1)
mut = mut_I_Flag.iloc[:,25:125].mean(1)
scipy.stats.wilcoxon(wt,mut)

df_Flag = pd.DataFrame({'Spt6':wt, 'Spt6Δ2-238':mut})
plt.figure(figsize=(6,3))
sns.violinplot(data=df_Flag, log_scale=2,
               linecolor='black', color='Blue')
plt.ylabel('log2(ChIP occupancy)')
plt.title('Flag')
plt.savefig('Flag_violin_log2.png', dpi=300)
plt.show();

# df.to_csv('df.csv',index=False,header=True)


wt = (wt_I_Flag/wt_I_Rpb1).mean(1)
mut = (mut_I_Flag/mut_I_Rpb1).mean(1)

scipy.stats.wilcoxon(wt,mut)

df_FlagvRpb1 = pd.DataFrame({'Spt6':wt, 'Spt6Δ2-238':mut})
plt.figure(figsize=(6,3))
sns.violinplot(data=df_FlagvRpb1, 
               log_scale=2,
               linecolor='black')
plt.ylabel('ChIP occupancy')
plt.title('Flag/Rpb1')
plt.savefig('FlagvRpb1_violin.png', dpi=300)
plt.show();

#df_FlagvRpb1.to_csv('df.csv',index=False,header=True)
