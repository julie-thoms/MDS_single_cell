#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 19:14:05 2023

@author: prisb
"""
import seaborn as sns

############ Time and patient plots #################3

tcells=['CD4+ T cells','CD8+ non-effector T cells','Treg','CD8+ effector T cells', 'Proliferative T cells']
others=["B cells", "B cell progenitors", "Stromal cells (CXCL12 and LEPR)", "na",  'Low quality High MT Low ribo',]
hsccells=list(set(combined_gex.obs['celltype'].value_counts().index.tolist()) - set(tcells) - set(others))

patientlist = ['P12', 'P11', 'P09', 'P08', 'P17', 'P02', 'P01', 'C4', 'P18', 'C2', 'C1', 'C5', 'C3', 'P03']
timepoints = ['C1D1', 'C1D8', 'C6D8', 'C7D1', 'C7D22', 'Progression', 'C12D29', 'Healthy']

##### Time based UMAPs per patient ##########

patientpath="/home/prisb/projects/henry_hashtag_experiment/Onedrive/MDS_AML_collaboration/patient_timepoint_UMAPS/"
for p in patientlist:
    fig, axs = plt.subplots(1, 8, sharex=True, sharey=True, figsize=(20,4))
    fig.suptitle(p, fontsize=14)
    for i, tp in enumerate(timepoints):
        ax = axs.ravel()[i]
        sc.pl.umap(combined_gex[combined_gex.obs["patient"].isin([p])], color="timepoint",groups=[tp], ax=ax, show=False, title=tp)
        leg = ax.legend()
        leg.remove()
        ax.set_xlabel('')
        ax.set_ylabel('')
    fig.tight_layout()
    fig.savefig(patientpath+p+"_timepoints.png")
    plt.close()
    


################# TIMEPOINTS #####################
timepoints = ['C1D1', 'C1D8', 'C6D8', 'C7D1', 'C7D22', 'C12D29', 'Progression', 'Healthy']
combined_gex.obs['celltype_group'] = ''
combined_gex.obs.loc[combined_gex.obs['celltype'].isin(tcells), 'celltype_group'] = 'T cells'
combined_gex.obs.loc[combined_gex.obs['celltype'].isin(hsccells), 'celltype_group'] = 'HSC cells'

#get a dataframe for timepoint
df = combined_gex.obs[['celltype', 'celltype_group', 'timepoint', 'patient']].copy()
df['c'] = 1 

#groupby
counts = df.groupby(list(df.columns[:-1])).sum()['c']
counts = counts[counts != 0]

tcells=['CD4+ T cells','CD8+ non-effector T cells','Treg','CD8+ effector T cells', 'Proliferative T cells']
others=["B cells", "B cell progenitors", "Stromal cells (CXCL12 and LEPR)", "na",  'Low quality High MT Low ribo',]
hsccells=list(set(combined_gex.obs['celltype'].value_counts().index.tolist()) - set(tcells) - set(others))

#counts dataframe
tcell_counts = counts.unstack(1, fill_value=0)['T cells'].unstack(0, fill_value=0)[tcells]
hsccell_counts = counts.unstack(1, fill_value=0)['HSC cells'].unstack(0, fill_value=0)[hsccells]

# fractions dataframe
tcell_fracs = ((tcell_counts.T / tcell_counts.sum(axis=1).values).T.loc[timepoints]).stack().reset_index()
hscell_fracs = ((hsccell_counts.T / hsccell_counts.sum(axis=1).values).T.loc[timepoints]).stack().reset_index()

#remove unused categories
tcell_fracs['celltype'] = tcell_fracs['celltype'].cat.remove_unused_categories()
hscell_fracs['celltype'] = hscell_fracs['celltype'].cat.remove_unused_categories()

#generate boxplots

fig, ax = plt.subplots(figsize=(25, 10)); sns.boxplot(data=tcell_fracs, x='timepoint', y=0, hue='celltype', boxprops={'alpha': 0.4},ax=ax); sns.stripplot(data=tcell_fracs, x='timepoint', y=0, hue='celltype', dodge=True,size=8, ax=ax)
fig, ax = plt.subplots(figsize=(25, 10)); sns.boxplot(data=hscell_fracs, x='timepoint', y=0, hue='celltype', boxprops={'alpha': 0.4},ax=ax); sns.stripplot(data=hscell_fracs, x='timepoint', y=0, hue='celltype', dodge=True,size=5, ax=ax)

############## RESPONSE ######################
outcomec6d28 = ['Non-Responder', 'Responder', 'Healthy']
outcomec12d29 = ['Responder', 'Progression', 'Healthy']

#get a dataframe for outcome_C6D28
df = combined_gex.obs[['celltype', 'celltype_group', 'outcome_C6D28', 'patient']].copy()
df['c'] = 1 

#groupby
counts = df.groupby(list(df.columns[:-1])).sum()['c']
counts = counts[counts != 0]

tcells=['CD4+ T cells','CD8+ non-effector T cells','Treg','CD8+ effector T cells', 'Proliferative T cells']
others=["B cells", "B cell progenitors", "Stromal cells (CXCL12 and LEPR)", "na",  'Low quality High MT Low ribo',]
hsccells=list(set(combined_gex.obs['celltype'].value_counts().index.tolist()) - set(tcells) - set(others))

#counts dataframe
tcell_counts = counts.unstack(1, fill_value=0)['T cells'].unstack(0, fill_value=0)[tcells]
hsccell_counts = counts.unstack(1, fill_value=0)['HSC cells'].unstack(0, fill_value=0)[hsccells]

# fractions dataframe
tcell_fracs = ((tcell_counts.T / tcell_counts.sum(axis=1).values).T.loc[outcomec6d28]).stack().reset_index()
hscell_fracs = ((hsccell_counts.T / hsccell_counts.sum(axis=1).values).T.loc[outcomec6d28]).stack().reset_index()

#remove unused categories
tcell_fracs['celltype'] = tcell_fracs['celltype'].cat.remove_unused_categories()
hscell_fracs['celltype'] = hscell_fracs['celltype'].cat.remove_unused_categories()

#generate boxplots

fig, ax = plt.subplots(figsize=(25, 10)); sns.boxplot(data=tcell_fracs, x='outcome_C6D28', y=0, hue='celltype', boxprops={'alpha': 0.4},ax=ax); sns.stripplot(data=tcell_fracs, x='outcome_C6D28', y=0, hue='celltype', dodge=True,size=15, ax=ax)
fig, ax = plt.subplots(figsize=(25, 10)); sns.boxplot(data=hscell_fracs, x='outcome_C6D28', y=0, hue='celltype', boxprops={'alpha': 0.4},ax=ax); sns.stripplot(data=hscell_fracs, x='outcome_C6D28', y=0, hue='celltype', dodge=True,size=10, ax=ax)

#get a dataframe for outcome_C12D29
df = combined_gex.obs[['celltype', 'celltype_group', 'outcome_C12D29', 'patient']].copy()
df['c'] = 1 

#groupby
counts = df.groupby(list(df.columns[:-1])).sum()['c']
counts = counts[counts != 0]

tcells=['CD4+ T cells','CD8+ non-effector T cells','Treg','CD8+ effector T cells', 'Proliferative T cells']
others=["B cells", "B cell progenitors", "Stromal cells (CXCL12 and LEPR)", "na",  'Low quality High MT Low ribo',]
hsccells=list(set(combined_gex.obs['celltype'].value_counts().index.tolist()) - set(tcells) - set(others))

#counts dataframe
tcell_counts = counts.unstack(1, fill_value=0)['T cells'].unstack(0, fill_value=0)[tcells]
hsccell_counts = counts.unstack(1, fill_value=0)['HSC cells'].unstack(0, fill_value=0)[hsccells]

# fractions dataframe
tcell_fracs = ((tcell_counts.T / tcell_counts.sum(axis=1).values).T.loc[outcomec12d29]).stack().reset_index()
hscell_fracs = ((hsccell_counts.T / hsccell_counts.sum(axis=1).values).T.loc[outcomec12d29]).stack().reset_index()

#remove unused categories
tcell_fracs['celltype'] = tcell_fracs['celltype'].cat.remove_unused_categories()
hscell_fracs['celltype'] = hscell_fracs['celltype'].cat.remove_unused_categories()

#generate boxplots

fig, ax = plt.subplots(figsize=(25, 10)); sns.boxplot(data=tcell_fracs, x='outcome_C12D29', y=0, hue='celltype', boxprops={'alpha': 0.4},ax=ax); sns.stripplot(data=tcell_fracs, x='outcome_C12D29', y=0, hue='celltype', dodge=True,size=15, ax=ax)
fig, ax = plt.subplots(figsize=(25, 10)); sns.boxplot(data=hscell_fracs, x='outcome_C12D29', y=0, hue='celltype', boxprops={'alpha': 0.4},ax=ax); sns.stripplot(data=hscell_fracs, x='outcome_C12D29', y=0, hue='celltype', dodge=True,size=10, ax=ax)

