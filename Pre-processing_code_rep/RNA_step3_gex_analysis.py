#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 20:06:14 2023

@author: prisb
"""
import numpy as np
import scanpy as sc
import pandas as pd
import warnings
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

#general file params
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=200, facecolor='white')
warnings.filterwarnings("ignore")


############### Load saved combined_gex####################
combined_gex = sc.read_h5ad("/home/prisb/projects/henry_hashtag_experiment/exmds/combined_gex.h5ad")


########### Investigative plots ##################

#Visualize individual experiments
dlist=combined_gex.obs.dataset.unique().tolist()
for i in dlist:
    fig = sc.pl.umap(combined_gex, color="dataset", groups=[i], return_fig=True)
    fig.savefig("/home/prisb/projects/henry_hashtag_experiment/figures/17_11_2023/"+str(i)+".png", bbox_inches="tight")
    plt.close()
    
#Visualize individual clusters
llist = combined_gex.obs.leiden.unique().tolist()
for i in llist:
    fig = sc.pl.umap(combined_gex, color="leiden", groups=[i], return_fig=True)
    fig.savefig("/home/prisb/projects/henry_hashtag_experiment/figures/17_11_2023/leiden"+str(i)+".png", bbox_inches="tight")
    plt.close()

# visualize marker genes
sc.pl.umap(combined_gex, color=['CD3D','CD3E','TRAC','CD2']) #T-cells
sc.pl.umap(combined_gex, color=['CD34','CD38','KIT','PROM1','FLT3']) #HSPC
sc.pl.umap(combined_gex, color=['MS4A1','CD79A','CD83']) #Bcells
sc.pl.umap(combined_gex, color=['FCER1A','CST3']) #Dendritic cells
sc.pl.umap(combined_gex, color=['FCGR3A']) #Monocytes
sc.pl.umap(combined_gex, color=['AZU1']) #Neutrophils

#Violin plots for cluster 3
sc.pl.rank_genes_groups_violin(combined_gex, groups='3', gene_names=['GZMK','GZMA','GZMB','GZMM'])
with rc_context({'figure.figsize':(9,3)}):
    sc.pl.violin(combined_gex, ['GZMB'], groupby='leiden')
    sc.pl.violin(combined_gex, ['GZMA'], groupby='leiden')
    sc.pl.violin(combined_gex, ['GZMK'], groupby='leiden')
    
sc.pl.rank_genes_groups_violin(combined_gex, groups='3', gene_names=['GZMK','GZMA','GZMB','GZMM'])

#New cluster names

new_cluster_names = {
    "0":"CD4+ T cells",
    "1":"TBC (Low quality?)",
    "2":"CD8+ non-effector T cells",
    "3":"CD8+ non-effector T cells",
    "4":"CD8+ effector T cells",
    "5":"DC progenitor",
    "6":"CD8+ effector T cells",
    "7":"CD8+ non-effector T cells",
    "8":"MEP",
    "9":"CD4+ T cells",
    "10":"TBC (Low quality?)",
    "11":"Myeloid progenitor",
    "12":"Myeloid progenitor",
    "13":"Low quality High MT Low ribo",
    "14":"HSC",
    "15":"MEP",
    "16":"LMPP",
    "17":"Myeloid progenitor",
    "18":"MEP",
    "19":"Granulocyte progenitor",
    "20":"Neutrophil progenitor",
    "21":"MEP",
    "22":"Treg",
    "23":"Proliferative progenitor",
    "24":"LMPP",
    "25":"MEP",
    "26":"MEP",
    "27":"CD8+ effector T cells",
    "28":"B cells",
    "29":"Erythroblast (RBC precursor)",
    "30":"Neutrophil progenitor",
    "31":"B cell progenitors",
    "32":"Proliferative T cells",
    "33":"Stromal cells (CXCL12 and LEPR)"
    }


combined_gex.obs['celltype'] = combined_gex.obs.leiden.astype("str").values
combined_gex.obs.celltype = combined_gex.obs.celltype.replace(new_cluster_names)
combined_gex.obs.celltype = combined_gex.obs.celltype.astype("category")

combined_gex[combined_gex.obs.leiden == "1"].obs.dataset.value_counts()
combined_gex[combined_gex.obs.leiden == "5"].obs.dataset.value_counts()
combined_gex[combined_gex.obs.leiden == "21"].obs.dataset.value_counts()

## Colour by broad category
## number each celltype and plot on UMAP which correspond to legend
## Why is cluster 23 a B-cell progenitor?
## CD34+ multi_lineage plot distribution to seperate the cells. Plot differential expression to see what markers are causing the divide?
## Violin plots for all the UMAPS
## Slides - What data is there and what data is not there. A table of patients, gex and citeseq.
## Slides - This UMAP with a coarse grained celltyping.
## Slides - Various gene expression UMAPs.



sns.violinplot(combined_gex[combined_gex.obs["celltype"] == "CD34+ multi-lineage (early progenitor)"].obsm["X_umap"][:, 0])
a = combined_gex[combined_gex.obs["celltype"] == "CD34+ multi-lineage (early progenitor)"]; a1 = a[a.obsm['X_umap'][:, 0] < 25]; a2 = a[a.obsm['X_umap'][:, 0] >= 25]


a = combined_gex[combined_gex.obs["celltype"] == "CD34+ multi-lineage (early progenitor)"]; a.obs['subtype'] = '1'; a.obs.loc[a.obsm['X_umap'][:, 0] < 25, 'subtype'] = '2'
pd.DataFrame(a.uns['rank_genes_groups']['names'])


## Exclude genes that start with RPS and RPL. Renormalize. 

adata = combined_gex.raw.to_adata()
np.expm1(adata.X.data)
np.expm1(adata.X.data).max()
adata.X.data[:] = np.expm1(adata.X.data)
adata.X
adata.X.sum(axis=1)
adata_noribo = adata[:, adata.var_names[~(adata.var_names.str.startswith('RPS') | adata.var_names.str.startswith('RPL'))]]
adata_noribo.X.sum(axis=1)
import matplotlib.pyplot as plt
plt.violinplot(adata_noribo.X.sum(axis=1))
adata_noribo.X.sum(axis=1)
np.asarray(adata_noribo.X.sum(axis=1))[:, 0]
a = np.asarray(adata_noribo.X.sum(axis=1))[:, 0]
plt.violinplot(a)
sc.pp.normalize_total(adata_noribo, target_sum=1e4)
a = np.asarray(adata_noribo.X.sum(axis=1))[:, 0]
plt.violinplot(a)
a
sc.pp.log1p(adata_noribo)
adata_noribo.X.sum(axis=1)
adata_noribo.X.data[:] = np.expm1(adata_noribo.X.data)
adata_noribo.X.sum(axis=1)
np.log1p(adata_noribo.X.data)
adata_noribo.X.data[:] = np.log1p(adata_noribo.X.data)
adata_noribo.X.sum(axis=1)
sc.tl.rank_genes_groups(adata_noribo, 'celltype')

adata_noribo.uns
adata_noribo.uns.keys()
adata_noribo.uns['rank_genes_groups']
adata_noribo.uns['rank_genes_groups']['names']
pd.DataFrame(adata_noribo.uns['rank_genes_groups']['names'])
pd.DataFrame(adata_noribo.uns['rank_genes_groups']['names']).to_csv('/home/prisb/projects/henry_hashtag_experiment/exmds/24_11_2023.rank_genes_groups.csv')
adata_noribo


### Transfer label #####
oldgex = sc.read_h5ad("/home/prisb/projects/henry_hashtag_experiment/exmds/old_files/combined_gex.h5ad")

testdf = oldgex.obs.loc[:,['barcode','dataset','celltype']]
combined_gex.obs["celltype"] = "na"

gexdf = combined_gex.obs.loc[:,["barcode","dataset","celltype"]]

for i in range(len(gexdf)):
    barcode=gexdf.iloc[i]["barcode"]
    ds=gexdf.iloc[i]["dataset"]
    print(str(i)+"\t"+barcode+"\t"+ds)
    if (testdf.loc[(testdf["barcode"] == barcode) & (testdf["dataset"] == ds)]).empty == True:
        pass
    else:
        ct=testdf.loc[(testdf["barcode"] == barcode) & (testdf["dataset"] == ds)]["celltype"][0]
        combined_gex.obs.loc[(combined_gex.obs.index == barcode) & (combined_gex.obs['dataset'] == ds),'celltype'] = ct
    
  
#tidying up the metadata    
combined_gex.obs['timepoint'] = combined_gex.obs['timepoint'].replace("C12D9", "C12D29")


#simple DE for responder and non-responders
#outcome_C6D28
celltypelist=['CD8+ non-effector T cells', 'CD4+ T cells', 'CD8+ effector T cells', 'Low quality High MT Low ribo', 'Treg', 'Proliferative T cells', 'B cell progenitors', 'Erythroblast (RBC precursor)', 'Granulocyte progenitor', 'TBC (Low quality?)', 'Neutrophil progenitor', 'Myeloid progenitor', 'MEP', 'Proliferative progenitor', 'LMPP', 'B cells', 'na', 'Stromal cells (CXCL12 and LEPR)', 'HSC', 'DC progenitor']

for k in celltypelist:
    adata=combined_gex[combined_gex.obs["celltype"].isin([k])]
    sc.tl.rank_genes_groups(adata, 'outcome_C6D28', groups=['Responder'], reference="Non-Responder", use_raw=False)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame(
        {group + '_' + key: result[key][group]
        for group in groups for key in ['names', 'pvals','logfoldchanges']})
    df.to_csv("/home/prisb/projects/henry_hashtag_experiment/Onedrive/MDS_AML_collaboration/outcomeC6D28_DE/RespondervsNonResponder_"+k+".csv", sep=",")
    del adata, result, groups, df


#outcome_C12D29
celltypelist=['CD8+ non-effector T cells', 'CD4+ T cells', 'CD8+ effector T cells', 'Low quality High MT Low ribo', 'Treg', 'Proliferative T cells', 'B cell progenitors', 'Erythroblast (RBC precursor)', 'Granulocyte progenitor', 'TBC (Low quality?)', 'Neutrophil progenitor', 'Myeloid progenitor', 'MEP', 'Proliferative progenitor', 'LMPP', 'B cells', 'na', 'Stromal cells (CXCL12 and LEPR)', 'HSC', 'DC progenitor']

for k in celltypelist:
    adata=combined_gex[combined_gex.obs["celltype"].isin([k])]
    sc.tl.rank_genes_groups(adata, 'outcome_C12D29', groups=['Responder'], reference="Progression", use_raw=False)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    df = pd.DataFrame(
        {group + '_' + key: result[key][group]
        for group in groups for key in ['names', 'pvals','logfoldchanges']})
    df.to_csv("/home/prisb/projects/henry_hashtag_experiment/Onedrive/MDS_AML_collaboration/outcomeC12D29_DE/RespondervsProgression_"+k+".csv", sep=",")
    del adata, result, groups, df   
    
### Accessing the gene expression values ###

cd34 = combined_gex[combined_gex[:,'CD34'].X > 0.1, :]

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


#### Differential Gene Expression with pyDEseq2 #####

from pydeseq2.ds import DeseqStats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference

adata = combined_gex.raw.to_adata()

#outcome_C6D28
for ct in adata.obs.celltype.unique():
    cell_subset = adata[adata.obs['celltype'] == ct]

    pbs = []
    for p in cell_subset.obs.patient.unique():
        samp_cell_subset = cell_subset[cell_subset.obs['patient'] == p]
            
        rep_adata = sc.AnnData(X = samp_cell_subset.X.sum(axis = 0),
                               var = samp_cell_subset.var[[]])
        
        rep_adata.obs_names = [p]
        rep_adata.obs['condition'] = samp_cell_subset.obs['outcome_C6D28'].iloc[0]
        
        pbs.append(rep_adata)
        
    pb = sc.concat(pbs)
    
    counts = pd.DataFrame(pb.X, columns = pb.var_names, index=pb.obs_names)  
    metadata=pb.obs
    
    #removes samples with a NaN in the condition, pyDEseq2 will throw an error if not
    samples_to_keep = ~metadata.condition.isna()
    counts_df = counts.loc[samples_to_keep]
    metadata = metadata.loc[samples_to_keep]
    
    #removes genes that have less than 10 counts
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    counts_df = counts_df[genes_to_keep]
    
    
    inference = DefaultInference(n_cpus=4)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        inference=inference)
    
    dds.deseq2()

    
    stat_res = DeseqStats(dds, 
                contrast=['condition','Responder','Non-Responder'],
                inference=inference)
    stat_res.summary()
    stat_res.results_df.to_csv('./outcome_C6D28/RespondervsNonResponder_'+str(ct)+'.csv')

del ct, pbs, p, pb, cell_subset, samp_cell_subset, rep_adata, counts, counts_df, metadata, samples_to_keep, genes_to_keep, inference, dds, stat_res

#outcome_C12D29
for ct in adata.obs.celltype.unique():
    cell_subset = adata[adata.obs['celltype'] == ct]

    pbs = []
    for p in cell_subset.obs.patient.unique():
        samp_cell_subset = cell_subset[cell_subset.obs['patient'] == p]
            
        rep_adata = sc.AnnData(X = samp_cell_subset.X.sum(axis = 0),
                               var = samp_cell_subset.var[[]])
        
        rep_adata.obs_names = [p]
        rep_adata.obs['condition'] = samp_cell_subset.obs['outcome_C12D29'].iloc[0]
        
        pbs.append(rep_adata)
        
    pb = sc.concat(pbs)
    
    counts = pd.DataFrame(pb.X, columns = pb.var_names, index=pb.obs_names)  
    metadata=pb.obs
    
    #removes samples with a NaN in the condition, pyDEseq2 will throw an error if not
    samples_to_keep = ~metadata.condition.isna()
    counts_df = counts.loc[samples_to_keep]
    metadata = metadata.loc[samples_to_keep]
    
    #removes genes that have less than 10 counts
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    counts_df = counts_df[genes_to_keep]
    
    
    inference = DefaultInference(n_cpus=4)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        inference=inference)
    
    dds.deseq2()

    
    stat_res = DeseqStats(dds, 
                contrast=['condition','Responder','Progression'],
                inference=inference)
    stat_res.summary()
    stat_res.results_df.to_csv('./outcome_C12D29/RespondervsProgression_'+str(ct)+'.csv')


