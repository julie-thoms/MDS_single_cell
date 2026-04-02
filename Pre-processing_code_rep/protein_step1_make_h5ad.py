#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 15:37:10 2023

@author: prisb
"""

import numpy as np
import scanpy as sc
import muon as mu
from muon import prot as pt
import seaborn as sns
import anndata as ad

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')


################### HSPC POOL1 ################################
mdata=mu.read_10x_mtx("/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/HSPC_pool1/")
mdata_raw = mu.read_10x_mtx("/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/HSPC_pool1/raw_feature_bc_matrix")

prot = mdata.mod['prot']

mdata_raw['rna'].obs["log10umi"] = np.array(np.log10(mdata_raw['rna'].X.sum(axis=1) + 1)).reshape(-1)
#Manually drawing the histogram because the predefined function in muon did show me enough xticks
dat = mdata_raw['rna']
obs_key = [i for i in ['log10umi'] if i in dat.obs.columns]
df = dat.obs.loc[:, obs_key]
obs_key = [i for i in ['log10umi'] if i in dat.obs.columns]
sns.histplot(df[df >= 0.5],bins=50, kde=True)

#Find the Isotype controls
ADTlist = mdata_raw['prot'].var_names.values
for i in range(len(ADTlist)):
    print(i,ADTlist[i])
isotypes = np.append(mdata_raw['prot'].var_names[43:47].values,mdata_raw['prot'].var_names[120:125].values)

#Preserve original counts
prot.layers['counts'] = prot.X

#DSB normalize
pt.pp.dsb(mdata, mdata_raw, empty_counts_range=(1.0, 3.0), isotype_controls=isotypes, random_state=1)
sc.pl.scatter(mdata['prot'], x="ADT.CD34", y="ADT.CD3", layers='counts')
sc.pl.scatter(mdata['prot'], x="ADT.CD34", y="ADT.CD3")

#Downstream Analysis
sc.tl.pca(prot)
sc.pp.neighbors(prot)
sc.tl.umap(prot, random_state=1)
sc.pl.umap(prot,color=['ADT.CD3','ADT.CD34'])

#write H5AD file
outfile = "/home/prisb/projects/henry_hashtag_experiment/exmds/protein_HSPC_pool1.h5ad"
prot.write_h5ad(outfile, compression="gzip")

########################### HSPC POOL 2 ############################
mdata2=mu.read_10x_mtx("/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/HSPC_pool2/")
mdata_raw2 = mu.read_10x_mtx("/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/HSPC_pool2/raw_feature_bc_matrix/")

prot2 = mdata2.mod['prot']

mdata_raw2['rna'].obs["log10umi"] = np.array(np.log10(mdata_raw2['rna'].X.sum(axis=1) + 1)).reshape(-1)
#Manually drawing the histogram because the predefined function in muon did show me enough xticks
dat2 = mdata_raw2['rna']
obs_key2 = [i for i in ['log10umi'] if i in dat2.obs.columns]
df2 = dat2.obs.loc[:, obs_key2]
sns.histplot(df2[df2 >= 0.5],bins=50, kde=True)

#Find the Isotype controls
ADTlist2 = mdata_raw2['prot'].var_names.values
for i in range(len(ADTlist2)):
    print(i,ADTlist2[i])
isotypes = np.append(mdata_raw2['prot'].var_names[43:47].values,mdata_raw2['prot'].var_names[120:125].values)

#Preserve original counts
prot2.layers['counts'] = prot2.X

#DSB normalize
pt.pp.dsb(mdata2, mdata_raw2, empty_counts_range=(1.0, 3.0), isotype_controls=isotypes, random_state=1)
sc.pl.scatter(mdata2['prot'], x="ADT.CD34", y="ADT.CD38", layers='counts')
sc.pl.scatter(mdata2['prot'], x="ADT.CD34", y="ADT.CD38")

#Downstream Analysis
sc.tl.pca(prot2)
sc.pp.neighbors(prot2)
sc.tl.umap(prot2, random_state=1)
sc.pl.umap(prot2,color=['ADT.CD3','ADT.CD34'])

#write H5AD file
outfile = "/home/prisb/projects/henry_hashtag_experiment/exmds/protein_HSPC_pool2.h5ad"
prot2.write_h5ad(outfile, compression="gzip")

##################### HSPC_pool3_repeat ############
mdata3=mu.read_10x_mtx("/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/HSPC_pool3_repeat/multi/per_sample_outs/")
mdata_raw3 = mu.read_10x_mtx("/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/HSPC_pool3_repeat/multi/raw_feature_bc_matrix/")

prot3 = mdata3.mod['prot']

mdata_raw3['rna'].obs["log10umi"] = np.array(np.log10(mdata_raw3['rna'].X.sum(axis=1) + 1)).reshape(-1)

#Manually drawing the histogram because the predefined function in muon did show me enough xticks
dat3 = mdata_raw3['rna']
obs_key3 = [i for i in ['log10umi'] if i in dat3.obs.columns]
df3 = dat3.obs.loc[:, obs_key3]
sns.histplot(df3[df3 >= 0.5],bins=50, kde=True)

#Find the Isotype controls
ADTlist3 = mdata_raw3['prot'].var_names.values
for i in range(len(ADTlist3)):
    print(i,ADTlist3[i])
isotypes = np.append(mdata_raw3['prot'].var_names[43:47].values,mdata_raw3['prot'].var_names[120:125].values)

#Preserve original counts
prot3.layers['counts'] = prot3.X

#DSB normalize
pt.pp.dsb(mdata3, mdata_raw3, empty_counts_range=(1.0, 3.0), isotype_controls=isotypes, random_state=1)
sc.pl.scatter(mdata3['prot'], x="ADT.CD34", y="ADT.CD38", layers='counts')
sc.pl.scatter(mdata3['prot'], x="ADT.CD34", y="ADT.CD38")

#Downstream Analysis
sc.tl.pca(prot3)
sc.pp.neighbors(prot3)
sc.tl.umap(prot3, random_state=1)
sc.pl.umap(prot3,color=['ADT.CD3','ADT.CD34'])

#write H5AD file
outfile = "/home/prisb/projects/henry_hashtag_experiment/exmds/protein_HSPC_pool3_repeat.h5ad"
prot3.write_h5ad(outfile, compression="gzip")

##################### HSPC_pool4 ############
mdata4=mu.read_10x_mtx("/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/HSPC_pool4/multi/per_sample_outs/")
mdata_raw4 = mu.read_10x_mtx("/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/HSPC_pool4/multi/raw_feature_bc_matrix/")

prot4 = mdata4.mod['prot']

mdata_raw4['rna'].obs["log10umi"] = np.array(np.log10(mdata_raw4['rna'].X.sum(axis=1) + 1)).reshape(-1)

#Manually drawing the histogram because the predefined function in muon did show me enough xticks
dat4 = mdata_raw4['rna']
obs_key4 = [i for i in ['log10umi'] if i in dat4.obs.columns]
df4 = dat4.obs.loc[:, obs_key4]
sns.histplot(df4[df4 >= 0.5],bins=50, kde=True)

#Find the Isotype controls
ADTlist4 = mdata_raw4['prot'].var_names.values
for i in range(len(ADTlist4)):
    print(i,ADTlist4[i])
isotypes = np.append(mdata_raw4['prot'].var_names[43:47].values,mdata_raw4['prot'].var_names[120:125].values)

#Preserve original counts
prot4.layers['counts'] = prot4.X

#DSB normalize
pt.pp.dsb(mdata4, mdata_raw4, empty_counts_range=(1.0, 3.0), isotype_controls=isotypes, random_state=1)
sc.pl.scatter(mdata4['prot'], x="ADT.CD34", y="ADT.CD38", layers='counts')
sc.pl.scatter(mdata4['prot'], x="ADT.CD34", y="ADT.CD38")

#Downstream Analysis
sc.tl.pca(prot4)
sc.pp.neighbors(prot4)
sc.tl.umap(prot4, random_state=1)
sc.pl.umap(prot4,color=['ADT.CD3','ADT.CD34'])

#write H5AD file
outfile = "/home/prisb/projects/henry_hashtag_experiment/exmds/protein_HSPC_pool4.h5ad"
prot4.write_h5ad(outfile, compression="gzip")

################## COMBINED PROTEIN ##########################

prot1 = sc.read_h5ad("./protein_HSPC_pool1.h5ad")
prot1.obs['dataset'] = "HSPC_pool1"


prot2 = sc.read_h5ad('./protein_HSPC_pool2.h5ad')
prot2.obs['dataset'] = "HSPC_pool2"

prot3 = sc.read_h5ad('./protein_HSPC_pool3_repeat.h5ad')
prot3.obs['dataset'] = 'HSPC_pool3_repeat'

prot4 = sc.read_h5ad('./protein_HSPC_pool4.h5ad')
prot4.obs['dataset'] = 'HSPC_pool4'

protdata = {"HSPC_pool1": prot1,
            "HSPC_pool2":prot2,
            'HSPC_pool3_repeat':prot3,
            'HSPC_pool4':prot4}

combined_protein = ad.concat(protdata, label="dataset", join="outer")
combined_protein.obs_names_make_unique()

fileout = "./combined_protein.h5ad"
combined_protein.write_h5ad(fileout, compression="gzip")
