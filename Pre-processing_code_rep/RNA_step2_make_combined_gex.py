#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 15:16:11 2024

@author: prisb
"""

import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import warnings
import anndata as ad

#general file params
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=200, facecolor='white')
warnings.filterwarnings("ignore")

#Load all libraries
folderpath="/home/prisb/projects/henry_hashtag_experiment/exmds/"

hashtag_c1d1=sc.read_h5ad(folderpath+"hashtag_c1d1.h5ad")
hashtag_c6d8=sc.read_h5ad(folderpath+"hashtag_c6d8.h5ad")
HSPC_pool1=sc.read_h5ad(folderpath+"HSPC_pool1_rna.h5ad")
HSPC_pool2=sc.read_h5ad(folderpath+"HSPC_pool2_rna.h5ad")
HSPC_pool3_repeat=sc.read_h5ad(folderpath+"HSPC_pool3_repeat_rna.h5ad")
HSPC_pool4=sc.read_h5ad(folderpath+"HSPC_pool4_rna.h5ad")
Pool1=sc.read_h5ad(folderpath+"Pool1_rna.h5ad")
Pool2=sc.read_h5ad(folderpath+"Pool2_rna.h5ad")
Pool3=sc.read_h5ad(folderpath+"Pool3_rna.h5ad")
Pool5_repeat=sc.read_h5ad(folderpath+"Pool5_repeat_rna.h5ad")
Pool6=sc.read_h5ad(folderpath+"Pool6_rna.h5ad")
Pool7=sc.read_h5ad(folderpath+"Pool7_rna.h5ad")

adatas ={"hashtag_c1d1":hashtag_c1d1, 
         "hashtag_c6d8":hashtag_c6d8,
         "HSPC_pool1":HSPC_pool1,
         "HSPC_pool2":HSPC_pool2,
         "HSPC_pool3_repeat":HSPC_pool3_repeat,
         "HSPC_pool4":HSPC_pool4,
         "Pool1":Pool1,
         "Pool2":Pool2,
         "Pool3":Pool3,
         "Pool5_repeat":Pool5_repeat,
         "Pool6":Pool6,
         "Pool7":Pool7}

combined_gex = ad.concat(adatas, label="dataset", join="outer")
combined_gex.obs_names_make_unique()

combined_gex.raw = combined_gex


sc.pp.normalize_total(combined_gex, target_sum=1e4)
sc.pp.log1p(combined_gex)
sc.pp.highly_variable_genes(combined_gex, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(combined_gex)
np.sum(combined_gex.var.highly_variable)


sc.tl.pca(combined_gex, svd_solver='arpack')
sc.pl.pca_variance_ratio(combined_gex, log=True)
sc.pp.neighbors(combined_gex, n_neighbors=10, n_pcs=26)
sc.tl.leiden(combined_gex, resolution=.75)
sc.tl.umap(combined_gex, spread=1., min_dist=.5, random_state=11)
sc.pl.umap(combined_gex, color="leiden", legend_loc="on data")
sc.pl.umap(combined_gex, color="dataset")


sc.tl.rank_genes_groups(combined_gex,'leiden', method='t-test_overestim_var')
sc.pl.rank_genes_groups(combined_gex, n_genes=25, sharey=False)
result = combined_gex.uns['rank_genes_groups']
groups = result['names'].dtype.names
df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']})



#### save combined_gex to file#####
outfile = "/home/prisb/projects/henry_hashtag_experiment/exmds/combined_gex.h5ad"

combined_gex.write_h5ad(outfile, compression="gzip")
