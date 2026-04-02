#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 10:53:45 2023

@author: prisb
"""

#Must call muon cannot do import muon as mu. For some reason it does not make muData objects properly that way

###Functions##
import muon as mu
import muon
import scanpy as sc
import matplotlib.pyplot as plt


sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

combined_gex = sc.read_h5ad("/home/prisb/projects/henry_hashtag_experiment/exmds/combined_gex.h5ad")
combined_protein = sc.read_h5ad("/home/prisb/projects/henry_hashtag_experiment/exmds/combined_protein.h5ad")

#Need this set because the index in the combined_gex was set to CategoricalIndex therefore couldn't be made unique and kept throwing an error.
combined_gex.obs_names = combined_gex.obs_names.tolist()
combined_gex.obs_names_make_unique()
combined_protein.obs_names_make_unique()


mdata = muon.MuData({"rna": combined_gex, "prot": combined_protein})
muon.pl.embedding(mdata, basis="rna:X_umap", color=["ADT.CD34"])

mdataoutfile =("./mdata_rna_protein.h5mu")
mdata.write(mdataoutfile)

##save figures for each antibody on the UMAP
antibody=combined_protein.var_names.tolist()

for i in antibody:
    fig = muon.pl.embedding(mdata, basis="rna:X_umap", color=[i], return_fig=True)
    fig.savefig("./UMAP_"+str(i)+".png", bbox_inches="tight")
    plt.close()
    
    