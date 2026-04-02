#!/usr/bin/env python3

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pyscenic
import anndata as ad
import pandas as pd
import warnings
import numpy as np
import loompy as lp
import matplotlib as mpl
warnings.filterwarnings("ignore")
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize

import json
import zlib
import base64

lf = lp.connect("data/SCENIC/aucell_10kbpdown.loom", mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
lf.close()

adata= sc.read_h5ad("data/SCENIC/combined_gex_scenic.h5ad")
rss_cellType = regulon_specificity_scores( auc_mtx, adata.obs['celltype'])
rss_cellType
cats = sorted(list(set(adata.obs['celltype'])))

topreg = []
for i,c in enumerate(cats):
    topreg.extend(
        list(rss_cellType.T[c].sort_values(ascending=False)[:5].index)
    )
topreg = list(set(topreg))

auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center', rotation=30)
    return f

colors = sns.color_palette('bright',n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap = [ colorsd[x] for x in adata.obs['celltype'] ]

sns.set()
sns.set(font_scale=1)
fig = palplot( colors, cats, size=1.0)
plt.savefig("MDS_cellType-heatmap-legend-top5.pdf", dpi=600, bbox_inches = "tight")
plt.close()

sns.set(font_scale=1.2)
g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
    yticklabels=False, vmin=-2, vmax=6, row_colors=colormap,
    cmap="YlGnBu", figsize=(21,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')    
g.ax_heatmap.set_xlabel('')
plt.savefig("MDS_heatmap.png", dpi=600, bbox_inches = "tight")
plt.close()

