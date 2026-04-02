import scanpy as sc 
import numpy as np
import pandas as pd
import muon as mu
import anndata as ad


mdata = mdata.read_h5mu("temp_15052025.h5mu")

mu.tl.mofa(mdata, outfile="mofa_model_15052025.hdf5", n_factors=15)

exit()