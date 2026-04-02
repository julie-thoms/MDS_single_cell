import numpy as np
import pandas as pd
import scanpy as sc
import warnings

#general file params
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=160, facecolor='white')
warnings.filterwarnings("ignore")

#Functions

#categorising filter 
def demuxlet_cat(x):
    y=x[:3]
    if y == 'SNG':
        return 'Single'
    if y == 'DBL':
        return 'Double'
    if y == 'AMB':
        return 'Ambiguous'
    
#patient assignment
genotypePatient_dict={
    "206902510049_R04C01": "P11",
    "206902510049_R05C01": "P12",
    "207071120017_R01C01":"C5",
    "207071120017_R01C02":"C4",
    "207071120017_R02C01":"P11",
    "207071120017_R03C01":"C1",
    "207071120017_R03C02":"C3",
    "207071120017_R05C01":"P12",
    "207071120017_R05C02":"P08",
    "207071120017_R06C01":"P17",
    "207071120017_R07C01":"P03",
    "207071120017_R08C01":"P18",
    "207071120017_R09C01":"P09",
    "207071120017_R09C02":"P02",
    "207071120017_R11C01":"P01",
    "207071120017_R12C01":"C2"
}

def demuxlet_patient(x):
    y=x[4:]
    if y in genotypePatient_dict:
        return genotypePatient_dict[y]
    else:
        return "ERROR"
    
# GEX only
def automateQC(libname):
    outfilename="/home/prisb/projects/henry_hashtag_experiment/exmds/"+libname+"_rna.h5ad"
    infilemtx="/home/prisb/projects/henry_hashtag_experiment/data/230703_scRNAseq_CITEseq_library/"+libname+"/"     
    adata=sc.read_10x_mtx(infilemtx)
    adata.var_names_make_unique()
    demuxlet=pd.read_table("/home/prisb/projects/henry_hashtag_experiment/data/demuxlet/230703_scRNAseq_CITEseq_library/"+libname+".BEST", sep="\t")
    demuxlet['category'] = demuxlet["BEST"].apply(demuxlet_cat)
    adata=adata[np.isin(adata.obs.index,demuxlet["BARCODE"])]
    adata_obs=pd.DataFrame(adata.obs)
    adata_obs["BARCODE"]=adata.obs.index
    demuxlet_ordered=adata_obs.merge(demuxlet, on = "BARCODE")
    demuxlet_ordered.index=demuxlet_ordered["BARCODE"]
    adata.obs=demuxlet_ordered
    adata.var["mt"] = adata.var_names.str.startswith("MT")
    adata.var["gene_expression"]=adata.var.feature_types == "Gene Expression"
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS","RPL"))
    sc.pl.highest_expr_genes(adata, n_top=20)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo"], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(adata, ["n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"], jitter=0.4, multi_panel=True)
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
    return adata, outfilename    





####### C1 D1 TIMEPOINT #######
#load in matrix
results_file='/home/prisb/projects/henry_hashtag_experiment/exmds/hashtag_c1d1.h5ad'
adata=sc.read_10x_mtx('/home/prisb/projects/henry_hashtag_experiment/data/220817_scRNAseq_hashtag_library/C1D1_repeat/per_sample_outs/count/sample_filtered_feature_bc_matrix/',
                      cache=True,
                      gex_only=False) #keeps other features types
adata.var_names_make_unique()


### Read in the demultiplexing and doublet detecting results
demuxafy = pd.read_table("/home/prisb/projects/henry_hashtag_experiment/data/demuxlet/220817_scRNAseq_hashtag_library/C1D1_repeat/C1D1_demuxlet.BEST", sep="\t")

#categorise the droplets
demuxafy['category'] = demuxafy['BEST'].apply(demuxlet_cat)


### Filter the AnnData object for droplet barcodes
adata = adata[np.isin(adata.obs.index,demuxafy["BARCODE"])]


### Order the demuxafy droplets in same order as the AnnData
adata_obs = pd.DataFrame(adata.obs)
adata_obs['BARCODE'] = adata.obs.index
demuxafy_ordered = adata_obs.merge(demuxafy, on = "BARCODE")
demuxafy_ordered.index = demuxafy_ordered["BARCODE"]

### Add demuxafy data to the AnnData
adata.obs = demuxafy_ordered

protein = adata[:, adata.var["feature_types"] == "Antibody Capture"].copy()
rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
adata.layers["log_transformed"] = np.log1p(adata.X)

#plot log transformed counts of antibody capture colored by demuxlet assignment
sc.pl.scatter(adata, x="patientP12", 
              y="patientP11",
              layers='log_transformed',
              color='category')
rna.var['mt'] = rna.var_names.str.startswith('MT-')
rna.var['gene_expression'] = rna.var.feature_types == 'Gene Expression'
adata.var['hashtag'] = adata.var.feature_types == 'Antibody Capture'
rna.var['ribo'] = rna.var_names.str.startswith(("RPS","RPL"))

#Gene expression
sc.pl.highest_expr_genes(rna, n_top=20, )
sc.pp.filter_cells(rna, min_genes=200)
sc.pp.filter_genes(rna, min_cells=3)
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(rna, 
             ['n_genes_by_counts', 
              'total_counts', 
              'pct_counts_mt',
              'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(rna, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(rna, x='total_counts', y='n_genes_by_counts')
sc.pl.scatter(rna, x='total_counts', y='pct_counts_ribo')

#slice rna to remove cells that have too many 
#mitochondrial genes expressed or too many total counts
#and keep only singlets

rna = rna[rna.obs.n_genes_by_counts < 8000, :]
rna = rna[rna.obs.pct_counts_mt < 20, :]
rna_single=rna[rna.obs.category == 'Single',:]
rna_doublet=rna[rna.obs.category != 'Single',:]

#add patient information to rna_single
rna_single.obs['patient'] = rna_single.obs['BEST'].apply(demuxlet_patient)
rna_single.obs["dataset_name"]="hashtag_C1D1"
rna_single.obs["timepoint"]="C1D1"
rna_single.obs["disease_state"] = "MDS"
rna_single.obs["outcome"] = "not available"
rna_single.obs["specific_outcome"] = "not available"


#write raw data to file
rna_single.write(results_file, compression='gzip')

########### C6 D8 TIME POINT##################
bresults_file='/home/prisb/projects/henry_hashtag_experiment/exmds/hashtag_c6d8.h5ad'

bdata=sc.read_10x_mtx('/home/prisb/projects/henry_hashtag_experiment/data/220817_scRNAseq_hashtag_library/C6D8/per_sample_outs/count/sample_filtered_feature_bc_matrix/',
                      cache=True)

bdata.var_names_make_unique()
bdemuxafy=pd.read_table("/home/prisb/projects/henry_hashtag_experiment/data/demuxlet/220817_scRNAseq_hashtag_library/C6D8/C6D8_demuxlet.BEST", sep="\t")
bdemuxafy['category']=bdemuxafy['BEST'].apply(demuxlet_cat)
bdata = bdata[np.isin(bdata.obs.index,bdemuxafy['BARCODE'])]
bdata_obs=pd.DataFrame(bdata.obs)
bdata_obs['BARCODE']=bdata.obs.index
bdemuxafy_ordered=bdata_obs.merge(bdemuxafy, on='BARCODE')
bdemuxafy_ordered.index=bdemuxafy_ordered['BARCODE']
bdata.obs=bdemuxafy_ordered
bprotein=bdata[:,bdata.var["feature_types"] == "Antibody Capture"].copy()
brna=bdata[:,bdata.var["feature_types"] == "Gene Expression"].copy()
bdata.layers['log_transformed']=np.log1p(bdata.X)
sc.pl.scatter(bdata, x='patientP12', 
              y='patientP11', 
              layers='log_transformed', 
              color='category')
brna.var['mt']=brna.var_names.str.startswith('MT')
brna.var['gene_expression']=brna.var.feature_types == 'Gene Expression'
brna.var['ribo']=brna.var_names.str.startswith(("RPS","RPL"))
sc.pl.highest_expr_genes(brna, n_top=20)
sc.pp.filter_cells(brna, min_genes=200)
sc.pp.filter_genes(brna,min_cells=3)
sc.pp.calculate_qc_metrics(brna, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(brna, ['n_genes_by_counts', 
                    'total_counts', 
                    'pct_counts_mt', 
                    'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)
sc.pl.scatter(brna, x='total_counts',y='pct_counts_mt')
sc.pl.scatter(brna, x='total_counts', y='n_genes_by_counts')
sc.pl.scatter(brna, x='total_counts', y='pct_counts_ribo')
brna=brna[brna.obs.n_genes_by_counts < 9000,:]
sc.pl.scatter(brna, x='total_counts', y='n_genes_by_counts')
brna=brna[brna.obs.pct_counts_mt < 20,:]
sc.pl.scatter(brna, x='total_counts',y='pct_counts_mt')
brna_single=brna[brna.obs.category == 'Single',:]
brna_doublet=brna[brna.obs.category != 'Single',:]

#add patient information to brna_single
brna_single.obs['patient'] = brna_single.obs['BEST'].apply(demuxlet_patient)
brna_single.obs["dataset_name"]="hashtag_C6D8"
brna_single.obs["timepoint"]="C6D8"
brna_single.obs["disease_state"] = "MDS"
brna_single.obs["outcome"] = "not available"
brna_single.obs["specific_outcome"] = "not available"

#Write to file
brna_single.write(bresults_file, compression='gzip')

###### 230703 scRNAseq and CITEseq libraries ############


############# Pool1 ###################

adata,results_file=automateQC("Pool1")
adata=adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="Pool1"

#Adding timepoint for each patient
adata_single.obs["timepoint"]="not available"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P08"].obs.index,"timepoint"] = "C1D1"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P01"].obs.index,"timepoint"] = "C7D1"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P17"].obs.index,"timepoint"] = "C1D1"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P02"].obs.index,"timepoint"] = "Progression"

adata_single.obs["disease_state"] = "MDS"

#Adding outcomes and specific outcomes for each patient
adata_single.obs["outcome_C6D8"] = "not available"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P08"].obs.index,"outcome_C6D8"] = "Non-responder"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P01"].obs.index,"outcome_C6D8"] = "Responder"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P17"].obs.index,"outcome_C6D8"] = "Responder"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P02"].obs.index,"outcome_C6D8"] = "Non-responder"

adata_single.obs["specific_outcome_C6D8"] = "not available"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P08"].obs.index,"specific_outcome_C6D8"] = "No response"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P01"].obs.index,"specific_outcome_C6D8"] = "CRbi"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P17"].obs.index,"specific_outcome_C6D8"] = "CR"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P02"].obs.index,"specific_outcome_C6D8"] = "No response"

adata_single.obs["outcome_C12D29"] = "not available"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P08"].obs.index,"outcome_C12D29"] = "Progression"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P01"].obs.index,"outcome_C12D29"] = "Responder"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P17"].obs.index,"outcome_C12D29"] = "Progression"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P02"].obs.index,"outcome_C12D29"] = "Progression"

adata_single.obs["specific_outcome_C12D29"] = "not available"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P08"].obs.index,"specific_outcome_C12D29"] ="Progression"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P01"].obs.index,"specific_outcome_C12D29"] ="CRh"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P17"].obs.index,"specific_outcome_C12D29"] ="Progression"
adata_single.obs.loc[adata_single[adata_single.obs["patient"] == "P02"].obs.index,"specific_outcome_C12D29"] ="Progression"

adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single

############# Pool2 ##################

adata,results_file=automateQC("Pool2")
adata=adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="Pool2"

obs_list=["timepoint","disease_state","outcome_C6D28","specific_outcome_C6D28","outcome_C12D29","specific_outcome_C12D29"]
for i in obs_list:
    adata_single.obs[str(i)] = "not available"
    
patientlist = adata_single.obs.patient.unique().tolist()
metadata={'P12':["C1D1","MDS","Responder","CR","Responder","CR"], 
    'P01':["C12D9","MDS","Responder","Crbi","Responder","CRh"], 
    'C4':["Healthy","Healthy","Healthy","Healthy","Healthy","Healthy"], 
    'P18':["C1D1","MDS","Responder","CR","Progression","Progression"], 
    'P17':["C7D1","MDS","Responder","CR","Progression","Progression"]}

for i in patientlist:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])

adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j

################ Pool3 #################

adata,results_file = automateQC("Pool3")
adata=adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="Pool3"

for a in obs_list:
    adata_single.obs[str(a)] = "not available"

patientlist = adata_single.obs.patient.unique().tolist()
metadata={
 'P17':["Progression","MDS","Responder","CR","Progression","Progression"],
 'P03':["C1D1","MDS","Responder","CR","Responder",",CR"],
 'P18':["C7D1","MDS","Responder","CR","Progression","Progression"],
 'P12':["C7D1","MDS","Responder","CR","Responder","CR"],
 'C2':["Healthy","Healthy","Healthy","Healthy","Healthy","Healthy"],
 'C1':["Healthy","Healthy","Healthy","Healthy","Healthy","Healthy"]}

for i in metadata:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])


adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j, a

################ Pool5_repeat #################

adata,results_file = automateQC("Pool5_repeat")
adata=adata[adata.obs.n_genes_by_counts < 4500, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="Pool5_repeat"

for a in obs_list:
    adata_single.obs[str(a)] = "not available"

patientlist = adata_single.obs.patient.unique().tolist()
metadata={
 'P03':['C12D29','MDS','Responder','CR','Responder','CR'],
 'P09':['C7D1','MDS','Non-Responder','No response','Responder','CRh'],
 'P02':['C1D1','MDS','Non-Responder','No response','Progression','Progression'],
 'P08':['C7D1','MDS','Non-Responder','No response','Progression','Progression']}


for i in metadata:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])

adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j, a

################ Pool6 #################

adata,results_file = automateQC("Pool6")
adata=adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="Pool6"

for a in obs_list:
    adata_single.obs[str(a)] = "not available"

patientlist = adata_single.obs.patient.unique().tolist()
metadata={
 'P08':["Progression","MDS","Non-Responder","No response","Progression","Progression"],
 'P09':["C12D29","MDS","Non-Responder","No response","Responder","CRh"],
 'P01':["C1D1","MDS","Responder","Crbi","Responder","CRh"],
 'P02':["C7D1","MDS","Non-Responder","No response","Progression","Progression"],
 'C5':["Healthy","Healthy","Healthy","Healthy","Healthy","Healthy"]}


for i in metadata:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])


adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j, a


################ Pool7 #################

adata,results_file = automateQC("Pool7")
adata=adata[adata.obs.n_genes_by_counts < 4000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="Pool7"

for a in obs_list:
    adata_single.obs[str(a)] = "not available"

patientlist = adata_single.obs.patient.unique().tolist()
metadata={
    'C3':["Healthy","Healthy","Healthy","Healthy","Healthy","Healthy"],
    'P18':["C7D22","MDS","Responder","CR",	"Progression","Progression"],
    'P09':["C1D8","MDS","Non-Responder","No response","Responder","CRh"],
    'P03':["C7D22","MDS","Responder","CR","Responder","CR"],
    'P12':["C7D22","MDS","Responder","CR","Responder","CR"]}


for i in metadata:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])

adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j, a

############## HSPC_pool1 GEX only ##################

adata, results_file = automateQC("HSPC_pool1")
adata=adata[adata.obs.n_genes_by_counts < 7000, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="HSPC_pool1"

for a in obs_list:
    adata_single.obs[str(a)] = "not available"

patientlist = adata_single.obs.patient.unique().tolist()
metadata={
  'P03':['C12D29','MDS','Responder','CR','Responder','CR'],
  'P08':['C1D1','MDS','Non-Responder','No response','Progression','Progression'],
  'P01':['C1D1','MDS','Responder','Crbi','Responder','CRh'],
  'P09':['C7D1','MDS','Non-Responder','No response','Responder','CRh'],
  'P18':['Progression','MDS','Responder','CR','Progression','Progression'],
  'P02':['C1D1','MDS','Non-Responder','No response','Progression','Progression'],
  'C4':['Healthy','Healthy','Healthy','Healthy','Healthy','Healthy'],
  'C3':['Healthy','Healthy','Healthy','Healthy','Healthy','Healthy']}


for i in metadata:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])

adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j, a

############## HSPC_pool2 GEX only ##################

adata, results_file = automateQC("HSPC_pool2")
adata=adata[adata.obs.n_genes_by_counts < 7000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="HSPC_pool2"

for a in obs_list:
    adata_single.obs[str(a)] = "not available"

patientlist = adata_single.obs.patient.unique().tolist()
metadata={
 'P17':['C1D1','MDS','Responder','CR','Progression','Progression'],
 'P02':['C7D1','MDS','Non-Responder','No response','Progression','Progression'],
 'P12':['C1D1','MDS','Responder','CR','Responder','CR'],
 'P08':['C7D1','MDS','Non-Responder','No response','Progression','Progression'],
 'P09':['C12D29','MDS','Non-Responder','No response','Responder','CRh'],
 'P01':['C7D1','MDS','Responder','Crbi','Responder','CRh']}


for i in metadata:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])

adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j, a

############## HSPC_pool3_repeat GEX only ##################

adata, results_file = automateQC("HSPC_pool3_repeat")
adata=adata[adata.obs.n_genes_by_counts < 8000, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="HSPC_pool3_repeat"

for a in obs_list:
    adata_single.obs[str(a)] = "not available"

patientlist = adata_single.obs.patient.unique().tolist()
metadata={
 'P03':['C1D1','MDS','Responder','CR','Responder','CR'],
 'P18':['C1D1','MDS','Responder','CR','Progression','Progression'],
 'P17':['C7D1','MDS','Responder','CR','Progression','Progression'],
 'P02':['Progression','MDS','Non-Responder','No response','Progression','Progression'],
 'P12':['C7D1','MDS','Responder','CR','Responder','CR'],
 'P08':['Progression','MDS','Non-Responder','No response','Progression','Progression'],
 'C1':['Healthy','Healthy','Healthy','Healthy','Healthy','Healthy']}


for i in metadata:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])

adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j, a

############## HSPC_pool4 GEX only ##################

adata, results_file = automateQC("HSPC_pool4")
adata=adata[adata.obs.n_genes_by_counts < 8000, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]
adata.obs["category"].value_counts()
adata_single=adata[adata.obs.category == "Single", :]
adata_single.obs['patient'] = adata_single.obs['BEST'].apply(demuxlet_patient)
adata_single.obs["dataset_name"]="HSPC_pool4"

for a in obs_list:
    adata_single.obs[str(a)] = "not available"

patientlist = adata_single.obs.patient.unique().tolist()
metadata={
 'P03':['C7D1','MDS','Responder','CR','Responder','CR'],
 'P09':['C1D1','MDS','Non-Responder','No response','Responder','CRh'],
 'P12':['C12D29','MDS','Responder','CR','Responder','CR'],
 'P01':['C12D29','MDS','Responder','Crbi','Responder','CRh'],
 'P18':['C7D1','MDS','Responder','CR','Progression','Progression'],
 'P17':['Progression','MDS','Responder','CR','Progression','Progression'],
 'C5':['Healthy','Healthy','Healthy','Healthy','Healthy','Healthy'],
 'C2':['Healthy','Healthy','Healthy','Healthy','Healthy','Healthy']}


for i in metadata:
    print(i)
    print(metadata[i])
    for j in range(len(metadata[i])):
        print(obs_list[j])
        print(metadata[i][j])
        adata_single.obs.loc[adata_single[adata_single.obs["patient"] == str(i)].obs.index,obs_list[j]]=str(metadata[i][j])

adata_single.write(results_file, compression="gzip")
del adata, results_file, adata_single, metadata, patientlist, i, j, a
