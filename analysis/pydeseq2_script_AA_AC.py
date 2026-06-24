import os
import pickle as pkl
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import scanpy as sc
import numpy as np
import argparse

parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("--adatafile")
parser.add_argument("--group")
parser.add_argument("--groupby")
parser.add_argument("--tissue")
parser.add_argument("--out")

# Read arguments from command line
args = parser.parse_args()
 
if args.adatafile:
    adatafile = args.adatafile
if args.group:
    group = args.group
if args.groupby:
    groupby = args.groupby
if args.out:
    out = args.out
if args.tissue:
    tissue = args.tissue
    
# Load adata and get total % afr ancestry
tcell=sc.read(adatafile)
if tissue=='All':
    tcell=tcell[(tcell.obs[groupby]==group)]
else:
    tcell=tcell[(tcell.obs[groupby]==group) & (tcell.obs['Type']==tissue)]

# Look at AA and AC only
tcell=tcell[tcell.obs['Race/ethnicity'].isin(['African-American','Afro-Caribbean'])]

# pseudobulk the sample
tcell_pseudo=sc.get.aggregate(tcell, by="PatientID",layer='counts', func="sum")

# Extract counts and metadata
counts_df = pd.DataFrame(tcell_pseudo.layers['sum'],index = tcell_pseudo.obs.index,columns=tcell_pseudo.var_names)
metadata =tcell.obs.drop_duplicates(subset='PatientID')[['Race/ethnicity','PatientID','Type']]
metadata.set_index("PatientID",inplace=True)
metadata=metadata.loc[counts_df.index]

metadata['RP'] = [0 if i.startswith(("PRO","JAM",'Jam')) else 1 for i in metadata.index]


# Add in covariates
sampinfo=pd.read_csv("/data1/home/mtsui/AA_Jam_single_nuc_analysis/all_sample_metadata.csv",index_col=0)
metadata=metadata.merge(sampinfo, how='left',left_index=True, right_index=True)
metadata=metadata[~metadata['Gleason'].isin(['Negative for malignancy','undetermined'])]
metadata.dropna(inplace=True)

# drop missing vals
counts_df = counts_df[counts_df.index.isin(metadata.index)]


genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
counts_df = counts_df[genes_to_keep]

inference = DefaultInference(n_cpus=2)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors=['Type',"Gleason",'Race/ethnicity'],
    continuous_factors=['Age of Patient'],
    ref_level=["Race/ethnicity",'African-American'], # set reference as AA
    refit_cooks=True,
    inference=inference,
)
dds.deseq2()

# compute stats
stat_res = DeseqStats(dds, inference=inference)
stat_res.summary()
res = stat_res.results_df.dropna()
res.to_csv(f"{out}/{group}_{tissue}_pydeseq2_AA_AC.csv")