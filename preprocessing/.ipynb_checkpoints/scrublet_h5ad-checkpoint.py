import anndata
import numpy as np
import scanpy as sc
import pandas as pd
import argparse
from cellbender.remove_background.downstream import anndata_from_h5

## Script to save raw cellbendered (and demultiplexed) files to h5ad
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("--file")

# Read arguments from command lineq
args = parser.parse_args()
 
if args.file:
    file = args.file

if "demuxlet" in file:
    # Read cellbendered/demuxletted h5ad
    adata = sc.read(filename=file)
else:
    adata = anndata_from_h5(file)

if "raw_feature_bc_matrix_cellbender" in file:
    sample = file.split("/")[-2]
else:
    sample = file.split("/")[-1].split("_cellbender")[0]

if all(x not in adata.obs.columns for x in ['orig.ident','PatientID']):
    adata.obs['orig.ident'] = sample
    adata.obs['PatientID'] = sample

# remove uns cause we don't need and 
# it causes write issues
del adata.uns

# save anndata
adata.write_h5ad(
    f"/wynton/scratch/mtsui/AA_integration/h5ad_cb_w_exp_cellcount/{sample}.h5ad",
    compression="lzf"
)


