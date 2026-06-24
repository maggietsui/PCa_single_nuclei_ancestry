import anndata
import numpy as np
import scanpy as sc
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
 
## Script to filter for cells and run scrublet

# Adding optional argument
parser.add_argument("--file")

# Read arguments from command lineq
args = parser.parse_args()
 
if args.file:
    file = args.file

adata = sc.read(filename=file)

adata= adata[adata.obs['cell_probability'] >= 0.95]
# Scrublet - identify doublets
# Note - real doublet rate should be around 40%
#sc.pp.scrublet(adata,expected_doublet_rate=0.20)
sc.pp.scrublet(adata,expected_doublet_rate=0.20,threshold=0.65)

# Drop na cells
adata = adata[(~adata.obs['predicted_doublet'].isna()) & (~adata.obs['doublet_score'].isna())]
adata.obs['predicted_doublet'] = adata.obs['predicted_doublet'].astype(bool)

newfile = file.split(".h5ad")[0] + "_scrublet.h5ad"

# save anndata
adata.write_h5ad(
    newfile,
    compression="lzf"
)

