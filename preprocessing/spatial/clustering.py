import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import numpy as np
import squidpy as sq
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("--adatafile")
parser.add_argument("--zarr")
# Read arguments from command line
args = parser.parse_args()
 
if args.adatafile:
    adatafile = args.adatafile   
    prefix=adatafile.split("/")[-1].split("_")[0]

if args.zarr:
    zarr = args.zarr   
    prefix=zarr.split("/")[-1].split(".")[0]

# Load in zarr
sdata = sd.read_zarr(zarr)

# Extract adata
adata = sdata.tables["table"]

sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

fig, axs = plt.subplots(1, 4, figsize=(15, 4))

axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)


axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,
    ax=axs[2],
)

axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)

plt.savefig("/data1/home/mtsui/Xenium/plots/EUR_qc_metrics.png", format="png",  bbox_inches="tight")

# Filtering

sc.pp.filter_cells(adata, min_counts=20)
perc_98=np.percentile(adata.obs['total_counts'], 98)
sc.pp.filter_cells(adata, max_counts=perc_98)
sc.pp.filter_cells(adata, min_genes=20)


# Normalization
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)


#regress out covariates
#covariates = [
#    'total_counts',
#    'n_genes_by_counts',
#    'cell_area'
#]

#sc.pp.regress_out(adata, covariates)


# Assign cores and remove cells that were floating between cores
cores = pd.read_csv("/data1/home/mtsui/Xenium/AFR_core_id_coords.csv",skiprows=2)
cores=cores.drop_duplicates()

coords = adata.obsm["spatial"]   # Nx2 array
x = coords[:, 0]
y = coords[:, 1]
adata.obs['core_id'] = np.NaN


for core in set(cores['Selection']):
    df = cores[cores['Selection']==core]
    # Your bounding box
    xmin, xmax = df['X'].min(), df['X'].max()
    ymin, ymax = df['Y'].min(), df['Y'].max()
    
    mask = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
    
    # Cell IDs (index in adata.obs)
    cells_in_box = adata.obs.index[mask]
    
    # Now extract any labels from adata.obs
    adata.obs.loc[cells_in_box, "core_id"] = df['Selection'].iloc[0]
    
adata=adata[~adata.obs['core_id'].isna()]

# Identify doublets in each core
for core in set(adata.obs['core_id']):
    curr_core=adata[adata.obs['core_id']==core]
    curr_core.write_h5ad(
        f"/data1/home/mtsui/Xenium/scrublet/{prefix}_{core}.h5ad",
        compression="lzf"
    )
    sc.pp.scrublet(curr_core,expected_doublet_rate=0.05)
    sc.pl.scrublet_score_distribution(curr_core)
    plt.savefig(f"/data1/home/mtsui/Xenium/plots/{prefix}_{core}_scrublet_distribution.png", format="png",  bbox_inches="tight")
    curr_core.obs['predicted_doublet'] = curr_core.obs['predicted_doublet'].astype(bool)
    curr_core.obs.to_csv(f"/data1/home/mtsui/Xenium/scrublet/{prefix}_{core}_scrublet_scores.csv")


# Remove doublets

import glob
scrub=[]
for file in glob.glob(f"/data1/home/mtsui/Xenium/scrublet/{prefix}*_scrublet_scores.csv"):
    df = pd.read_csv(file,index_col=0)
    scrub.append(df)
scrub=pd.concat(scrub)
scrub = scrub[scrub['cell_id'].isin(adata.obs['cell_id'])]

adata.obs= adata.obs.merge(scrub[['cell_id','predicted_doublet']].set_index("cell_id"), left_on ='cell_id', right_on ='cell_id',how='left')
adata.obs['predicted_doublet']=adata.obs['predicted_doublet'].astype(bool)
adata=adata[adata.obs['predicted_doublet']==False]

# Remove clusters that have median count < 100
median_count = adata.obs.groupby("leiden")['n_counts'].median()
adata.obs.index = adata.obs.index.astype(str)
adata=adata[adata.obs['leiden'].isin(median_count[median_count > 100].index)]

# PCA and clustering
#sc.pp.highly_variable_genes(adata, n_top_genes=2500)
#adata = sc.read(adatafile)
sc.pp.pca(adata, n_comps=100)
sc.pl.pca_variance_ratio(adata, n_pcs=100, log=True)
plt.savefig(f"/data1/home/mtsui/Xenium/plots/{prefix}_pca_variance.png", format="png",  bbox_inches="tight")

sc.pp.neighbors(adata,n_pcs=100)
sc.tl.umap(adata)

sc.tl.leiden(adata, resolution = 0.7)

sc.pl.umap(adata,color=['leiden'])
plt.savefig(f"/data1/home/mtsui/Xenium/plots/{prefix}_umap.png", format="png",  bbox_inches="tight")

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

#adata.write_h5ad(
#    f"/data1/home/mtsui/Xenium/{prefix}_filtered.h5ad",
#    compression="lzf"
#)

