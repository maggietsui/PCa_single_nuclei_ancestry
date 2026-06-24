import anndata
import numpy as np
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import glob

#print("CUDA is available: ")
#print(torch.cuda.is_available())
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

########################################################
# Read in all datasets
########################################################

# Read all anndata (cellbendered, non-QCed, scrublet) into this list
anndata_list = []

for file in glob.glob("/wynton/scratch/mtsui/AA_integration/h5ad_cb_w_exp_cellcount/*_scrublet.h5ad"):
#    if file in remove:
#        continue
    adata = sc.read(filename=file)
    adata.var_names_make_unique()
    print(file, " # cells, # genes: ", adata.shape)
    anndata_list.append(adata)
    
    
# Concatenate all the anndatas
adata = anndata.concat(anndata_list, join="inner", label=None)

# Remove some more
#adata = adata[~adata.obs['PatientID'].isin(['PCa12N','PCa14N','PCa26T','PRO506'])]
adata.obs_names_make_unique()

# Add batch col (by lane)
adata.obs['orig.ident'] = adata.obs['orig.ident'].astype('category')
adata.obs['batch'] = adata.obs['orig.ident'].cat.codes

print("# cells, # genes after merge:", adata.shape)

########################################################
# Filtering
########################################################
# Filter the doublets
adata = adata[adata.obs['predicted_doublet']==False]

# renove batches with <10 cells
counts=adata.obs['batch'].value_counts()
remove = counts[counts<10].index
adata = adata[~adata.obs['batch'].isin(remove)]

# Standard filtering
# (nFeature_RNA > 500) & (nCount_RNA >1000) & (nFeature_RNA < 7500) & (nCount_RNA < 50000)
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_cells(adata, max_genes=7500)
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, max_counts=50000)

# Filtered out all genes with <100 cell expression.
sc.pp.filter_genes(adata, min_cells=100)

# remove mitochondrial genes
# since single nuc should not capture them
# and if it does they're probably noisy
adata = adata[:,~adata.var_names.str.startswith("MT-")]
adata = adata[:,~adata.var_names.str.startswith("ENSG")]

print("# cells, # genes after filtering: ", adata.shape)

# save anndata
adata.write_h5ad(
    "/wynton/scratch/mtsui/AA_integration/h5ad_cb_w_exp_cellcount/all_sn_merged_filtered.h5ad",
    compression="lzf"
)

#adata = sc.read(filename="/wynton/scratch/mtsui/AA_integration/h5ad_cb_w_exp_cellcount/all_sn_merged_filtered.h5ad")
########################################################
# Preprocessing
########################################################

# Before normalizing, save cellbender raw into counts
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.raw = adata # keep full dimension safe

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000, # Use 2000, it errored with 3000
    layer="counts",
    batch_key="batch",
    span=1,
    subset=False,
)
########################################################
# SCVI integration
########################################################

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")


#model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
#model = scvi.model.SCVI(adata, n_layers=2, n_latent=10, gene_likelihood="nb")
model = scvi.model.SCVI(adata, n_layers=2, n_latent=10, gene_likelihood="nb")
#model.train(accelerator='gpu',devices=1,use_gpu = True)
model.train()
model.save("/wynton/scratch/mtsui/AA_integration/h5ad_cb_w_exp_cellcount/scvi_model/")

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

########################################################
# Clustering and UMAP
########################################################
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)


sc.tl.leiden(adata,resolution = 0.5)


# save anndata
adata.write_h5ad(
    "/wynton/scratch/mtsui/AA_integration/h5ad_cb_w_exp_cellcount/all_sn_integrated_nlayers2_nlatent10_final_demux.h5ad",
    compression="lzf"
)

########################################################
# DEGs
########################################################
#adata = sc.read(filename="/wynton/scratch/mtsui/AA_integration/h5ad_new/all_sn_integrated_nlayers2_nlatent10_final.h5ad")
# get DEGs
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

# save anndata
adata.write_h5ad(
    f"/wynton/scratch/mtsui/AA_integration/h5ad_cb_w_exp_cellcount/all_sn_integrated_nlayers2_nlatent10_final_demux.h5ad",
    compression="lzf"
)