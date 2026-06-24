### Celda identifies coordinated gene modules or cell states within a single cell dataset
library(sceasy)
library(reticulate)
use_condaenv('scanpy')
library(Seurat)
library(celda)
library(ggplot2)
options(future.globals.maxSize = 16000000000)
library(scran)
#library(singleCellTK)
#options(buildtools.check = function(action) TRUE)

##### Read in Tumor cell datasets and run inital celda
setwd("/c4/home/mtsui/AA_Jam_single_nuc_analysis/celda_tcell")
h5ad_file = "/c4/home/mtsui/AA_Jam_single_nuc_analysis/celltype_subset/tcells_demux_integrated_demux.h5ad"

### Convert anndata to seurat, then seurat to SCE
SCE_epi_seur = sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")
SCE_epi = sceasy::convertFormat(SCE_epi_seur, from="seurat", to="sce")

rm(SCE_epi_seur)

useAssay <- "counts"
altExpName <- "featureSubset"

saveRDS(SCE_epi, "./SCE_epi.rds")
#SCE_epi = readRDS("/c4/home/mtsui/AA_Jam_single_nuc_analysis/celda/SCE_epi.rds")
SCE_epi <- SCE_epi[rowSums(counts(SCE_epi)) > 0, ]
SCE_epi <-SCE_epi[, colSums(counts(SCE_epi)) > 0]

# Select HVG
hvg = getTopHVGs(SCE_epi, n=3000)
##SCE_epi <- selectFeatures(SCE_epi, minCount = 5, minCell = 20, useAssay = useAssay, altExpName = altExpName)
assay <- SummarizedExperiment::assay(SCE_epi, i = useAssay)
minCount =3
minCell=3
sceSubset <- SCE_epi[Matrix::rowSums(assay >= minCount) >= minCell, ]
sceSubset <- sceSubset[rownames(sceSubset) %in% hvg, ]

# Set the altExp
SingleCellExperiment::altExp(SCE_epi, altExpName) <- sceSubset
S4Vectors::metadata(SCE_epi)[["select_features"]] <- list(
            xClass = "SingleCellExperiment",
            minCount = minCount,
            minCell = minCell,
            useAssay = useAssay,
            altExpName = altExpName)
counts(altExp(SCE_epi)) <- as.matrix(counts(altExp(SCE_epi)))


### Test L (number of modules)
moduleSplit <- recursiveSplitModule(SCE_epi, useAssay = useAssay, altExpName = altExpName, initialL = 10, maxL = 70)
pdf("./GridSearchPerplexity.pdf")
plotGridSearchPerplexity(moduleSplit, altExpName = altExpName, sep = 10)
dev.off()

pdf("./RPC.pdf")
plotRPC(moduleSplit, altExpName = altExpName)
dev.off()

saveRDS(moduleSplit,"./moduleSplit.rds")
saveRDS(SCE_epi,"./SCE_epi.rds")

### Pt2: try K (cluster number)
moduleSplit <- readRDS("./moduleSplit.rds")
L = 50
temp <- subsetCeldaList(moduleSplit, list(L = L))
SCE_epi <- recursiveSplitCell(SCE_epi, useAssay = useAssay, altExpName = altExpName, initialK = 3, maxK = 30, yInit = celdaModules(temp))
pdf("./K_value_GridSearch.pdf")
plotGridSearchPerplexity(SCE_epi)
dev.off()

pdf("./K_value_RPC.pdf")
plotRPC(SCE_epi)
dev.off()

saveRDS(SCE_epi,"./sce_L_k.rds")

SCE_epi = readRDS("./sce_L_k.rds")
SCE_epi = celda_CG(SCE_epi,L = 50, K = 20, sampleLabel=SCE_epi$PatientID)
SCE_epi <- celdaUmap(SCE_epi, useAssay = useAssay, altExpName = altExpName,seed = 12345)

saveRDS(SCE_epi,"./sce_L_k.rds")

pdf("./plots/celda_UMAP.pdf")
plotDimReduceCluster(x = SCE_epi, reducedDimName = "celda_UMAP")
dev.off()

pdf("./plots/celda_all_modules.pdf")
plotDimReduceModule(x = SCE_epi, reducedDimName = "celda_UMAP", rescale = TRUE)
dev.off()


#ta <- featureModuleTable(SCE_epi)
#write.table(as.matrix(ta),'feature_module_table_reordered.tsv',quote = F,sep = '\t',row.names = F)


#for (module_number in 1:50){
#	sce = celda::splitModule(SCE_epi,module = module_number)
#	plotDimReduceModule(
#  	sce,
#  	reducedDimName = 'celda_UMAP',
#  	modules = module_number
#	)
#	ggsave(paste0("Module_",module_number,"_UMAP.pdf"))
#}
