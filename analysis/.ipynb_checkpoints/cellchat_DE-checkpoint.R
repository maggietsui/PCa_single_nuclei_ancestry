## Differential ligand receptor interactions between two groups of patients 
library(CellChat)
library(anndata)
library(ggplot2)
options(future.globals.maxSize = 16000000000)

setwd('/wynton/scratch/mtsui/cellchat')
prefix='immune_stroma'
ea_rds = "/wynton/scratch/mtsui/cellchat/cellchat_objs/cellchat_European-American_immune_stroma_for_cellchat.rds"
aa_rds = "/wynton/scratch/mtsui/cellchat/cellchat_objs/cellchat_African-American_immune_stroma_for_cellchat.rds"
##########################################
# read in the two cellchat objs to compare
cellchat.ea <- readRDS(ea_rds)
cellchat.aa <- readRDS(aa_rds)

# Since AA contains celltypes not found in EA, lift EA object
group.new = levels(cellchat.aa@idents)
cellchat.ea <- liftCellChat(cellchat.ea, group.new)

cellchat.ea <- aggregateNet(cellchat.ea)
cellchat.aa <- aggregateNet(cellchat.aa)

object.list <- list(AC = cellchat.ea, AA = cellchat.aa)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

##########################################
# Plot differential number of interactions across populations
##########################################
png(filename=paste0("plots/",prefix,"_AA_AC_diff_interacts_strength.png"),width=200,height=200,units="mm",res=300)

gg1 <- netVisual_heatmap(cellchat)
# Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

dev.off()

##########################################
# Plot all interaction dotplots
##########################################

png(filename=paste0("plots/",prefix,"_AA_AC_diff_sig_pathways.png"),width=1400,height=200,units="mm",res=300)
gg1 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in AA", angle.x = 45, remove.isolate = T)
# Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in AA", angle.x = 45, remove.isolate = T)
gg1 + gg2

dev.off()

##########################################
# Plot cell populations with significant changes
##########################################

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets


gg <- list()
for (i in 1:length(object.list)) {
    object.list[[i]] = netAnalysis_computeCentrality(
        object = object.list[[i]],
        slot.name = "netP",
        thresh = 0.05
    )
    
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) + xlim(0,1) + ylim(0,1)
}
png(filename=paste0("plots/",prefix,"_AA_EA_celltypes_w_sig_changes.png"),width=250,height=100,units="mm",res=300)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
dev.off()

##########################################
# Plot changes in flow for each pathway
##########################################
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
png(filename=paste0("plots/",prefix,"_AA_EA_info_flow_each_pathway.png"),width=250,height=100,units="mm",res=300)
gg1 + gg2
dev.off()

png(filename=paste0("plots/",prefix,"_AA_AC_top_pathways_bubble.png"),width=250,height=100,units="mm",res=300)
netVisual_bubble(cellchat, sources.use = c("FIB-1","CAF-4","CAF-3",'SM-1','SM-2'), targets.use = c("FIB-1","CAF-1","CAF-2","CAF-4","CAF-3",'SM-1','SM-2'),  comparison = c(1, 2), angle.x = 45)
dev.off()

png(filename=paste0("plots/",prefix,"_AA_AC_BMP_sender_receiver_heatmap.png"),width=250,height=100,units="mm",res=300)
cellchat.ea=netAnalysis_computeCentrality(
        object = cellchat.ea,
        slot.name = "netP",
        thresh = 0.05
    )
netAnalysis_signalingRole_network(cellchat.ea, "PTN", width = 10, height = 2,
                        font.size = 10)
dev.off()

save(cellchat, file = "cellchat_merged_AA_EA_immune_epi.RData")


