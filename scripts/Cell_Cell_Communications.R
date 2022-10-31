## CellChat   for github

options(stringsAsFactors = F)

library(ggplot2)
library(BuenColors)
library(Seurat)
library(CellChat)
library(patchwork)
library(reshape2)
library(cowplot)
library(pheatmap)
library(viridis)

WAT.s.sub <- readRDS("WAT_scRNA_MTX_sc_20220921.RDS")
WAT.s.sub <- NormalizeData(WAT.s.sub)


cellchat.sc <- createCellChat(object = WAT.s.sub,assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)


# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling


# set the used database in the object
cellchat.sc@DB <- CellChatDB.use


# subset the expression data of signaling genes for saving computation cost
cellchat.sc <- subsetData(cellchat.sc) # This step is necessary even if using the whole database
# future::plan("multiprocess", workers = 4) # do parallel


cellchat.sc <- identifyOverExpressedGenes(cellchat.sc)
cellchat.sc <- identifyOverExpressedInteractions(cellchat.sc)


cellchat.sc <- computeCommunProb(cellchat.sc)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.sc <- filterCommunication(cellchat.sc, min.cells = 10)

df.net <- subsetCommunication(cellchat.sc)

cellchat.sc <- computeCommunProbPathway(cellchat.sc)


cellchat.sc <- aggregateNet(cellchat.sc)


# define the color and group size of scWAT
whole_color_use_sc <- c(colorRampPalette(colors = c("#ecd9f1", "#967bce"))(length(grep("lyC|B_cell",unique(cellchat.sc@meta$ident),value = T))), #ly purple
                        colorRampPalette(colors = c("#C6DBEF", "#075a84"))(length(grep("myC|Mast",unique(cellchat.sc@meta$ident),value = T))), #my bule
                        colorRampPalette(colors = c("#F3E55C", "#E8602D"))(length(grep("sfC",unique(cellchat.sc@meta$ident),value = T))), #fb orange
                        "#C0C0C0", # grey
                        colorRampPalette(colors = c("#a6dbbb", "#359566"))(length(grep("vC",unique(cellchat.sc@meta$ident),value = T))) ##vc green
)
groupSize.sc <- as.numeric(table(cellchat.sc@idents))
netVisual_circle(cellchat.sc@net$weight, 
                 vertex.weight = groupSize.sc, 
                 weight.scale = T, 
                 label.edge= F, vertex.size.max = 5,
                 title.name = "Interaction weights/strength",
                 arrow.size = 0.05, color.use = whole_color_use_sc)








cellchat_df.sc <- cellchat.sc@net$weight

a <- c(intersect(rownames(cellchat_df.sc)[rowSums(cellchat_df.sc)==0],colnames(cellchat_df.sc)[colSums(cellchat_df.sc)==0]))



cellchat_df.sc <- cellchat_df.sc[!(rownames(cellchat_df.sc) %in% a), !(colnames(cellchat_df.sc) %in% a)]

annotation_col <- data.frame(row=names(table(cellchat.sc@idents)),cell_type = "a",row.names = names(table(cellchat.sc@idents)))

annotation_col$cell_type[grep("lyC|B_cell",annotation_col$row)] <- "lyC"
annotation_col$cell_type[grep("myC|Mast",annotation_col$row)] <- "myC"
annotation_col$cell_type[grep("sfC",annotation_col$row)] <- "sfC"
annotation_col$cell_type[grep("adipose",annotation_col$row)] <- "adipose"
annotation_col$cell_type[grep("vC",annotation_col$row)] <- "vC"

annotation_col <- data.frame(cell_type = annotation_col$cell_type,row.names = rownames(annotation_col))

annotation_row <- annotation_col

ann_colors <- list(cell_type = c(lyC = colorRampPalette(colors = c("#ecd9f1", "#967bce"))(3)[2], 
                                 myC = colorRampPalette(colors = c("#C6DBEF", "#075a84"))(3)[2], 
                                 sfC = colorRampPalette(colors = c("#F3E55C", "#E8602D"))(3)[2], 
                                 adipose = "#C0C0C0", 
                                 vC = colorRampPalette(colors = c("#a6dbbb", "#359566"))(3)[2]
)
)


pheatmap(cellchat_df.sc,
         cluster_cols = T,
         cluster_rows = T,
         # cellwidth = 8.5,
         # cellheight = 8.5,
         cellwidth = 7.5,
         cellheight = 7.5,
         border_color = NA,
         fontsize_number = 5,
         cutree_cols = 3,
         cutree_rows = 3,
         clustering_method = "ward.D2",
         # display_numbers = T,
         # angle_col = 45,
         treeheight_col = 20,
         treeheight_row = 20,
         fontsize = 8,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         color = viridis(n = 51, alpha = 1, begin = 0, end = 1, option = "viridis")
)



temp <- cellchat_df.sc[rownames(cellchat_df.sc) %in% c("sfC16_late_CPA","adipose","sfC06_FAP","sfC10_FAP","sfC11_MSL","sfC15_FAP","myC02_LAM","myC15_Mox","myC14_Classical_Mo","vC11_Blood_EC_vein","myC00_M2","myC01_M2","myC13_M2"), 
                       colnames(cellchat_df.sc) %in% c("sfC16_late_CPA","adipose","sfC06_FAP","sfC10_FAP","sfC11_MSL","sfC15_FAP","myC02_LAM","myC15_Mox","myC14_Classical_Mo","vC11_Blood_EC_vein","myC00_M2","myC01_M2","myC13_M2")]

temp_color.sc <- c(colorRampPalette(colors = c("#C6DBEF", "#075a84"))(length(grep("myC|Mast",unique(cellchat.sc@meta$ident),value = T)))[c(1,2,3,14,15,16)], #my bule
                   colorRampPalette(colors = c("#F3E55C", "#E8602D"))(length(grep("sfC",unique(cellchat.sc@meta$ident),value = T)))[c(7,11,12,16,17)], #fb orange
                   "#C0C0C0", # grey
                   colorRampPalette(colors = c("#a6dbbb", "#359566"))(length(grep("vC",unique(cellchat.sc@meta$ident),value = T)))[c(12)]
)

netVisual_circle(temp, 
                 vertex.weight = as.numeric(table(cellchat.sc@idents)[names(table(cellchat.sc@idents)) %in% c("sfC16_late_CPA","adipose","sfC06_FAP","sfC10_FAP","sfC11_MSL","sfC15_FAP","myC02_LAM","myC15_Mox","myC14_Classical_Mo","vC11_Blood_EC_vein","myC00_M2","myC01_M2","myC13_M2")]), 
                 weight.scale = T, 
                 edge.weight.max = max(temp),title.name = "Weight",vertex.size.max = 5,
                 color.use = temp_color.sc
)



for (i in cellchat.sc@netP$pathways) {
  a <- netVisual_aggregate(object = cellchat.sc, 
                           signaling = i, 
                           layout = "circle",
                           vertex.weight = groupSize.sc,
                           vertex.size.max = 5,
                           color.use = whole_color_use_sc)
  aa <- ggdraw(a)
  ggsave(
    filename = paste0("sc_",i,".pdf"),
    plot = aa,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 7,
    height = 7,
    units = c("in"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )
}


pathways <- c("TNF","KIT","CXCL","CCL","IL16","IL2","CSF","GRN","GALECTIN","PAR")[c("TNF","KIT","CXCL","CCL","IL16","IL2","CSF","GRN","GALECTIN","PAR") %in% cellchat.sc@netP$pathways]

for (pathway in pathways) {
  p <- netAnalysis_contribution(cellchat.sc, signaling = pathway)
  p <- ggdraw(p)
  ggsave(
    filename = paste0("sc_",pathway,"_contribution_bar.pdf"),
    plot = p,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 7,
    height = 7,
    units = c("in"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )
  pairLR <- extractEnrichedLR(cellchat.sc, signaling = pathway, geneLR.return = FALSE)
  
  
  netVisual(object = cellchat.sc, 
            signaling = pathway, 
            layout = "circle",
            vertex.size.max = 5,
            color.use = whole_color_use_sc,
            weight.scale = T,
            out.format = "pdf")
}


# Compute the network centrality scores
cellchat.sc <- netAnalysis_computeCentrality(cellchat.sc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.sc, 
                                  signaling = "CCL", 
                                  width = 32, 
                                  height = 10, 
                                  color.use = whole_color_use_sc,
                                  font.size = 10)




netAnalysis_signalingRole_scatter(cellchat.sc,color.use = whole_color_use_sc) +
  xlim(c(0,9.1)) +
  ylim(c(0,9.1)) +
  geom_abline(slope = 1) +
  theme(legend.position=c(0.1,0.9))


cellchat.sc <- computeNetSimilarity(cellchat.sc, type = "functional")
cellchat.sc <- netEmbedding(cellchat.sc, type = "functional",umap.method = "uwot")
cellchat.sc <- netClustering(cellchat.sc, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat.sc, type = "functional", label.size = 3.5)+
  theme(legend.position=c(0.9,0.3))


















































WAT.o.sub <- readRDS("WAT_scRNA_MTX_om_20220921.RDS")
WAT.o.sub <- NormalizeData(WAT.o.sub)


cellchat.om <- createCellChat(object = WAT.o.sub,assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)


# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling


# set the used database in the object
cellchat.om@DB <- CellChatDB.use


# subset the expression data of signaling genes for saving computation cost
cellchat.om <- subsetData(cellchat.om) # This step is necessary even if using the whole database


cellchat.om <- identifyOverExpressedGenes(cellchat.om)
cellchat.om <- identifyOverExpressedInteractions(cellchat.om)


cellchat.om <- computeCommunProb(cellchat.om)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.om <- filterCommunication(cellchat.om, min.cells = 10)


df.net <- subsetCommunication(cellchat.om)


cellchat.om <- computeCommunProbPathway(cellchat.om)


cellchat.om <- aggregateNet(cellchat.om)

whole_color_use_om <- c(colorRampPalette(colors = c("#ecd9f1", "#967bce"))(length(grep("lyC|B_cell",unique(cellchat.om@meta$ident),value = T))), #ly purple
                        colorRampPalette(colors = c("#C6DBEF", "#075a84"))(length(grep("myC|Mast",unique(cellchat.om@meta$ident),value = T))), #my bule
                        colorRampPalette(colors = c("#F3E55C", "#E8602D"))(length(grep("ofC",unique(cellchat.om@meta$ident),value = T))), #fb orange
                        "#C0C0C0", # grey
                        colorRampPalette(colors = c("#a6dbbb", "#359566"))(length(grep("vC",unique(cellchat.om@meta$ident),value = T))) ##vc green
)

groupSize.om <- log(log(as.numeric(table(cellchat.om@idents))))
netVisual_circle(cellchat.om@net$weight, 
                 vertex.weight = groupSize.om, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength",
                 arrow.size = 0.05, 
                 vertex.size.max = 5,
                 color.use = whole_color_use_om)




cellchat_df.om <- cellchat.om@net$weight


a <- c(intersect(rownames(cellchat_df.om)[rowSums(cellchat_df.om)==0],colnames(cellchat_df.om)[colSums(cellchat_df.om)==0]))


cellchat_df.om <- cellchat_df.om[!(rownames(cellchat_df.om) %in% a), !(colnames(cellchat_df.om) %in% a)]


annotation_col <- data.frame(row=names(table(cellchat.om@idents)),cell_type = "a",row.names = names(table(cellchat.om@idents)))

annotation_col$cell_type[grep("lyC|B_cell",annotation_col$row)] <- "lyC"
annotation_col$cell_type[grep("myC|Mast",annotation_col$row)] <- "myC"
annotation_col$cell_type[grep("ofC",annotation_col$row)] <- "ofC"
annotation_col$cell_type[grep("adipose",annotation_col$row)] <- "adipose"
annotation_col$cell_type[grep("vC",annotation_col$row)] <- "vC"


annotation_col <- data.frame(cell_type = annotation_col$cell_type,row.names = rownames(annotation_col))

annotation_row <- annotation_col

ann_colors <- list(cell_type = c(lyC = colorRampPalette(colors = c("#ecd9f1", "#967bce"))(3)[2], 
                                 myC = colorRampPalette(colors = c("#C6DBEF", "#075a84"))(3)[2], 
                                 ofC = colorRampPalette(colors = c("#F3E55C", "#E8602D"))(3)[2], 
                                 adipose = "#C0C0C0", 
                                 vC = colorRampPalette(colors = c("#a6dbbb", "#359566"))(3)[2]
)
)

pheatmap(cellchat_df.om,
         cluster_cols = T,
         cluster_rows = T,
         # cellwidth = 8.5,
         # cellheight = 8.5,
         cellwidth = 7.5,
         cellheight = 7.5,
         border_color = NA,
         fontsize_number = 5,
         cutree_cols = 3,
         cutree_rows = 3,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         clustering_method = "ward.D",
         treeheight_col = 20,
         treeheight_row = 20,
         fontsize = 8,
         color = viridis(n = 51, alpha = 1, begin = 0, end = 1, option = "viridis")
)



temp <- cellchat_df.om[rownames(cellchat_df.om) %in% c("ofC06_MSL","ofC11_MSL","ofC03_MSL","adipose","ofC09_FAP","ofC12_MSL","myC09_M2","myC13_M2","myC07_M2" ,"myC10_MMe","ofC14_FAP","myC00_M2","myC01_M2","lyC09_CD16neg_NK","vC02_Blood_EC_artery"), 
                       colnames(cellchat_df.om) %in% c("ofC06_MSL","ofC11_MSL","ofC03_MSL","adipose","ofC09_FAP","ofC12_MSL","myC09_M2","myC13_M2","myC07_M2" ,"myC10_MMe","ofC14_FAP","myC00_M2","myC01_M2","lyC09_CD16neg_NK","vC02_Blood_EC_artery")]


whole_color_use_om <- c(colorRampPalette(colors = c("#ecd9f1", "#967bce"))(length(grep("lyC|B_cell",unique(cellchat.om@meta$ident),value = T))), #ly purple
                        colorRampPalette(colors = c("#C6DBEF", "#075a84"))(length(grep("myC|Mast",unique(cellchat.om@meta$ident),value = T))), #my bule
                        colorRampPalette(colors = c("#F3E55C", "#E8602D"))(length(grep("ofC",unique(cellchat.om@meta$ident),value = T))), #fb orange
                        "#C0C0C0", # grey
                        colorRampPalette(colors = c("#a6dbbb", "#359566"))(length(grep("vC",unique(cellchat.om@meta$ident),value = T))) ##vc green
)

temp_color.om <- c(colorRampPalette(colors = c("#ecd9f1", "#967bce"))(length(grep("lyC|B_cell",unique(cellchat.om@meta$ident),value = T)))[10], #ly purple
                   colorRampPalette(colors = c("#C6DBEF", "#075a84"))(length(grep("myC|Mast",unique(cellchat.om@meta$ident),value = T)))[c(1,2,8,10,11,14)], #my bule
                   colorRampPalette(colors = c("#F3E55C", "#E8602D"))(length(grep("ofC",unique(cellchat.om@meta$ident),value = T)))[c(4,7,10,12,13,15)], #fb orange
                   "#C0C0C0", # grey
                   colorRampPalette(colors = c("#a6dbbb", "#359566"))(length(grep("vC",unique(cellchat.om@meta$ident),value = T)))[3] ##vc green
)

netVisual_circle(temp, 
                 vertex.weight = as.numeric(table(cellchat.om@idents)[names(table(cellchat.om@idents)) %in% c("ofC06_MSL","ofC11_MSL","ofC03_MSL","adipose","ofC09_FAP","ofC12_MSL","myC09_M2","myC13_M2","myC07_M2" ,"myC10_MMe","ofC14_FAP","myC00_M2","myC01_M2","lyC09_CD16neg_NK","vC02_Blood_EC_artery")]), 
                 weight.scale = T, 
                 edge.weight.max = max(temp),title.name = "Weight",
                 vertex.size.max = 5,
                 color.use = temp_color.om
)


pathways <- c("TNF","KIT","CXCL","CCL","IL16","IL2","CSF","GRN","GALECTIN","PAR")[c("TNF","KIT","CXCL","CCL","IL16","IL2","CSF","GRN","GALECTIN","PAR") %in% cellchat.om@netP$pathways]

for (pathway in pathways) {
  p <- netAnalysis_contribution(cellchat.om, signaling = pathway)
  p <- ggdraw(p)
  ggsave(
    filename = paste0("om_",pathway,"_contribution_bar.pdf"),
    plot = p,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 7,
    height = 7,
    units = c("in"),
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )

  netVisual(object = cellchat.om, 
            signaling = pathway, 
            layout = "circle",
            vertex.size.max = 5,
            color.use = whole_color_use_om,
            weight.scale = T,
            out.format = "pdf")
}





for (i in cellchat.om@netP$pathways) {
  a <- netVisual_aggregate(object = cellchat.om, 
                           signaling = i, 
                           layout = "circle",
                           vertex.weight = groupSize.om,
                           vertex.size.max = 5,
                           color.use = whole_color_use_om)
  aa <- ggdraw(a)
  ggsave(
    filename = paste0("om_",i,".pdf"),
    plot = aa,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 7,
    height = 7,
    units = "in",
    dpi = 300,
    limitsize = TRUE,
    bg = NULL
  )
}



# Compute the network centrality scores
cellchat.om <- netAnalysis_computeCentrality(cellchat.om, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways



netAnalysis_signalingRole_scatter(cellchat.om,color.use = whole_color_use_om) +
  xlim(c(0,3.7)) +
  ylim(c(0,3.7)) +
  geom_abline(slope = 1) +
  theme(legend.position=c(0.1,0.9))




cellchat.om <- computeNetSimilarity(cellchat.om, type = "functional")
cellchat.om <- netEmbedding(cellchat.om, type = "functional",umap.method = "uwot")
cellchat.om <- netClustering(cellchat.om, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat.om, type = "functional", label.size = 3.5)+
  theme(legend.position=c(0.9,0.2))
