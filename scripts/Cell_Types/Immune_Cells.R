#Analysis Immune cells
#setup path variables
DIR_WD <- getwd()
DIR_ROOT <- file.path(getwd(), "..")  # 
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results")
DIR_FIG <- file.path(DIR_RES, "figures")

#load packages
packages <- c("Seurat", "patchwork", "harmony", "DT", "dplyr", "writexl", "viridis", "ggplot2")
invisible(lapply(packages, library, character.only=TRUE))

#define colors
color_func_blue <- colorRampPalette(colors = c("#C6DBEF", "#075a84"))
color_func_orange <- colorRampPalette(colors = c("#F3E55C", "#E8602D"))
color_func_green <- colorRampPalette(colors = c("#a6dbbb", "#359566"))
color_func_purp <- colorRampPalette(colors = c("#ecd9f1", "#967bce"))

#load integrated data
immune <- readRDS(file.path(DIR_DATA,"INTEGRATED/Immune.rds"))
immune <- RunUMAP(immune, dims = 1:10, reduction = "scvi", n.components = 2)
immune <- FindNeighbors(object = immune, verbose = T, reduction = "scvi", dims = 1:10)
res <- 0.6
immune<- FindClusters(object = immune, verbose = T, algorithm = 1, resolution = res)

#plot combined UMAP (ED Figure 3a)
blau <- color_func_blue(6) #monocytes, macrophages, DCs
orange <- color_func_orange(2) #T and NK cells
gruen <- color_func_green(2) #B and plasma cells
lila <- color_func_purp(2) #Mast cells

farben <- c(orange[2], blau[1:6],lila[2], gruen[2])
#reduce resolution for plotting general overview to have less clusters
res <- 0.2
immune<- FindClusters(object = immune, verbose = T, algorithm = 1, resolution = res)
p <- DimPlot(object = immune, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = farben,
             label = F, label.box = T,
             label.size = 8, 
             repel = F, 
             pt.size = 1) + 
  NoAxes(); p

#explore marker genes #(ED Figure 3b)
im_genes <- c("MS4A1","CPA3","MRC1", "ITGAX" ,"FABP4","CD3D")
VlnPlot(immune, features = im_genes, pt.size = F,same.y.lims = F, stack = T, combine =T, cols = viridis(6) )+scale_y_discrete(limits=rev)

#calculate proportions & marker genes
res <- 0.6
immune<- FindClusters(object = immune, verbose = T, algorithm = 1, resolution = res)
cluster_prop <- immune[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)
cluster_prop$percent <- cluster_prop$n / sum(cluster_prop$n)*100

cluster_prop$seurat_clusters <- factor(cluster_prop$seurat_clusters, levels = seq(from=0, to=18, by=1))
cluster_prop$n <- NULL
write.table(cluster_prop, file.path(DIR_RES,"Immune_Percent.txt"), sep="\t", row.names = F)
ggplot(data = cluster_prop, aes(x=seurat_clusters, y=1, size=percent))+ geom_count(shape=21, fill=farben)+ theme_classic()

top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(immune, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = logf_thr,
)

se_markers$pct.diff <- se_markers$pct.1 - se_markers$pct.2
se_markers_top <- se_markers %>% 
  group_by(as.numeric(cluster)) %>% 
  dplyr::top_n(n = top_n_genes, wt = avg_log2FC)

se_markers_top <- se_markers_top[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]
datatable(se_markers_top, rownames = F, caption = paste("Top", top_n_genes, "marker genes for the computed Seurat clusters"))
fname <- paste0("MARKER/IMMUNE_res06.xlsx")

cluster_prop <- immune[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

method_spot_count <- aggregate(cluster_prop$n, by=list(method = cluster_prop$method), FUN=sum)
N_snSeq <- method_spot_count$x[method_spot_count$method == "snSeq"]
N_scSeq <- method_spot_count$x[method_spot_count$method == "scSeq"]
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/N_snSeq*100, cluster_prop$n/N_scSeq*100)

write.csv(cluster_prop, file.path(DIR_RES,"clustering_countspots.csv"), row.names = F)
datatable(cluster_prop, rownames = F, caption = paste("Cluster spot count and proportions"))


logf_thr <- 0.1 #set cutoff to speed up marker gene calculation
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)
  markers_seu_clusters[[c]] <- FindMarkers(immune, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.1)
}
#write marker genes to Excel file
sheets <- list()
sheets["top_markers_all_clusters"] <- list(as.data.frame(se_markers_top))
for (i in 1:length(markers_seu_clusters)){
  c_name <- names(markers_seu_clusters[i])
  c_data <- markers_seu_clusters[[i]]
  c_names <- colnames(c_data)
  c_data$gene <- rownames(c_data)
  c_data <- c_data[, c("gene", c_names)]
  sheets[c_name] <- list(c_data)
}

write_xlsx(
  x = sheets,
  path = file.path(DIR_RES, fname),
  col_names = TRUE,
  format_headers = TRUE
)

#proportions by depots
cluster_prop <- immune[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, tissue) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

tissue_spot_count <- aggregate(cluster_prop$n, by=list(depot = cluster_prop$tissue), FUN=sum)
N_om <- tissue_spot_count$x[tissue_spot_count$tissue == "om"]
N_sc <- tissue_spot_count$x[tissue_spot_count$tissue == "sc"]
N_pvat <- tissue_spot_count$x[tissue_spot_count$tissue == "pvat"]
cluster_prop$percent <- ifelse(cluster_prop$tissue == "sc", cluster_prop$n/N_sc*100, ifelse(cluster_prop$tissue == "om", cluster_prop$n/N_om*100, cluster_prop$n/N_pvat*100))


#subset T, NK and NKT cells
#based on res=0.6
T_cells <- subset(x = immune, subset = seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 4 | seurat_clusters == 16)
DimPlot(T_cells)
T_cells <- RunUMAP(T_cells, dims = 1:10, reduction = "scvi", n.components = 2)
T_cells <- FindNeighbors(object = T_cells, verbose = T, reduction = "scvi", dims = 1:10)
#test resolution
for (res in c(0, seq(0.1, 1, 0.1), 1.2) ){
  print(res)
  T_cells <- FindClusters(object = T_cells, verbose = T, algorithm = 1, resolution = res)
}
dims_test <- paste0("SCT_snn_res.", c(0.1, 0.3, 0.5, 0.7, 0.9, 1))
DimPlot(object = T_cells, dims = c(1,2), reduction = "umap", group.by = dims_test, 
        pt.size = 0.5,
        label = T, label.size = 3, ncol = 3) & NoLegend()
#decided on resolution that split between NK and NKT cells after validation
T_cells<- FindClusters(object = T_cells, verbose = T, algorithm = 1, resolution = 0.65)
#reset color variables
blau <- color_func_blue(4)
orange <- color_func_orange(4)
gruen <- color_func_green(2)
lila <- color_func_purp(2)
farben <- c(blau[1], lila[2], orange[1],orange[2], blau[2], orange[4], gruen[1], "lightgrey", blau[3], gruen[2], blau[4])
p <- DimPlot(object = T_cells, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = farben,
             label = F, label.box = T,
             label.size = 8, 
             repel = F, 
             pt.size = 2) + 
  NoAxes(); p
#plot UMAPs for key marker genes
FeaturePlot(T_cells, "IL7R", pt.size = 3, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(T_cells, "IL2RA", pt.size = 2, cols = viridis(5), max.cutoff = 3, order = T)
FeaturePlot(T_cells, "ZBTB16", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)

DimPlot(immune, group.by = "tissue")
saveRDS(T_cells, file.path(DIR_RES,"INTEGRATED/T_Cells.rds"))

#generate T and NK cell marker genes
top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(T_cells, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = logf_thr,
)

se_markers$pct.diff <- se_markers$pct.1 - se_markers$pct.2
se_markers_top <- se_markers %>% 
  group_by(as.numeric(cluster)) %>% 
  dplyr::top_n(n = top_n_genes, wt = avg_log2FC)


cluster_prop <- T_cells[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)
cluster_prop$percent <- cluster_prop$n / sum(cluster_prop$n)*100
cluster_prop$n <- NULL
write.table(cluster_prop, file.path(DIR_RES, "Tcells_Percent.txt"), sep="\t", row.names = F)

se_markers_top <- se_markers_top[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]
datatable(se_markers_top, rownames = F, caption = paste("Top", top_n_genes, "marker genes for the computed Seurat clusters"))
fname <- paste0("MARKER/Tcells_res065.xlsx")

cluster_prop <- T_cells[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

method_spot_count <- aggregate(cluster_prop$n, by=list(method = cluster_prop$method), FUN=sum)
N_snSeq <- method_spot_count$x[method_spot_count$method == "snSeq"]
N_scSeq <- method_spot_count$x[method_spot_count$method == "scSeq"]
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/N_snSeq*100, cluster_prop$n/N_scSeq*100)
write.csv(cluster_prop, file.path(DIR_RES, "clustering_countspots_Tcells.csv"), row.names = F)
datatable(cluster_prop, rownames = F, caption = paste("Cluster spot count and proportions"))
logf_thr <- 0.1
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)

  markers_seu_clusters[[c]] <- FindMarkers(T_cells, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.1)
}

sheets <- list()
sheets["top_markers_all_clusters"] <- list(as.data.frame(se_markers_top))
for (i in 1:length(markers_seu_clusters)){
  c_name <- names(markers_seu_clusters[i])
  c_data <- markers_seu_clusters[[i]]
  c_names <- colnames(c_data)
  c_data$gene <- rownames(c_data)
  c_data <- c_data[, c("gene", c_names)]
  
  sheets[c_name] <- list(c_data)
}
write_xlsx(
  x = sheets,
  path = file.path(DIR_RES, fname),
  col_names = TRUE,
  format_headers = TRUE
)
#analyzing cluster proportion, by depot and method
cluster_prop <- T_cells[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)
cluster_prop$percent <- cluster_prop$n / sum(cluster_prop$n)*100

cluster_prop$seurat_clusters <- factor(cluster_prop$seurat_clusters, levels = seq(from=0, to=10, by=1))
ggplot(data = cluster_prop, aes(x=seurat_clusters, y=1, size=percent))+ geom_count(shape=21, fill=orange)+ theme_classic()

cluster_prop <- T_cells[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, tissue) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

tissue_spot_count <- aggregate(cluster_prop$n, by=list(depot = cluster_prop$tissue), FUN=sum)
N_om <- tissue_spot_count$x[tissue_spot_count$tissue == "om"]
N_sc <- tissue_spot_count$x[tissue_spot_count$tissue == "sc"]
N_pvat <- tissue_spot_count$x[tissue_spot_count$tissue == "pvat"]
cluster_prop$percent <- ifelse(cluster_prop$tissue == "sc", cluster_prop$n/N_sc*100, ifelse(cluster_prop$tissue == "om", cluster_prop$n/N_om*100, cluster_prop$n/N_pvat*100))
#relative proportions between individual clusters (i.e. for circle plots as shown in ED Fig 2e)
cluster_prop1 <- cluster_prop
cluster_prop1 <- NULL
for (i in 0:11) {
  temp <- cluster_prop[which(cluster_prop$seurat_clusters == i),]
  temp$percent2 <- temp$percent/ sum(temp$percent) *100
  cluster_prop1 <- rbind(cluster_prop1, temp)
}
cluster_prop1$group = "depot"
cluster_prop <- T_cells[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)
cluster_prop$percent <- cluster_prop$n / sum(cluster_prop$n)*100
cluster_prop1 <- merge(cluster_prop1, cluster_prop, by="seurat_clusters", all.x = TRUE, suffixes = c("","_cluster"))
cluster_prop1$seurat_clusters <- factor(cluster_prop1$seurat_clusters, levels = seq(from=0, to=10, by=1))
cluster_prop1$tissue <- factor(cluster_prop1$tissue, levels = c("sc", "om", "pvat")) 
#plot proportion by depot as seen in Fig 2b
ggplot(cluster_prop1, aes(y=percent, x=tissue, fill=seurat_clusters)) + geom_bar(position="stack", stat="identity")+ theme_classic()+ 
  scale_fill_manual(values=farben) + scale_x_discrete(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))
write.table(cluster_prop1, file.path(DIR_RES,"Tcells_proportions.txt"), sep="\t", row.names = F, col.names = T)

#recluster myeloid cells
myeloid <- subset(x = immune, subset = seurat_clusters == 0 | seurat_clusters == 3 | seurat_clusters == 6 | seurat_clusters == 9 |
                    seurat_clusters == 10 | seurat_clusters == 13 | seurat_clusters == 5 | seurat_clusters == 7 |
                    seurat_clusters == 11 | seurat_clusters == 14)
myeloid <- RunUMAP(myeloid, dims = 1:10, reduction = "scvi", n.components = 2)
myeloid <- FindNeighbors(object = myeloid, verbose = T, reduction = "scvi", dims = 1:10)
#test different resolutions
for (res in c(0, seq(0.1, 1, 0.1), 1.2) ){
  print(res)
  myeloid <- FindClusters(object = myeloid, verbose = T, algorithm = 1, resolution = res)
}
dims_test <- paste0("SCT_snn_res.", c(0.1, 0.3, 0.5, 0.7, 0.9, 1))
DimPlot(object = myeloid, dims = c(1,2), reduction = "umap", group.by = dims_test, 
        pt.size = 0.5,
        label = T, label.size = 3, ncol = 3) & NoLegend()
myeloid<- FindClusters(object = myeloid, verbose = T, algorithm = 1, resolution = 0.6)
#at res= 0.6 observed key clusters such as LAMs, MMe, DCs that otherwise did not split at lower resolutions

blau <- color_func_blue(7) #M2-like macrophages
orange <- color_func_orange(3) #Monocytes and DC2
gruen <- color_func_green(4) #M1-like macrophages
lila <- color_func_purp(3) #unclear, novel(?), M2-like macrophages as highlighted in the manucript text
farben <- c(blau[1:2], gruen[1], orange[1], blau[3], orange[2],gruen[2],blau[4], lila[2], blau[5],gruen[3], blau[6], lila[3],blau[7], orange[3], gruen[4])
p <- DimPlot(object = myeloid, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = farben,
             label = T, label.box = T,
             label.size = 8, 
             repel = F, 
             pt.size = 2) + 
  NoAxes(); p
FeaturePlot(myeloid, "CD14", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
#Violin plot for Fig. 6c
mox <- c("TNF", "IL1B", "CXCL1","CXCL2","CXCL3","CXCL8","CCL2","CCL3","CCL4")
VlnPlot(myeloid, features = mox, pt.size = F,same.y.lims = F, stack = T, combine =T, cols = viridis(9) )+scale_y_discrete(limits=rev)
#visualy confirming integration of Mox cells from methods/ cohorts
mox_se <- myeloid[,myeloid$seurat_clusters==15]
DimPlot(mox_se)
DimPlot(mox_se, group.by = "orig.ident")


#plot UMAPs highlighting marker gene expression for Fig. 2d 
FeaturePlot(object = myeloid, features = "TIMP1",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "CALD1",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "FOXO1",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "ACACB",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "BCOR",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "CRISPLD2",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "MRC1",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "ITGAX",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "PPARG",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "CD1C",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "CD9",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = myeloid, features = "FCGR3A",pt.size = 2, cols = viridis(10), max.cutoff = 3, order = T)

saveRDS(myeloid, file.path(DIR_RES,"INTEGRATED/Myeloids.rds"))

#Marker genes Myeloid cells
cluster_prop <- myeloid[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, tissue) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)
tissue_spot_count <- aggregate(cluster_prop$n, by=list(depot = cluster_prop$tissue), FUN=sum)
N_om <- tissue_spot_count$x[tissue_spot_count$tissue == "om"]
N_sc <- tissue_spot_count$x[tissue_spot_count$tissue == "sc"]
N_pvat <- tissue_spot_count$x[tissue_spot_count$tissue == "pvat"]
cluster_prop$percent <- ifelse(cluster_prop$tissue == "sc", cluster_prop$n/N_sc*100, ifelse(cluster_prop$tissue == "om", cluster_prop$n/N_om*100, cluster_prop$n/N_pvat*100))
cluster_prop1 <- cluster_prop
cluster_prop1 <- NULL
#analog T cells to compare within cluster between depots
for (i in 0:15) {
  temp <- cluster_prop[which(cluster_prop$seurat_clusters == i),]
  temp$percent2 <- temp$percent/ sum(temp$percent) *100
  cluster_prop1 <- rbind(cluster_prop1, temp)
}
cluster_prop1$group = "depot"
cluster_prop <- myeloid[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)
cluster_prop$percent <- cluster_prop$n / sum(cluster_prop$n)*100
cluster_prop1 <- merge(cluster_prop1, cluster_prop, by="seurat_clusters", all.x = TRUE, suffixes = c("","_cluster"))
cluster_prop1$tissue <- factor(cluster_prop1$tissue, levels = c("sc", "om", "pvat")) 
#plot bar chart Fig 2f
ggplot(cluster_prop1, aes(y=percent, x=(tissue), fill=seurat_clusters)) + geom_bar(position="stack", stat="identity", color="black")+ theme_classic()+ 
  scale_fill_manual(values=farben) + scale_x_discrete(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))
write.table(cluster_prop1, file.path(DIR_RES,"Myeloid_proportions.txt"), sep="\t", row.names = F, col.names = T)

#Marker genes myeloid cells
top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(myeloid, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = logf_thr,
)

se_markers$pct.diff <- se_markers$pct.1 - se_markers$pct.2
se_markers_top <- se_markers %>% 
  group_by(as.numeric(cluster)) %>% 
  dplyr::top_n(n = top_n_genes, wt = avg_log2FC)

cluster_prop <- myeloid[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)
cluster_prop$percent <- cluster_prop$n / sum(cluster_prop$n)*100
cluster_prop$n <- NULL
write.table(cluster_prop, file.path(DIR_RES,"Myeloid_Percent.txt"), sep="\t", row.names = F)

se_markers_top <- se_markers_top[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]
datatable(se_markers_top, rownames = F, caption = paste("Top", top_n_genes, "marker genes for the computed Seurat clusters"))
fname <- paste0("MARKER/MYELOID_res6.xlsx")

cluster_prop <- myeloid[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)


method_spot_count <- aggregate(cluster_prop$n, by=list(method = cluster_prop$method), FUN=sum)
N_snSeq <- method_spot_count$x[method_spot_count$method == "snSeq"]
N_scSeq <- method_spot_count$x[method_spot_count$method == "scSeq"]
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/N_snSeq*100, cluster_prop$n/N_scSeq*100)
write.csv(cluster_prop, file.path(DIR_RES,"clustering_countspots_myeloid.csv"), row.names = F)
datatable(cluster_prop, rownames = F, caption = paste("Cluster spot count and proportions"))

logf_thr <- 0.1
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)
  
  markers_seu_clusters[[c]] <- FindMarkers(myeloid, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.1)
}

sheets <- list()
sheets["top_markers_all_clusters"] <- list(as.data.frame(se_markers_top))
for (i in 1:length(markers_seu_clusters)){
  c_name <- names(markers_seu_clusters[i])
  c_data <- markers_seu_clusters[[i]]
  c_names <- colnames(c_data)
  c_data$gene <- rownames(c_data)
  c_data <- c_data[, c("gene", c_names)]
  
  sheets[c_name] <- list(c_data)
}

DIR_RES <- getwd()

write_xlsx(
  x = sheets,
  path = file.path(DIR_RES, fname),
  col_names = TRUE,
  format_headers = TRUE
)

#reclustering and analysis of Mast cells and B cells as shown in ED Fig 3 e and f
#Mast cells
mast <- subset(x = immune, subset = seurat_clusters == 12)
mast <- RunUMAP(mast, dims = 1:10, reduction = "scvi", n.components = 2)
mast <- FindNeighbors(object = mast, verbose = T, reduction = "scvi", dims = 1:10)
for (res in c(0, seq(0.1, 1, 0.1), 1.2) ){
  print(res)
  mast <- FindClusters(object = mast, verbose = T, algorithm = 1, resolution = res)
}
dims_test <- paste0("SCT_snn_res.", c(0.1, 0.3, 0.5, 0.7, 0.9, 1))
DimPlot(object = mast, dims = c(1,2), reduction = "umap", group.by = dims_test, 
        pt.size = 0.5,
        label = T, label.size = 3, ncol = 3) & NoLegend()
mast<- FindClusters(object = mast, verbose = T, algorithm = 1, resolution = 0.6)

#key marker highlighted in UMAPs ED Fig 3f
FeaturePlot(object = mast, features = "CMA1",pt.size = 3, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = mast, features = "KIT",pt.size = 3, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = mast, features = "SIGLEC8",pt.size = 3, cols = viridis(10), max.cutoff = 3, order = T)

# no clear subclusters, change everything to same purple color
p <- DimPlot(object = mast, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = c(rep(color_func_purp(2)[2],4)),
             label = F, label.box = T,
             label.size = 8, 
             repel = F, 
             pt.size = 3) + 
  NoAxes(); p
saveRDS(mast, file.path(DIR_RES,"INTEGRATED/mast_cells.rds"))

#Marker genes for mast cells
top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(mast, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = logf_thr,
)

se_markers$pct.diff <- se_markers$pct.1 - se_markers$pct.2
se_markers_top <- se_markers %>% 
  group_by(as.numeric(cluster)) %>% 
  dplyr::top_n(n = top_n_genes, wt = avg_log2FC)

se_markers_top <- se_markers_top[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]
datatable(se_markers_top, rownames = F, caption = paste("Top", top_n_genes, "marker genes for the computed Seurat clusters"))
fname <- paste0("MARKER/Mast_cells_res06.xlsx")

cluster_prop <- mast[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

method_spot_count <- aggregate(cluster_prop$n, by=list(subject = cluster_prop$method), FUN=sum)
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/14595*100, cluster_prop$n/16971*100)

write.csv(cluster_prop, file.path(DIR_RES, "clustering_countspots_mastcells.csv"), row.names = F)
datatable(cluster_prop, rownames = F, caption = paste("Cluster spot count and proportions"))


logf_thr <- 0.1
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)
  
  markers_seu_clusters[[c]] <- FindMarkers(mast, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.1)
}

sheets <- list()
sheets["top_markers_all_clusters"] <- list(as.data.frame(se_markers_top))
for (i in 1:length(markers_seu_clusters)){
  c_name <- names(markers_seu_clusters[i])
  c_data <- markers_seu_clusters[[i]]
  c_names <- colnames(c_data)
  c_data$gene <- rownames(c_data)
  c_data <- c_data[, c("gene", c_names)]
  
  sheets[c_name] <- list(c_data)
}

write_xlsx(
  x = sheets,
  path = file.path(DIR_RES, fname),
  col_names = TRUE,
  format_headers = TRUE
)




#subset B
B_cells <- subset(x = immune, subset = seurat_clusters == 15)
B_cells <- RunUMAP(B_cells, dims = 1:10, reduction = "scvi", n.components = 2)
B_cells <- FindNeighbors(object = B_cells, verbose = T, reduction = "scvi", dims = 1:10)
for (res in c(0, seq(0.1, 1, 0.1), 1.2) ){
  print(res)
  B_cells <- FindClusters(object = B_cells, verbose = T, algorithm = 1, resolution = res)
}
dims_test <- paste0("SCT_snn_res.", c(0.1, 0.3, 0.5, 0.7, 0.9, 1))
DimPlot(object = B_cells, dims = c(1,2), reduction = "umap", group.by = dims_test, 
        pt.size = 0.5,
        label = T, label.size = 3, ncol = 3) & NoLegend()
B_cells<- FindClusters(object = B_cells, verbose = T, algorithm = 1, resolution = 0.5)
#split into 2 distinct populations
p <- DimPlot(object = B_cells, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = c(color_func_blue(2)[2],color_func_orange(2)[2]),
             label = F, label.box = T,
             label.size = 8, 
             repel = F, 
             pt.size = 4) + 
  NoAxes(); p
#Expression UMAPs as shown in ED Fig 3e
FeaturePlot(object = B_cells, features = "HLA-DRA",pt.size = 3, cols = viridis(10), max.cutoff = 3, order = T)
FeaturePlot(object = B_cells, features = "LYN",pt.size = 3, cols = viridis(10), max.cutoff = 3, order = T)

saveRDS(B_cells, file.path(DIR_RES,"INTEGRATED/B_cells.rds"))

#Marker genes B cells
top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(B_cells, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = logf_thr,
)

se_markers$pct.diff <- se_markers$pct.1 - se_markers$pct.2
se_markers_top <- se_markers %>% 
  group_by(as.numeric(cluster)) %>% 
  dplyr::top_n(n = top_n_genes, wt = avg_log2FC)

se_markers_top <- se_markers_top[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]
datatable(se_markers_top, rownames = F, caption = paste("Top", top_n_genes, "marker genes for the computed Seurat clusters"))
fname <- paste0("MARKER/Bcells_res06.xlsx")

cluster_prop <- B_cells[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

method_spot_count <- aggregate(cluster_prop$n, by=list(subject = cluster_prop$method), FUN=sum)
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/14595*100, cluster_prop$n/16971*100)

write.csv(cluster_prop, file.path(DIR_RES,"clustering_countspots_Bcells.csv"), row.names = F)
datatable(cluster_prop, rownames = F, caption = paste("Cluster spot count and proportions"))


logf_thr <- 0.1
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)
  
  markers_seu_clusters[[c]] <- FindMarkers(B_cells, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.1)
}

sheets <- list()
sheets["top_markers_all_clusters"] <- list(as.data.frame(se_markers_top))
for (i in 1:length(markers_seu_clusters)){
  c_name <- names(markers_seu_clusters[i])
  c_data <- markers_seu_clusters[[i]]
  c_names <- colnames(c_data)
  c_data$gene <- rownames(c_data)
  c_data <- c_data[, c("gene", c_names)]
  
  sheets[c_name] <- list(c_data)
}

write_xlsx(
  x = sheets,
  path = file.path(DIR_RES, fname),
  col_names = TRUE,
  format_headers = TRUE
)