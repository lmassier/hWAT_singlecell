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
fibroblasts <- readRDS(file.path(DIR_DATA,"INTEGRATED/Fibroblasts.rds"))
fibroblasts <- RunUMAP(fibroblasts, dims = 1:10, reduction = "scvi", n.components = 2)
fibroblasts <- FindNeighbors(object = fibroblasts, verbose = T, reduction = "scvi", dims = 1:10)
fibroblasts<- FindClusters(object = fibroblasts, verbose = T, algorithm = 1, resolution = 0.4)
saveRDS(fibroblasts, "NEW_INTEGRATED/fibroblasts_8000.rds")

#integration not good between sc and om, due to MSL cells and other differences
#split to Sc and Om for downstream analysis
#omental FAPs and MSL
om <- fibroblasts[,fibroblasts$tissue == "om"]
om <- RunUMAP(om, dims = 1:10, reduction = "scvi", n.components = 2)
om <- FindNeighbors(object = om, verbose = T, reduction = "scvi", dims = 1:10)
om<- FindClusters(object = om, verbose = T, algorithm = 1, resolution = 0.5)
DimPlot(om, label = T)

saveRDS(om, file.path(DIR_DATA,"INTEGRATED/fibroblasts_om.rds"))

#UMAP Fig 4i
p <- DimPlot(object = om, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = c(color_func_blue(4),color_func_orange(4), color_func_green(4),color_func_purp(3)),
             label = F, label.box = T,
             label.size = 8, 
             repel = F, 
             pt.size = 2) + 
  NoAxes(); p

#check distribution of known marker genes
#ED Fig. 4d
library(viridis)
FeaturePlot(om, "APOD", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(om, "CD55", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(om, "EZR", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(om, "CD74", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(om, "MT2A", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(om, "DPP4", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)

#calculate marker genes
top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(om, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = logf_thr,
)

se_markers$pct.diff <- se_markers$pct.1 - se_markers$pct.2
se_markers_top <- se_markers %>% 
  group_by(as.numeric(cluster)) %>% 
  dplyr::top_n(n = top_n_genes, wt = avg_log2FC)

se_markers_top <- se_markers_top[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]
fname <- paste0("MARKER/FIBRO_OM05.xlsx")

cluster_prop <- om[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

#visualize method distribution
method_spot_count <- aggregate(cluster_prop$n, by=list(method = cluster_prop$method), FUN=sum)
N_snSeq <- method_spot_count$x[method_spot_count$method == "snSeq"]
N_scSeq <- method_spot_count$x[method_spot_count$method == "scSeq"]
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/N_snSeq*100, cluster_prop$n/N_scSeq*100)

cluster_prop <- om[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

method_spot_count <- aggregate(cluster_prop$n, by=list(subject = cluster_prop$method), FUN=sum)
cluster_prop$percent <- cluster_prop$n/sum(cluster_prop$n)
ggplot(cluster_prop, aes(x=1, y=percent, fill=seurat_clusters))+ geom_bar(stat = "identity", position = "stack")+
  theme_classic() + scale_y_continuous(expand=c(0,0))+ scale_fill_manual(values = c(color_func_blue(4),color_func_orange(4), color_func_green(4),color_func_purp(3)))


logf_thr <- 0.2
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)
  
  markers_seu_clusters[[c]] <- FindMarkers(om, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.2)
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
library(writexl)
write_xlsx(
  x = sheets,
  path = file.path(DIR_RES, fname),
  col_names = TRUE,
  format_headers = TRUE
)

#perivascular WAT, only one study!
pvat <- fibroblasts[,fibroblasts$tissue == "pvat"]
pvat <- RunUMAP(pvat, dims = 1:10, reduction = "scvi", n.cpvatponents = 2)
pvat <- FindNeighbors(object = pvat, verbose = T, reduction = "scvi", dims = 1:10)
pvat<- FindClusters(object = pvat, verbose = T, algorithm = 1, resolution = 0.5)
saveRDS(pvat, file.path(DIR_DATA,"INTEGRATED/fibroblasts_pvat.rds"))
#ED Fig 4f upper panel
p <- DimPlot(object = pvat, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = c(color_func_blue(2),color_func_orange(2), color_func_green(2),color_func_purp(2)),
             label = F, label.box = T,
             label.size = 8, 
             repel = F, 
             pt.size = 2) + 
  NoAxes(); p
#ED Fig 4f, lower panel
FeaturePlot(pvat, "PDGFRA", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(pvat, "PDGFRB", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(pvat, "CD34", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(pvat, "CD74", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(pvat, "MT2A", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(pvat, "FAM155A", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)

#Marker genes pvat
cluster_prop <- pvat[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

cluster_prop$percent <- cluster_prop$n/sum(cluster_prop$n)
ggplot(cluster_prop, aes(x=1, y=percent, fill=seurat_clusters))+ geom_bar(stat = "identity", position = "stack")+
  theme_classic() + scale_y_continuous(expand=c(0,0))+ scale_fill_manual(values = c(color_func_blue(2),color_func_orange(2), color_func_green(2),color_func_purp(2)))

top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(pvat, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = logf_thr,
)

se_markers$pct.diff <- se_markers$pct.1 - se_markers$pct.2
se_markers_top <- se_markers %>% 
  group_by(as.numeric(cluster)) %>% 
  dplyr::top_n(n = top_n_genes, wt = avg_log2FC)

se_markers_top <- se_markers_top[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]
fname <- paste0("MARKER/FIBRO_pvat05.xlsx")

cluster_prop <- pvat[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

method_spot_count <- aggregate(cluster_prop$n, by=list(method = cluster_prop$method), FUN=sum)
N_snSeq <- method_spot_count$x[method_spot_count$method == "snSeq"]
N_scSeq <- method_spot_count$x[method_spot_count$method == "scSeq"]
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/N_snSeq*100, cluster_prop$n/N_scSeq*100)

logf_thr <- 0.1
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)
  
  markers_seu_clusters[[c]] <- FindMarkers(pvat, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.1)
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
library(writexl)
write_xlsx(
  x = sheets,
  path = file.path(DIR_RES, fname),
  col_names = TRUE,
  format_headers = TRUE
)

#subcutaneous
subc <- subset(fibroblasts, subset = tissue == "sc")

subc <- RunUMAP(subc, dims = 1:10, reduction = "scvi", n.components = 2)
subc <- FindNeighbors(object = subc, verbose = T, reduction = "scvi", dims = 1:10)
subc<- FindClusters(object = subc, verbose = T, algorithm = 1, resolution = 0.6)

#define cluster specific colors for sc WAT FAPs based on knowledge from adipogenesis data (stem cell, transient, committed, not regulated)
blau <- color_func_blue(3)
orange <- color_func_orange(2)
gruen <- color_func_green(9)
lila <- color_func_purp(5)
farben <- c(blau[2],gruen[5], orange[1],gruen[1],lila[2:3],gruen[2],lila[4:5], gruen[3:4],gruen[6],gruen[7], blau[1],gruen[8:9], blau[3])
#UMAP Fig 4a
p <- DimPlot(object = subc, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = farben,
             label = F, label.box = T,
             label.size = 8, 
             repel = F, 
             pt.size = 2, order = T) + 
  NoAxes(); p

#Figure 4b, ED Fig 4a
FeaturePlot(subc, features = c("PDGFRA","PDGFRB","CD34","LY6"), pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(subc, "CD74", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
FeaturePlot(subc, "EZR", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)

#Running monocle3 for sc FAPs based on knowledge that CD55/PI16 cluster represent stem cells across organs and species
#following http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html
set.seed(1234)

subc.cds <- as.cell_data_set(subc)
subc.cds <- cluster_cells(cds = subc.cds)
subc.cds <- learn_graph(subc.cds, use_partition = TRUE)
subc.cds <- order_cells(subc.cds, reduction_method = "UMAP")

plot_cells(
  cds = subc.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, trajectory_graph_color = "black",
  trajectory_graph_segment_size = 1.5, label_branch_points = F, label_leaves = F, label_roots = F
)
subc <- AddMetaData(
  object = subc,
  metadata = subc.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)
#Fig 4d
FeaturePlot(subc, "Pseudotime", pt.size = 0.1) & scale_color_viridis_c()

#Find Markers sc FAPs
top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(subc, 
                             only.pos = TRUE, 
                             min.pct = 0.1, 
                             logfc.threshold = logf_thr,
)

se_markers$pct.diff <- se_markers$pct.1 - se_markers$pct.2
se_markers_top <- se_markers %>% 
  group_by(as.numeric(cluster)) %>% 
  dplyr::top_n(n = top_n_genes, wt = avg_log2FC)


se_markers_top <- se_markers_top[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]

fname <- paste0("MARKER/FIBRO_SC06.xlsx")

cluster_prop <- subc[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)


method_spot_count <- aggregate(cluster_prop$n, by=list(method = cluster_prop$method), FUN=sum)
N_snSeq <- method_spot_count$x[method_spot_count$method == "snSeq"]
N_scSeq <- method_spot_count$x[method_spot_count$method == "scSeq"]
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/N_snSeq*100, cluster_prop$n/N_scSeq*100)

logf_thr <- 0.1
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)
  
  markers_seu_clusters[[c]] <- FindMarkers(subc, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.1)
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
library(writexl)
write_xlsx(
  x = sheets,
  path = file.path(DIR_RES, fname),
  col_names = TRUE,
  format_headers = TRUE
)
