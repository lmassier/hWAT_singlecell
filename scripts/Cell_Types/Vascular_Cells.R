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
vascular <- readRDS(file.path(DIR_DATA,"MERGED/Vascular.rds"))
vascular <- RunUMAP(vascular, dims = 1:10, reduction = "scvi", n.components = 2)
vascular <- FindNeighbors(object = vascular, verbose = T, reduction = "scvi", dims = 1:10)
vascular<- FindClusters(object = vascular, verbose = T, algorithm = 1, resolution = 0.6)
saveRDS(vascular, file.path(DIR_RES,"INTEGRATED/Vascular.rds"))

#visualize UMAP by method/cohort/depot
DimPlot(vascular, group.by = c("method", "orig.ident","tissue"), shuffle = T)
#assign colors
blau <- color_func_blue(5)
lila <- color_func_purp(3)
gruen <- color_func_green(2)
reds <- color_func_orange(2)
farben <- c(blau[1:4], gruen[1],reds[1:2],gruen[2], lila, blau[5]) 
p <- DimPlot(object = vascular, 
             reduction = "umap", dims = c(1,2), 
             group.by = "seurat_clusters",
             cols = farben,
             label = T, label.box = T, 
             label.size = 8, 
             repel = F, 
             pt.size = 2) + 
  NoAxes(); p

#top marker genes per cluster for Violin Plot Fig 3b
vascular_gene <- c("PCDH10" ,"CD48", "TYROBP" ,"TIMP1" ,"ACTA2", "LYVE1", "APOD" ,"NCKAP5", "LDLR", "EFNB2", "ACKR1", "CA4")
VlnPlot(vascular, features = vascular_gene, pt.size = F,same.y.lims = F, stack = T, combine =T, cols = viridis(12) )+scale_y_discrete(limits=rev)

#Marker genes and proportion
cluster_prop <- vascular[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)
cluster_prop$percent <- cluster_prop$n / sum(cluster_prop$n)*100
cluster_prop$n <- NULL
write.table(cluster_prop, file.path(DIR_RES,"Vascular_Percent.txt"), sep="\t", row.names = F)
cluster_prop$seurat_clusters <- factor(cluster_prop$seurat_clusters, levels = seq(from=0, to=11, by=1))
ggplot(data = cluster_prop, aes(x=seurat_clusters, y=1, size=percent))+ geom_count(shape=21)+ theme_classic()


cluster_prop <- vascular[[]]
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
for (i in 0:11) {
  temp <- cluster_prop[which(cluster_prop$seurat_clusters == i),]
  temp$percent2 <- temp$percent/ sum(temp$percent) *100
  cluster_prop1 <- rbind(cluster_prop1, temp)
}
cluster_prop1$group = "depot"
cluster_prop1$tissue <- factor(cluster_prop1$tissue, levels = c("sc", "om", "pvat")) 
#car chart Fig 3f
ggplot(cluster_prop1, aes(y=percent, x=tissue, fill=seurat_clusters)) + geom_bar(position="stack", stat="identity", color="black")+ 
  theme_classic() + scale_x_discrete(expand = c(0,0))+ scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = farben)


cluster_prop <- vascular[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

method_spot_count <- aggregate(cluster_prop$n, by=list(method = cluster_prop$method), FUN=sum)
N_snSeq <- method_spot_count$x[method_spot_count$method == "snSeq"]
N_scSeq <- method_spot_count$x[method_spot_count$method == "scSeq"]
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/N_snSeq*100, cluster_prop$n/N_scSeq*100)
cluster_prop2 <- cluster_prop
cluster_prop2 <- NULL
for (i in 0:11) {
  temp <- cluster_prop[which(cluster_prop$seurat_clusters == i),]
  temp$percent2 <- temp$percent/ sum(temp$percent) *100
  cluster_prop2 <- rbind(cluster_prop2, temp)
}
cluster_prop2$group = "method"

cluster_prop2$seurat_clusters <- factor(cluster_prop2$seurat_clusters)
#ED Fig 2e
ggplot(cluster_prop2, aes(x="", y=percent2, fill=method,percent_cluster2 )) +
  geom_bar(stat="identity", color="white") +
  coord_polar("y", start=0) + facet_wrap(vars(seurat_clusters), nrow = 1)+theme_void()+ scale_fill_manual(values = c("#85448e", "#398753"))

#visualize if clusters are found in each cohort (not in current manuscript)
cluster_prop <- vascular[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, orig.ident) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

cluster_prop$percent <- ""
cluster_prop2 <- cluster_prop[NULL,]

method_spot_count <- aggregate(cluster_prop$n, by=list(subject = cluster_prop$orig.ident), FUN=sum)
for (i in unique(method_spot_count$subject)) {
  temp <- cluster_prop[cluster_prop$orig.ident == i,]
  temp$percent <- temp$n/(method_spot_count[method_spot_count$subject ==i,]$x)*100
  cluster_prop2 <- rbind(cluster_prop2, temp)
}
ggplot(data = cluster_prop2, aes(x=orig.ident, y = percent, fill= orig.ident)) + geom_bar(stat = "identity")+  
  facet_wrap( vars(seurat_clusters), nrow = 1, ncol = 15, scales = "free")+ 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  
  scale_y_continuous(expand = c(0,0)) 

#FindMarker genes
top_n_genes <- 20
logf_thr <- 0.15

se_markers <- FindAllMarkers(vascular, 
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
fname <- paste0("MARKER/Vascular_res06.xlsx")

cluster_prop <- vascular[[]]
cluster_prop <- cluster_prop %>% 
  group_by(seurat_clusters, method) %>%
  dplyr::count()
cluster_prop <- as.data.frame(cluster_prop)

method_spot_count <- aggregate(cluster_prop$n, by=list(method = cluster_prop$method), FUN=sum)
N_snSeq <- method_spot_count$x[method_spot_count$method == "snSeq"]
N_scSeq <- method_spot_count$x[method_spot_count$method == "scSeq"]
cluster_prop$percent <- ifelse(cluster_prop$method == "snSeq", cluster_prop$n/N_snSeq*100, cluster_prop$n/N_scSeq*100)

#write.csv(cluster_prop, file.path(DIR_RES,"clustering_countspots_vascular.csv", row.names = F)
datatable(cluster_prop, rownames = F, caption = paste("Cluster spot count and proportions"))

logf_thr <- 0.1
cluster_calc_de <- unique(sort(as.numeric(as.character(subset(cluster_prop, n>3)$seurat_clusters))))
markers_seu_clusters <- list()
for( i in cluster_calc_de ) {
  print(paste("Finding markers for cluster", i))
  c <- paste0("cluster_", i)
  
  markers_seu_clusters[[c]] <- FindMarkers(vascular, ident.1 = as.character(i), logfc.threshold = logf_thr, min.pct = 0.1)
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

