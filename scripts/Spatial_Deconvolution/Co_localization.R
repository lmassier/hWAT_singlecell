# co-localization patterns in STx


options(stringsAsFactors = F)

library(ggplot2)
library(Seurat)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)


## load Seurat object

st_s42 <- readRDS("~/Spatial_Deconvolution/st_s42.RDS")
st_s42$row <- st_s42@images$slice1@coordinates$row
st_s42$col <- st_s42@images$slice1@coordinates$col
st_s42$imagerow <- st_s42@images$slice1@coordinates$imagerow
st_s42$imagecol <- st_s42@images$slice1@coordinates$imagecol

st_s44 <- readRDS("~/Spatial_Deconvolution/st_s44.RDS")
st_s44$row <- st_s44@images$slice1@coordinates$row
st_s44$col <- st_s44@images$slice1@coordinates$col
st_s44$imagerow <- st_s44@images$slice1@coordinates$imagerow
st_s44$imagecol <- st_s44@images$slice1@coordinates$imagecol

st_s46 <- readRDS("~/Spatial_Deconvolution/st_s46.RDS")
st_s46$row <- st_s46@images$slice1@coordinates$row
st_s46$col <- st_s46@images$slice1@coordinates$col
st_s46$imagerow <- st_s46@images$slice1@coordinates$imagerow
st_s46$imagecol <- st_s46@images$slice1@coordinates$imagecol

st_s48 <- readRDS("~/Spatial_Deconvolution/st_s48.RDS")
st_s48$row <- st_s48@images$slice1@coordinates$row
st_s48$col <- st_s48@images$slice1@coordinates$col
st_s48$imagerow <- st_s48@images$slice1@coordinates$imagerow
st_s48$imagecol <- st_s48@images$slice1@coordinates$imagecol

st_s49 <- readRDS("~/Spatial_Deconvolution/st_s49.RDS")
st_s49$row <- st_s49@images$slice1@coordinates$row
st_s49$col <- st_s49@images$slice1@coordinates$col
st_s49$imagerow <- st_s49@images$slice1@coordinates$imagerow
st_s49$imagecol <- st_s49@images$slice1@coordinates$imagecol
st_s49 <- SCTransform(st_s49, assay = "Spatial", verbose = FALSE)
st_s49 <- RunPCA(st_s49, verbose = FALSE)
st_s49 <- RunUMAP(st_s49, dims = 1:30, verbose = FALSE)
st_s49 <- FindNeighbors(st_s49, dims = 1:30, verbose = FALSE)
st_s49 <- FindClusters(st_s49, verbose = FALSE)

st_s50 <- readRDS("~/Spatial_Deconvolution/st_s50.RDS")
st_s50$row <- st_s50@images$slice1@coordinates$row
st_s50$col <- st_s50@images$slice1@coordinates$col
st_s50$imagerow <- st_s50@images$slice1@coordinates$imagerow
st_s50$imagecol <- st_s50@images$slice1@coordinates$imagecol

st_s51 <- readRDS("~/Spatial_Deconvolution/st_s51.RDS")
st_s51$row <- st_s51@images$slice1@coordinates$row
st_s51$col <- st_s51@images$slice1@coordinates$col
st_s51$imagerow <- st_s51@images$slice1@coordinates$imagerow
st_s51$imagecol <- st_s51@images$slice1@coordinates$imagecol

st_s52 <- readRDS("~/Spatial_Deconvolution/st_s52.RDS")
st_s52$row <- st_s52@images$slice1@coordinates$row
st_s52$col <- st_s52@images$slice1@coordinates$col
st_s52$imagerow <- st_s52@images$slice1@coordinates$imagerow
st_s52$imagecol <- st_s52@images$slice1@coordinates$imagecol

st_s54 <- readRDS("~/Spatial_Deconvolution/st_s54.RDS")
st_s54$row <- st_s54@images$slice1@coordinates$row
st_s54$col <- st_s54@images$slice1@coordinates$col
st_s54$imagerow <- st_s54@images$slice1@coordinates$imagerow
st_s54$imagecol <- st_s54@images$slice1@coordinates$imagecol

st_s55 <- readRDS("~/Spatial_Deconvolution/st_s55.RDS")
st_s55$row <- st_s55@images$slice1@coordinates$row
st_s55$col <- st_s55@images$slice1@coordinates$col
st_s55$imagerow <- st_s55@images$slice1@coordinates$imagerow
st_s55$imagecol <- st_s55@images$slice1@coordinates$imagecol



## Input deconvolution results

## cell2location
for (i in c(c(42,44,46,48:52,54,55))) {
  cell2location_temp <- read.csv(paste0("~/zjw/20220528WAT/result/s",i,"_20220811/cell2location_result.csv"),header = 1,row.names = 1)
  cell2location_temp <- cell2location_temp/rowSums(cell2location_temp)
  cell2location_temp$barcode <- rownames(cell2location_temp)
  
  a <-  paste0('cell2location_temp <- merge(cell2location_temp,data.frame(barcode=rownames(st_s',i,'@meta.data)),by="barcode",all=T)')
  eval(parse(text=a))
  print(a)
  
  cell2location_temp[is.na(cell2location_temp)] <- 0
  rownames(cell2location_temp) <- cell2location_temp$barcode
  cell2location_temp <- cell2location_temp[order(cell2location_temp$barcode),-1]
  colnames(cell2location_temp) <-paste0("cell2location_",colnames(cell2location_temp)) %>% gsub(pattern = "q05cell_abundance_w_sf_",replacement = "")
  
  cell2location_temp <- cell2location_temp[,order(colnames(cell2location_temp))]
  a <- paste0("st_s",i,"@meta.data <- cbind(st_s",i,"@meta.data, cell2location_temp)")
  eval(parse(text=a))
  print(a)
}



df <- data.frame()
for (i in c(42,44,46,48:52,54,55)) {
  temp <- get(paste0('st_s',i))
  temp <- temp@meta.data[,-c(1:7)]
  rownames(temp) <- paste0('st_s',i,"_",rownames(temp))
  df <- rbind(df,
              temp)
}

df[df<0.03] <- 0


df_nonzero <- df[,!colSums(df)==0]
cor_mtx <- cor(df_nonzero[,grep("cell2location",colnames(df_nonzero))],method = "pearson")
cor_mtx[!cor_mtx==0] <- 0
cor_mtx[is.na(cor_mtx)] <- 0
for (i in c(42,44,46,48:52,54,55)) {
  temp <- df_nonzero[grep(i,rownames(df_nonzero)),]
  temp <- temp[,grep("cell2location",colnames(temp))]
  temp <- cor(temp,method = "pearson")
  temp[is.na(temp)] <- 0
  cor_mtx <- cor_mtx+temp
}
cor_mtx <- cor_mtx/10
rownames(cor_mtx) <- rownames(cor_mtx) %>% gsub(pattern = "cell2location_",replacement = "")
colnames(cor_mtx) <- colnames(cor_mtx) %>% gsub(pattern = "cell2location_",replacement = "")



annotation_col <- data.frame(row=rownames(cor_mtx),cell_type = "a",row.names = rownames(cor_mtx))

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


a <- pheatmap(cor_mtx,
              cluster_cols = T,
              cluster_rows = T,
              cellwidth = 8.5,
              cellheight = 8.5,
              clustering_method = "ward.D2",
              cutree_cols = 4,
              cutree_rows = 4,
              border_color = NA,
              fontsize_number = 5,
              treeheight_col = 0,
              annotation_col = annotation_col,
              annotation_row = annotation_row,
              annotation_colors = ann_colors,
              fontsize = 8,
              color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51),
              breaks = c(seq(-0.399568829,0,0.399568829/25),seq(0,0.698076379,0.698076379/25)[-1])
)

cor_mtx_temp1 <- cor_mtx[a$tree_row$order[1:16],
                         a$tree_row$order[1:16]]


a <- melt(cor_mtx_temp1)
a <- a[abs(a$value)>(0.1) & abs(a$value)<1,]
a <- a[order(a$Var1),]
a <- a[order(a$value),]
a <- a[seq(1,nrow(a),2),]
a$versus <- paste0(a$Var1,"_vs_",a$Var2)


df_temp <- data.frame()
for (i in c(42,44,46,48:52,54,55)) {
  temp <- df_nonzero[grep(i,rownames(df_nonzero)),]
  temp <- temp[,grep("cell2location",colnames(temp))]
  temp <- cor(temp,method = "pearson")
  temp[is.na(temp)] <- 0
  temp <- melt(temp)
  temp$subject <- paste0("s",i)
  df_temp <- rbind(df_temp,temp)
}

df_temp$Var1 <- df_temp$Var1 %>% gsub(pattern = "cell2location_",replacement = "")
df_temp$Var2 <- df_temp$Var2 %>% gsub(pattern = "cell2location_",replacement = "")

df_temp$versus <- paste0(df_temp$Var1,"_vs_",df_temp$Var2)

df_temp <- dcast(df_temp,formula = versus~subject,value.var = "value")
rownames(df_temp) <- df_temp$versus
df_temp <- df_temp[,-1]

df_temp <- df_temp[rownames(df_temp) %in% a$versus,]

df_temp <- t(df_temp)

pheatmap(df_temp,
         cluster_cols = T,
         cluster_rows = F,
         border_color = NA,
         fontsize_number = 7,
         fontsize = 10,
         color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51),
         breaks = seq(-0.9405891,0.9405891,0.9405891*2/50),
)


## Zoom in 
st_s55_temp_myeloid <- subset(st_s55, subset = imagerow > 4200)
st_s55_temp_myeloid <- subset(st_s55_temp_myeloid, subset = imagerow < 4800)
st_s55_temp_myeloid <- subset(st_s55_temp_myeloid, subset = imagecol > 7000)
st_s55_temp_myeloid <- subset(st_s55_temp_myeloid, subset = imagecol < 7600)

a <- st_s55_temp_myeloid@meta.data
a <- a[,grep("cell2location_myC",colnames(a))]
colnames(a) <- gsub(pattern = "cell2location_",replacement = "",x = colnames(a))
order_c <- names(colSums(a)[order(colSums(a),decreasing = T)])
a <- melt(a)
a$variable <- factor(a$variable,levels = order_c)

aa <- colorRampPalette(colors = c("#C6DBEF", "#075a84"))(17)[1:16]
names(aa) <- unique(a$variable)
aa[levels(a$variable)]

ggboxplot(data = a,x = "variable",y = "value",xlab = "",ylab = "Percentage",fill = "variable",palette = aa[levels(a$variable)],legend = "right",add.params = list(binwidth = 0.02))+
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=90,hjust = 1,vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(10,10),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))

st_s52_temp_myeloid <- subset(st_s52, subset = imagerow > 3700)
st_s52_temp_myeloid <- subset(st_s52_temp_myeloid, subset = imagerow < 4350)
st_s52_temp_myeloid <- subset(st_s52_temp_myeloid, subset = imagecol > 3650)
st_s52_temp_myeloid <- subset(st_s52_temp_myeloid, subset = imagecol < 4300)


a <- st_s52_temp_myeloid@meta.data
a <- a[,grep("cell2location_myC",colnames(a))]
colnames(a) <- gsub(pattern = "cell2location_",replacement = "",x = colnames(a))
order_c <- names(colSums(a)[order(colSums(a),decreasing = T)])
a <- melt(a)
a$variable <- factor(a$variable,levels = order_c)

aa <- colorRampPalette(colors = c("#C6DBEF", "#075a84"))(17)[1:16]
names(aa) <- unique(a$variable)
aa[levels(a$variable)]

ggboxplot(data = a,x = "variable",y = "value",xlab = "",ylab = "Percentage",fill = "variable",palette = aa[levels(a$variable)],legend = "right",add.params = list(binwidth = 0.02))+
  theme(axis.text.x = element_text(face = "plain",size = 11,angle=90,hjust = 1,vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(10,10),
        axis.text.y = element_text(face = "plain",size = 11),
        axis.title.x = element_text(face = "plain",size = 13),
        axis.title.y = element_text(face = "plain",size = 13))



# Blend visualization

st_s49@reductions$umap@cell.embeddings[,1] <- (st_s49$col)
st_s49@reductions$umap@cell.embeddings[,2] <- (-st_s49$row)



FeaturePlot(st_s49, features = c("cell2location_sfC08_FAP", "cell2location_sfC12_FAP"), order = T, min.cutoff = "q1", max.cutoff = "q99", blend = T, cols = c("darkblue", "green", "magenta"), blend.threshold = 0,pt.size = 1.7) &DarkTheme() &NoAxes()
FeaturePlot(st_s49, features = c("cell2location_sfC08_FAP", "cell2location_myC02_LAM"), order = T, min.cutoff = "q1", max.cutoff = "q99", blend = T, cols = c("darkblue", "green", "magenta"), blend.threshold = 0,pt.size = 1.7) &DarkTheme() &NoAxes()
FeaturePlot(st_s49, features = c("cell2location_sfC12_FAP", "cell2location_vC01_Blood_EC_vein"), order = T, min.cutoff = "q1", max.cutoff = "q99", blend = T, cols = c("darkblue", "green", "magenta"), blend.threshold = 0,pt.size = 1.7) &DarkTheme() &NoAxes()
