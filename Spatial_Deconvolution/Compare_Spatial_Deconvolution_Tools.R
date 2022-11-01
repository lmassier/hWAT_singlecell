# Calculate the correlation between different deconvolution tools and different celltypes


options(stringsAsFactors = F)

library(ggplot2)
library(BuenColors)
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
  cell2location_temp <- read.csv(paste0("~/Spatial_Deconvolution/s",i,"/cell2location_result.csv"),header = 1,row.names = 1)
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


## DestVI
for (i in c(c(42,44,46,48:52,54,55))) {
  DestVI_temp <- read.csv(paste0("~/Spatial_Deconvolution/s",i,"/DestVI_result.csv"),row.names = 1)
  DestVI_temp <- DestVI_temp/rowSums(DestVI_temp)
  DestVI_temp$barcode <- rownames(DestVI_temp)
  
  a <-  paste0('DestVI_temp <- merge(DestVI_temp,data.frame(barcode=rownames(st_s',i,'@meta.data)),by="barcode",all=T)')
  eval(parse(text=a))
  print(a)
  
  DestVI_temp[is.na(DestVI_temp)] <- 0
  rownames(DestVI_temp) <- DestVI_temp$barcode
  DestVI_temp <- DestVI_temp[order(DestVI_temp$barcode),-1]
  colnames(DestVI_temp) <-paste0("DestVI_",colnames(DestVI_temp)) %>% gsub(pattern = "X",replacement = "")
  
  DestVI_temp <- DestVI_temp[,order(colnames(DestVI_temp))]
  a <- paste0("st_s",i,"@meta.data <- cbind(st_s",i,"@meta.data, DestVI_temp)")
  eval(parse(text=a))
  print(a)
}


## RCTD
for (i in c(c(42,44,46,48:52,54,55))) {
  RCTD_temp <- read.csv(paste0("~/Spatial_Deconvolution/s",i,"/RCTD_result.csv"),header = 1,row.names = 1)
  RCTD_temp <- RCTD_temp/rowSums(RCTD_temp)
  RCTD_temp$barcode <- rownames(RCTD_temp)
  
  a <-  paste0('RCTD_temp <- merge(RCTD_temp,data.frame(barcode=rownames(st_s',i,'@meta.data)),by="barcode",all=T)')
  eval(parse(text=a))
  print(a)
  
  RCTD_temp[is.na(RCTD_temp)] <- 0
  rownames(RCTD_temp) <- RCTD_temp$barcode
  RCTD_temp <- RCTD_temp[order(RCTD_temp$barcode),-1]
  colnames(RCTD_temp) <-paste0("RCTD_",colnames(RCTD_temp)) %>% gsub(pattern = "X",replacement = "")
  
  RCTD_temp <- RCTD_temp[,order(colnames(RCTD_temp))]
  a <- paste0("st_s",i,"@meta.data <- cbind(st_s",i,"@meta.data, RCTD_temp)")
  eval(parse(text=a))
  print(a)
}


## SPOTlight
for (i in c(c(42,44,46,48:52,54,55))) {
  SPOTlight_temp <- read.csv(paste0("~/Spatial_Deconvolution/s",i,"/SPOTlight_result.csv"),header = 1,row.names = 1)
  a <- paste0('rownames(SPOTlight_temp) <- rownames(st_s',i,'@meta.data)')
  eval(parse(text=a))
  print(a)
  SPOTlight_temp <- SPOTlight_temp/rowSums(SPOTlight_temp)
  SPOTlight_temp$barcode <- rownames(SPOTlight_temp)
  
  a <-  paste0('SPOTlight_temp <- merge(SPOTlight_temp,data.frame(barcode=rownames(st_s',i,'@meta.data)),by="barcode",all=T)')
  eval(parse(text=a))
  print(a)
  
  SPOTlight_temp[is.na(SPOTlight_temp)] <- 0
  rownames(SPOTlight_temp) <- SPOTlight_temp$barcode
  SPOTlight_temp <- SPOTlight_temp[order(SPOTlight_temp$barcode),-1]
  colnames(SPOTlight_temp) <-paste0("SPOTlight_",colnames(SPOTlight_temp)) %>% gsub(pattern = "X",replacement = "")
  
  SPOTlight_temp <- SPOTlight_temp[,order(colnames(SPOTlight_temp))]
  a <- paste0("st_s",i,"@meta.data <- cbind(st_s",i,"@meta.data, SPOTlight_temp)")
  eval(parse(text=a))
  print(a)
}


## stereoscope
for (i in c(c(42,44,46,48:52,54,55))) {
  stereoscope_temp <- read.csv(paste0("~/Spatial_Deconvolution/s",i,"/stereoscope_result.csv"),header = 1,row.names = 1)
  stereoscope_temp <- stereoscope_temp/rowSums(stereoscope_temp)
  stereoscope_temp$barcode <- rownames(stereoscope_temp)
  
  a <-  paste0('stereoscope_temp <- merge(stereoscope_temp,data.frame(barcode=rownames(st_s',i,'@meta.data)),by="barcode",all=T)')
  eval(parse(text=a))
  print(a)
  
  stereoscope_temp[is.na(stereoscope_temp)] <- 0
  rownames(stereoscope_temp) <- stereoscope_temp$barcode
  stereoscope_temp <- stereoscope_temp[order(stereoscope_temp$barcode),-1]
  colnames(stereoscope_temp) <-paste0("stereoscope_",colnames(stereoscope_temp)) %>% gsub(pattern = "X",replacement = "")
  
  stereoscope_temp <- stereoscope_temp[,order(colnames(stereoscope_temp))]
  a <- paste0("st_s",i,"@meta.data <- cbind(st_s",i,"@meta.data, stereoscope_temp)")
  eval(parse(text=a))
  print(a)
}


## Tangram
for (i in c(c(42,44,46,48:52,54,55))) {
  Tangram_temp <- read.csv(paste0("~/Spatial_Deconvolution/s",i,"/Tangram_result.csv"),header = 1,row.names = 1)
  Tangram_temp <- Tangram_temp/rowSums(Tangram_temp)
  Tangram_temp$barcode <- rownames(Tangram_temp)
  
  a <-  paste0('Tangram_temp <- merge(Tangram_temp,data.frame(barcode=rownames(st_s',i,'@meta.data)),by="barcode",all=T)')
  eval(parse(text=a))
  print(a)
  
  Tangram_temp[is.na(Tangram_temp)] <- 0
  rownames(Tangram_temp) <- Tangram_temp$barcode
  Tangram_temp <- Tangram_temp[order(Tangram_temp$barcode),-1]
  colnames(Tangram_temp) <-paste0("Tangram_",colnames(Tangram_temp)) %>% gsub(pattern = "X",replacement = "")
  
  
  Tangram_temp <- Tangram_temp[,order(colnames(Tangram_temp))]
  Tangram_temp <- Tangram_temp/rowSums(Tangram_temp)
  a <- paste0("st_s",i,"@meta.data <- cbind(st_s",i,"@meta.data, Tangram_temp)")
  eval(parse(text=a))
  print(a)
}



## Calculate correlation martix

df <- data.frame()
for (i in c(42,44,46,48:52,54,55)) {
  temp <- get(paste0('st_s',i))
  temp <- temp@meta.data[,-c(1:7)]
  rownames(temp) <- paste0('st_s',i,"_",rownames(temp))
  df <- rbind(df,
              temp)
}

df[df<0.03] <- 0

temp <- strsplit(colnames(df),split = "_")
temp1 <- lapply(temp, function(x) x[1])
temp1 <- unlist(temp1)
temp2 <- lapply(temp, function(x) x[2])
temp2 <- unlist(temp2)

colnames(df) <- paste0(temp2,"_",temp1)

df_temp <- cor(df,method = "pearson")
df_temp[!df_temp==0] <- 0
df_temp[is.na(df_temp)] <- 0
for (i in c(42,44,46,48:52,54,55)) {
  temp <- df[grep(i,rownames(df)),]
  temp <- cor(temp,method = "pearson")
  temp[is.na(temp)] <- 0
  df_temp <- df_temp + temp
}
df_temp <- df_temp/10

pheatmap(df_temp,
         cluster_cols = T,
         cluster_rows = T,
         cellwidth = 9,
         cellheight = 9,
         border_color = NA,
         treeheight_col = 0,
         treeheight_row = 0, 
         angle_col = 270,
         color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51),
         breaks = seq(-0.8897437727,0.8897437727,0.8897437727*2/50),
)



## Example of Myeloid cell

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

SpatialFeaturePlot(st_s52,features = c("cell2location_Myeloid"),pt.size.factor = 1.7,image.alpha = 0,min.cutoff = "q1", max.cutoff = "q99") +
  scale_fill_gradientn(colours = SpatialColors(n = 100)[51:100])

SpatialFeaturePlot(st_s55,features = c("cell2location_Myeloid"),pt.size.factor = 1.7,image.alpha = 0,min.cutoff = "q1", max.cutoff = "q99") +
  scale_fill_gradientn(colours = SpatialColors(n = 100)[51:100])
