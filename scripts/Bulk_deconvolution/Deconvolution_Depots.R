#Analysis Immune cells
#setup path variables
DIR_WD <- getwd()
DIR_ROOT <- file.path(getwd(), "..")  # 
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results")
DIR_FIG <- file.path(DIR_RES, "figures")

#load packages
packages <- c("readxl", "tidyverse")
invisible(lapply(packages, library, character.only=TRUE))

#define colors
color_func_blue <- colorRampPalette(colors = c("#C6DBEF", "#075a84"))
color_func_orange <- colorRampPalette(colors = c("#F3E55C", "#E8602D"))
color_func_green <- colorRampPalette(colors = c("#a6dbbb", "#359566"))
color_func_purp <- colorRampPalette(colors = c("#ecd9f1", "#967bce"))

#define function
read_excel_netsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#read xlsx file with marker genes, one cluster per excel sheet
#sheet one contains top markers for each cluster
fap_marker <- read_excel_netsheets(file.path(DIR_RES,"MARKER/FIBRO.xlsx"))
adipo_marker <- read_excel_netsheets(file.path(DIR_RES,"MARKER/ADIPO_res02.xlsx"))
immuno_marker <- read_excel_netsheets(file.path(DIR_RES,"MARKER/IMMUNE_res06.xlsx"))
endo_marker <- read_excel_netsheets(file.path(DIR_RES,"MARKER/ENDO_res06.xlsx"))
myeloid <- read_excel_netsheets(file.path(DIR_RES,"MARKER/MYELOID_res06.xlsx"))
tcell <- read_excel_netsheets(file.path(DIR_RES,"MARKER/Tcells_res065.xlsx"))
#remove first sheet and set to only significantly & positive (upregulated) genes
#need to reduce/clean 
fap_marker <- fap_marker[2:length(fap_marker)]
for (i in 1:length(fap_marker)) {
  fap_marker[[i]] <- fap_marker[[i]][which(fap_marker[[i]]$avg_log2FC > 0 & fap_marker[[i]]$p_val_adj < 0.05),]
  fap_marker[[i]]$cluster <- paste0("FAP_",i-1)
}
adipo_marker <- adipo_marker[2:length(adipo_marker)]
for (i in 1:length(adipo_marker)) {
  adipo_marker[[i]] <- adipo_marker[[i]][which(adipo_marker[[i]]$avg_log2FC > 0 & adipo_marker[[i]]$p_val_adj < 0.05),]
  adipo_marker[[i]]$cluster <- paste0("Adipo_",i-1)
}
immuno_marker <- immuno_marker[2:length(immuno_marker)]
for (i in 1:length(immuno_marker)) {
  immuno_marker[[i]] <- immuno_marker[[i]][which(immuno_marker[[i]]$avg_log2FC > 0 & immuno_marker[[i]]$p_val_adj < 0.05),]
  immuno_marker[[i]]$cluster <- paste0("Immune_",i-1)
}
endo_marker <- endo_marker[2:length(endo_marker)]
for (i in 1:length(endo_marker)) {
  endo_marker[[i]] <- endo_marker[[i]][which(endo_marker[[i]]$avg_log2FC > 0 & endo_marker[[i]]$p_val_adj < 0.05),]
  endo_marker[[i]]$cluster <- paste0("Endo_",i-1)
}
myeloid <- myeloid[2:length(myeloid)]
for (i in 1:length(myeloid)) {
  myeloid[[i]] <- myeloid[[i]][which(myeloid[[i]]$avg_log2FC > 0 & myeloid[[i]]$p_val_adj < 0.05),]
  myeloid[[i]]$cluster <- paste0("Myeloid_",i-1)
}
tcell <- tcell[2:length(tcell)]
for (i in 1:length(tcell)) {
  tcell[[i]] <- tcell[[i]][which(tcell[[i]]$avg_log2FC > 0 & tcell[[i]]$p_val_adj < 0.05),]
  tcell[[i]]$cluster <- paste0("Tcell_",i-1)
}
alle <- bind_rows(c(fap_marker, adipo_marker, immuno_marker, endo_marker, myeloid, tcell))
#subset to necessary column
all_lowres <- alle %>% dplyr::select(cluster, gene, avg_log2FC)

#read and process bulk data
#EMIF = Arner et al.
#EPI = Krieg et al.
#emif expression set
emif <- readRDS(file.path(DIR_DATA,"EMIFeset.rds"))
#emif pheno data
E_pheno <- read.delim(file.path(DIR_DATA,"EMIF_pheno.txt"))
#emif IDs
E_IDs <- read.delim(file.path(DIR_DATA,"EMIF_Subjects.txt"))
E_pheno2 <- merge(E_pheno, E_IDs, by.x = "DNA", by.y = "ID")
row.names(E_pheno2) <- E_pheno2$NAME
E_pheno2 <- E_pheno2[order(rownames(E_pheno2)),]
dir.create(file.path(DIR_RES,"Deconvolution"))
#use BisqueRNA to predict cell type proportions
res <- BisqueRNA::MarkerBasedDecomposition(emif, all_lowres, unique_markers = FALSE, weighted=TRUE, min_gene = 12, w_col = "avg_log2FC")
deconv <- as.data.frame(res$bulk.props)
deconv <- t(deconv)
write.table(deconv, row.names = T, col.names = NA, file = file.path(DIR_RES,"Deconvolution/EMIF_DepotsII.txt")
#deconv <- read.table(file.path(DIR_RES,"Deconvolution/EMIF_DepotsII.txt", header = T , row.names = 1)
deconv <- as.data.frame(deconv)
all(rownames(deconv)== rownames(E_pheno2)) #TRUE
depots_lowres <- cbind(deconv, E_pheno2)
emif_lowres <- depots_lowres
#go through cell type classes to generate violin plots
#FAP
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(23:36)) {
  depots_temp[1:160,1] <- depots_lowres[,i]
  depots_temp$cluster <- colnames(depots_lowres)[[i]]
  depots_temp$tissue <- depots_lowres$TISSUE
  depots_temp$ID <- depots_lowres$DNA
  depots <- rbind(depots, depots_temp)
}
prefix <- "FAP_"
suffix <- seq(0:13)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
depots2 <- depots[depots$cluster == "FAP_0" |depots$cluster == "FAP_3" | depots$cluster == "FAP_4"| depots$cluster == "FAP_7" ,]
depots2$tissue <- factor(depots2$tissue, levels =c("SC","OM"))
ggplot(data = depots2, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))

#Immune
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(37:55)) {
  depots_temp[1:160,1] <- depots_lowres[,i]
  depots_temp$cluster <- colnames(depots_lowres)[[i]]
  depots_temp$tissue <- depots_lowres$TISSUE
  depots_temp$ID <- depots_lowres$DNA
  depots <- rbind(depots, depots_temp)
}
prefix <- "Immune_"
suffix <- seq(0:19)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
ggplot(data = depots, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))
#Endo
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(11:22)) {
  depots_temp[1:160,1] <- depots_lowres[,i]
  depots_temp$cluster <- colnames(depots_lowres)[[i]]
  depots_temp$tissue <- depots_lowres$TISSUE
  depots_temp$ID <- depots_lowres$DNA
  depots <- rbind(depots, depots_temp)
}
prefix <- "Endo_"
suffix <- seq(0:11)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
#For Fig. 3g; Arner et al.
depots2 <- depots[depots$cluster == "Endo_6" | depots$cluster == "Endo_8" ,]
depots2$tissue <- factor(depots2$tissue, levels =c("SC","OM"))
ggplot(data = depots2, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))



#Tcells
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(72:82)) {
  depots_temp[1:160,1] <- emif_lowres[,i]
  depots_temp$cluster <- colnames(emif_lowres)[[i]]
  depots_temp$tissue <- emif_lowres$TISSUE
  depots_temp$ID <- emif_lowres$DNA
  depots <- rbind(depots, depots_temp)
}

prefix <- "Tcell_"
suffix <- seq(0:10)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
#For Fig. 2c; Arner et al.
depots <- depots[depots$cluster == "Tcell_0" | depots$cluster == "Tcell_3",]
depots$tissue <- factor(depots$tissue, levels = c("SC","OM"))
ggplot(data = depots, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))

#myeloid
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(56:71)) {
  depots_temp[1:160,1] <- emif_lowres[,i]
  depots_temp$cluster <- colnames(emif_lowres)[[i]]
  depots_temp$tissue <- emif_lowres$TISSUE
  depots_temp$ID <- emif_lowres$DNA
  depots <- rbind(depots, depots_temp)
}

prefix <- "Myeloid_"
suffix <- seq(0:15)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
#For Fig. 2g; Arner et al.
depots <- depots[depots$cluster == "Myeloid_4" | depots$cluster == "Myeloid_8" | depots$cluster == "Myeloid_7"| depots$cluster == "Myeloid_12",]
depots$tissue <- factor(depots$tissue, levels = c("SC","OM"))
ggplot(data = depots, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))

#Krieg et al.
#read eset
depots<- readRDS(file.path(DIR_DATA,"EPIeset.rds"))
#read phenotypes
depots_pheno <- read.table(file.path(DIR_DATA,"EPI_pheno.csv"), header = T)
rownames(depots_pheno) <- paste0("X",depots_pheno$ID)
#remove unnecessary depots (mesenteric, epiploic)
om_sc <- c("X1950sc","X1950vis","X1959sc","X1959vis","X2029sc","X2029vis","X2039sc","X2039vis","X2040sc","X2040vis","X2045sc","X2045vis","X2059sc",
           "X2059vis","X2060vis","X2065sc","X2065vis","X2087sc","X2087vis","X2118sc","X2118vis","X2119sc","X2119vis","X2120sc","X2125sc","X2125vis","X2127sc",
           "X2127vis","X2128sc","X2128vis","X2181sc","X2181vis","X2229sc","X2229vis","X2301sc","X2301vis","X2333sc","X2333vis","X2364sc","X2364vis","X2459sc",
           "X2459vis","X2464sc","X2464vis","X2470sc","X2470vis")
depots_pheno <- depots_pheno[row.names(depots_pheno) %in% om_sc,]
depots_pheno <- depots_pheno[order(rownames(depots_pheno)),]

#run deconvolution
res <- BisqueRNA::MarkerBasedDecomposition(depots, all_lowres, unique_markers = FALSE, weighted=TRUE, min_gene = 8, w_col = "avg_log2FC")
deconv <- as.data.frame(res$bulk.props)
deconv <- t(deconv)
write.table(deconv, row.names = T, col.names = NA, file = file.path(DIR_RES,"Deconvolution/EPI_all_lowresII.txt")
#deconv <- read.table(file.path(DIR_RES,"Deconvolution/EPI_all_lowresII.txt", header = T , row.names = 1)
deconv <- as.data.frame(deconv)
all(rownames(deconv)== rownames(depots_pheno)) #TRUE
depots_lowres <- cbind(deconv, depots_pheno)
#FAP
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(23:36)) {
  depots_temp[1:46,1] <- depots_lowres[,i]
  depots_temp$cluster <- colnames(depots_lowres)[[i]]
  depots_temp$tissue <- depots_lowres$tissue
  depots_temp$ID <- depots_lowres$DNA
  depots <- rbind(depots, depots_temp)
}
prefix <- "FAP_"
suffix <- seq(0:13)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
#FAPs not in current manuscript version, as split into sc and om to begin with
depots2 <- depots[depots$cluster == "FAP_0" |depots$cluster == "FAP_3" | depots$cluster == "FAP_4"| depots$cluster == "FAP_7" ,]
depots2$tissue <- factor(depots2$tissue, levels =c("sc","om"))
ggplot(data = depots2, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))

#Immune
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(37:55)) {
  depots_temp[1:46,1] <- depots_lowres[,i]
  depots_temp$cluster <- colnames(depots_lowres)[[i]]
  depots_temp$tissue <- depots_lowres$tissue
  depots_temp$ID <- depots_lowres$DNA
  depots <- rbind(depots, depots_temp)
}
prefix <- "Immune_"
suffix <- seq(0:19)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
ggplot(data = depots, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))
#Endo
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(11:22)) {
  depots_temp[1:46,1] <- depots_lowres[,i]
  depots_temp$cluster <- colnames(depots_lowres)[[i]]
  depots_temp$tissue <- depots_lowres$tissue
  depots_temp$ID <- depots_lowres$DNA
  depots <- rbind(depots, depots_temp)
}
prefix <- "Endo_"
suffix <- seq(0:11)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
#Fig 3g; Krieg et al. (lower panel)
depots2 <- depots[depots$cluster == "Endo_6" | depots$cluster == "Endo_8" ,]
depots2$tissue <- factor(depots2$tissue, levels =c("sc","om"))
ggplot(data = depots2, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))

#tcells
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(72:82)) {
  depots_temp[1:46,1] <- depots_lowres[,i]
  depots_temp$cluster <- colnames(depots_lowres)[[i]]
  depots_temp$tissue <- depots_lowres$tissue
  depots_temp$ID <- depots_lowres$DNA
  depots <- rbind(depots, depots_temp)
}

prefix <- "Tcell_"
suffix <- seq(0:10)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
#Fig 2c, Krieg et al. (right panel)
depots <- depots[depots$cluster == "Tcell_0" | depots$cluster == "Tcell_3",]
depots$tissue <- factor(depots$tissue, levels = c("sc","om"))
ggplot(data = depots, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))

#myeloid
depots <- data.frame(
  deconv = integer(),
  cluster = character(),
  tissue= character(), 
  ID = character(),
  stringsAsFactors = F
)
depots_temp <- depots
for (i in c(56:71)) {
  depots_temp[1:46,1] <- depots_lowres[,i]
  depots_temp$cluster <- colnames(depots_lowres)[[i]]
  depots_temp$tissue <- depots_lowres$tissue
  depots_temp$ID <- depots_lowres$DNA
  depots <- rbind(depots, depots_temp)
}

prefix <- "Myeloid_"
suffix <- seq(0:15)-1
depots$cluster <- factor(depots$cluster, levels =paste0(prefix, suffix))
#Fig 2g; Krieg et al. (right panel)
depots <- depots[depots$cluster == "Myeloid_4" | depots$cluster == "Myeloid_8" | depots$cluster == "Myeloid_7"| depots$cluster == "Myeloid_12",]
depots$tissue <- factor(depots$tissue, levels = c("sc","om"))
ggplot(data = depots, aes(x=tissue, y=deconv, fill=tissue))+ geom_violin()+theme_classic()+ facet_wrap(vars(cluster), nrow = 1)+
  scale_fill_manual(values = c("#1E6990",  "#EA6F31"))
