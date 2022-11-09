#Analysis FAP II
#setup path variables
DIR_WD <- getwd()
DIR_ROOT <- file.path(getwd(), "..")  # 
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results")
DIR_FIG <- file.path(DIR_RES, "figures")

#define functions
library(readxl)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#load marker genes
all_FAP <- read_excel_allsheets(file.path(DIR_RES,"MARKER/FIBRO.xlsx"))
sc_FAP <- read_excel_allsheets(file.path(DIR_RES,"MARKER/FIBRO_SC06.xlsx"))
om_FAP <- read_excel_allsheets(file.path(DIR_RES,"MARKER/FIBRO_OM05.xlsx"))
pvat_FAP <- read_excel_allsheets(file.path(DIR_RES,"MARKER/FIBRO_pvat05.xlsx"))

#compare marker genes across depot specific clusters
FC <- 0.5
for (i in 1:length(all_FAP)) {
  all_FAP[[i]] <- all_FAP[[i]][which(all_FAP[[i]]$avg_log2FC > FC & all_FAP[[i]]$p_val_adj < 0.05),]
  #all_FAP[[i]] <- all_FAP[[i]] %>% slice_max(avg_log2FC, n=10)
}
for (i in 1:length(sc_FAP)) {
  sc_FAP[[i]] <- sc_FAP[[i]][which(sc_FAP[[i]]$avg_log2FC > FC & sc_FAP[[i]]$p_val_adj < 0.05),]
  #sc_FAP[[i]] <- sc_FAP[[i]] %>% slice_max(avg_log2FC, n=10)
}
for (i in 1:length(om_FAP)) {
  om_FAP[[i]] <- om_FAP[[i]][which(om_FAP[[i]]$avg_log2FC > FC & om_FAP[[i]]$p_val_adj < 0.05),]
  #om_FAP[[i]] <- om_FAP[[i]] %>% slice_max(avg_log2FC, n=10)
}
for (i in 1:length(pvat_FAP)) {
  pvat_FAP[[i]] <- pvat_FAP[[i]][which(pvat_FAP[[i]]$avg_log2FC > FC & pvat_FAP[[i]]$p_val_adj < 0.05),]
  #pvat_FAP[[i]] <- pvat_FAP[[i]] %>% slice_max(avg_log2FC, n=10)
}

#needs to be reduced/cleaned
output <-data.frame(Cohort_1=integer(),
                    Cohort_2=integer(),
                    Number=integer(),
                    Percent_1=integer(),
                    Percent_2=integer(),
                    Total_1= integer(),
                    Total_2= integer(),
                    stringsAsFactors=FALSE)
output2 <- output
#switched to comparisons of interest
for (i in 2:length(all_FAP)) {
  for (j in 2:length(sc_FAP)) {
    overlaps <- intersect(all_FAP[[i]]$gene, sc_FAP[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(all_FAP[[i]]) * 100
    perc2 <- numb/ nrow(sc_FAP[[j]]) * 100
    output2[1,1] <- paste0("ALL_FAP",i-2)
    output2[1,2] <- paste0("SC_FAP",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output2[1,6] <- nrow(all_FAP[[i]])
    output2[1,7] <- nrow(sc_FAP[[j]])
    output <- rbind(output, output2)
  }
  for (j in 2:length(om_FAP)) {
    overlaps <- intersect(all_FAP[[i]]$gene, om_FAP[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(all_FAP[[i]]) * 100
    perc2 <- numb/ nrow(om_FAP[[j]]) * 100
    output2[1,1] <- paste0("ALL_FAP",i-2)
    output2[1,2] <- paste0("OM_FAP",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output2[1,6] <- nrow(all_FAP[[i]])
    output2[1,7] <- nrow(om_FAP[[j]])
    output <- rbind(output, output2)
  }
  for (j in 2:length(pvat_FAP)) {
    overlaps <- intersect(all_FAP[[i]]$gene, pvat_FAP[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(all_FAP[[i]]) * 100
    perc2 <- numb/ nrow(pvat_FAP[[j]]) * 100
    output2[1,1] <- paste0("ALL_FAP",i-2)
    output2[1,2] <- paste0("PVAT_FAP",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output2[1,6] <- nrow(all_FAP[[i]])
    output2[1,7] <- nrow(pvat_FAP[[j]])
    output <- rbind(output, output2)
  }
}

for (i in 2:length(sc_FAP)) {
  for (j in 2:length(om_FAP)) {
    overlaps <- intersect(sc_FAP[[i]]$gene, om_FAP[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sc_FAP[[i]]) * 100
    perc2 <- numb/ nrow(om_FAP[[j]]) * 100
    output2[1,1] <- paste0("SC_FAP",i-2)
    output2[1,2] <- paste0("OM_FAP",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output2[1,6] <- nrow(sc_FAP[[i]])
    output2[1,7] <- nrow(om_FAP[[j]])
    output <- rbind(output, output2)
  }
}  
for (j in 2:length(pvat_FAP)) {
  overlaps <- intersect(sc_FAP[[i]]$gene, pvat_FAP[[j]]$gene)
  numb <- length(overlaps)
  perc <- numb/ nrow(sc_FAP[[i]]) * 100
  perc2 <- numb/ nrow(pvat_FAP[[j]]) * 100
  output2[1,1] <- paste0("SC_FAP",i-2)
  output2[1,2] <- paste0("PVAT_FAP",j-2)
  output2[1,3] <- numb
  output2[1,4] <- perc
  output2[1,5] <- perc2
  output2[1,6] <- nrow(sc_FAP[[i]])
  output2[1,7] <- nrow(pvat_FAP[[j]])
  output <- rbind(output, output2)
}

for (i in 2:length(om_FAP)) {
  
  for (j in 2:length(pvat_FAP)) {
    overlaps <- intersect(om_FAP[[i]]$gene, pvat_FAP[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(om_FAP[[i]]) * 100
    perc2 <- numb/ nrow(pvat_FAP[[j]]) * 100
    output2[1,1] <- paste0("OM_FAP",i-2)
    output2[1,2] <- paste0("PVAT_FAP",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output2[1,6] <- nrow(om_FAP[[i]])
    output2[1,7] <- nrow(pvat_FAP[[j]])
    output <- rbind(output, output2)
  }
}
#calculate jaccard index for comparisons
output$jaccard <- output$Number/(output$Total_1+output$Total_2- output$Number)
output$Cohort_1 <- factor(output$Cohort_1, levels = unique(output$Cohort_1))
output$Cohort_2 <- factor(output$Cohort_2, levels = unique(output$Cohort_2))

library(ggplot2)
library(viridis)

#ggplot(data = output, aes(x=Cohort_1, y=Cohort_2, fill=jaccard))+ geom_tile() +scale_fill_stepsn(colors=viridis(10), n.breaks=9)+
#  theme_classic()+ scale_x_discrete(expand=c(0,0))+ scale_y_discrete(expand=c(0,0))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

x <- split(output, output$Cohort_1)
y <- split(output, output$Cohort_2)

OMvsSC <- data.frame(row.names = y$OM_FAP0_FAP0$Cohort_1)
for (i in 1:15) {
  OMvsSC[1:17,i] <- y[[i]]$jaccard
}
colnames(OMvsSC) <- x$SC_FAP0$Cohort_2
row.names(OMvsSC) <- y$OM_FAP0$Cohort_1
library(pheatmap)
#Figure 4j
pheatmap(OMvsSC, color = viridis(10), border_color = NA)

#check meso and fip marker from mouse studies (ED Fig. 4e)
hep_meso <- c("Lgals2", "Krt18", "Gpm6a", "Upk1b", "Msln", "Igfbp2", "Bnc1", "Tspan8", "Isl1", "Chst4", "Palmd", "Stmn2", "Pkhd1l1", "Lrrn4", "Cxadr", "Lgals7", "Krt8", "Slpi", "Krt19", "Upk3b")
hep_fip <- c("Ly6c1", "Lsp1", "Has2", "Sfrp4", "Stmn4", "Efhd1", "Pi16", "Nov", "Gngt2", "Ptgs2", "Fndc1", "Cmah", "Smpd3", "Anxa3", "Uchl1", "Limch1", "Dact2", "Edn1", "Adamst16", "Chstl1")
library(nichenetr) # used nichenetr to convert mouse to human
hep_meso <- convert_mouse_to_human_symbols(hep_meso)
hep_fip <- convert_mouse_to_human_symbols(hep_fip)
hep_meso <- hep_meso[!is.na(hep_meso)]
hep_fip <- hep_fip[!is.na(hep_fip)] # 4 lost
grad_nik <- colorRampPalette(brewer.pal(n = 7, name = "PRGn"))
DoHeatmap(subset(SE_om, downsample = 100), features = hep_meso, size = 3)+ scale_fill_stepsn(colours = viridis(4))
DoHeatmap(subset(om, downsample = 100), features = hep_fip, size = 3)+ scale_fill_stepsn(colours = viridis(4))
FeaturePlot(SE_om, features = "CD34", pt.size = 2, cols = viridis(15), max.cutoff = 3, order = T)
SE_sc<- FindClusters(object = SE_sc, verbose = T, algorithm = 1, resolution = 0.6)
DoHeatmap(subset(SE_sc, downsample = 100), features = hep_meso, size = 3)
DoHeatmap(subset(SE_sc, downsample = 100), features = hep_fip, size = 3)+ scale_fill_stepsn(colours = viridis(8))
