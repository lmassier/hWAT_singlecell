#redo network v3
#prepare network data
# one node for each cluster in each data set
# edges if overlapping genes, weight = amount of genes
# start by reading data and generating nodes
library(readxl)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

DIR_WD <- getwd()
DIR_ROOT <- file.path(getwd(), "..")  # 
DIR_DATA <- file.path(DIR_ROOT, "data")
DIR_RES <- file.path(DIR_ROOT, "results")
DIR_FIG <- file.path(DIR_RES, "figures")

#load scSeq marker
setwd(file.path(DIR_WD,"scSEQ_SVF/"))
vijay_sc <- read_excel_allsheets("./grundberg/grundberg_sc_marker.xlsx") 
vijay_om <- read_excel_allsheets("./grundberg/grundberg_om_marker.xlsx")
jaitin <- read_excel_allsheets("./Jaitin/jaitin_marker.xlsx")
merrick <- read_excel_allsheets("./merrick/merrick_marker.xlsx")
karunakaran <- read_excel_allsheets("./Karunakaran/karuna_marker.xlsx")
acosta <- read_excel_allsheets("./acosta/acosta_marker.xlsx")
hildreth <- read_excel_allsheets("./Hildreth/hildreth_marker.xlsx")
#load spatial marker
spatial <- read_excel_allsheets("./../spatial/baseline.markers_ALL.xlsx")
#load snSeq marker
setwd(file.path(DIR_WD,"snSeq_AT/"))
sun <- read_excel_allsheets("./Sun_Wolfrum/sun_marker.xlsx")
angueira <- read_excel_allsheets("./Angueira_Seale/seale_marker.xlsx")
lipidlab <- read_excel_allsheets("./lipidlab/lipidlab_marker.xlsx") #Massier et al. #1
leipzig_sc <- read_excel_allsheets("./blueher/leipzig_sc_marker.xlsx") #Massier et al. #2
leipzig_om <- read_excel_allsheets("./blueher/leipzig_om_marker.xlsx") #Massier et al. #4
rosen_sn_om <- read_excel_allsheets("./rosen/sn_omental_marker.xlsx")
rosen_sn_sc <- read_excel_allsheets("./rosen/sn_other_marker.xlsx")
rosen_sc <- read_excel_allsheets("./rosen/sc_rosen_marker.xlsx") 
lipid_1s <- read_excel_allsheets("./lipid_sn_sc/marker.xlsx") #Massier et al. #3
data_list <- list(vijay_om, vijay_sc, jaitin, merrick, karunakaran, acosta, hildreth, spatial, sun, angueira, lipidlab, leipzig_om, leipzig_sc, rosen_sn_om, rosen_sn_sc, rosen_sc, lipid_1s)

#filter to only retain significantly upregulated genes per cluster
#tested different FC cutoffs
FC <- 0
for (i in 1:length(vijay_om)) {
  vijay_om[[i]] <- vijay_om[[i]][which(vijay_om[[i]]$avg_log2FC > FC & vijay_om[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(vijay_sc)) {
  vijay_sc[[i]] <- vijay_sc[[i]][which(vijay_sc[[i]]$avg_log2FC > FC & vijay_sc[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(jaitin)) {
  jaitin[[i]] <- jaitin[[i]][which(jaitin[[i]]$avg_log2FC > FC & jaitin[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(merrick)) {
  merrick[[i]] <- merrick[[i]][which(merrick[[i]]$avg_log2FC > FC & merrick[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(karunakaran)) {
  karunakaran[[i]] <- karunakaran[[i]][which(karunakaran[[i]]$avg_log2FC > FC & karunakaran[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(acosta)) {
  acosta[[i]] <- acosta[[i]][which(acosta[[i]]$avg_log2FC > FC & acosta[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(hildreth)) {
  hildreth[[i]] <- hildreth[[i]][which(hildreth[[i]]$avg_log2FC > FC & hildreth[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(spatial)) {
  spatial[[i]] <- spatial[[i]][which(spatial[[i]]$avg_logFC > FC & spatial[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(sun)) {
  sun[[i]] <- sun[[i]][which(sun[[i]]$avg_log2FC > FC & sun[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(angueira)) {
  angueira[[i]] <- angueira[[i]][which(angueira[[i]]$avg_log2FC > FC & angueira[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(rosen_sn_om)) {
  rosen_sn_om[[i]] <- rosen_sn_om[[i]][which(rosen_sn_om[[i]]$avg_log2FC > FC & rosen_sn_om[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(rosen_sn_sc)) {
  rosen_sn_sc[[i]] <- rosen_sn_sc[[i]][which(rosen_sn_sc[[i]]$avg_log2FC > FC & rosen_sn_sc[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(rosen_sc)) {
  rosen_sc[[i]] <- rosen_sc[[i]][which(rosen_sc[[i]]$avg_log2FC > FC & rosen_sc[[i]]$p_val_adj < 0.05),]
}
for (i in 1:length(lipid_1s)) {
  lipid_1s[[i]] <- lipid_1s[[i]][which(lipid_1s[[i]]$avg_log2FC > FC & lipid_1s[[i]]$p_val_adj < 0.05),]
}

#get all possible intersections
#can probably be written much shorter
output <-data.frame(Cohort_1=integer(),
                    Cohort_2=integer(),
                    Number=integer(),
                    Percent_1=integer(),
                    Percent_2=integer(),
                    stringsAsFactors=FALSE)
output2 <- output
#vijay sc VS ALL
for (i in 2:length(vijay_sc)) {
  for (j in 2:length(vijay_om)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, vijay_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(vijay_om[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("B_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(jaitin)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, jaitin[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(jaitin[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("C_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(merrick)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, merrick[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(merrick[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("D_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(karunakaran)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, karunakaran[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(karunakaran[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("E_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(acosta)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, acosta[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(acosta[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("F_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 1:length(hildreth)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, hildreth[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(hildreth[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("G_",j-1)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(sun)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, sun[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(sun[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("H_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(angueira)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, angueira[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(angueira[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("I_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(vijay_sc[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_sc[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("A_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#vijay om VS ALL (- vijay sc)
for (i in 2:length(vijay_om)) {
  
  for (j in 2:length(jaitin)) {
    overlaps <- intersect(vijay_om[[i]]$gene, jaitin[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(jaitin[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("C_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(merrick)) {
    overlaps <- intersect(vijay_om[[i]]$gene, merrick[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(merrick[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("D_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(karunakaran)) {
    overlaps <- intersect(vijay_om[[i]]$gene, karunakaran[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(karunakaran[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("E_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(acosta)) {
    overlaps <- intersect(vijay_om[[i]]$gene, acosta[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(acosta[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("F_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 1:length(hildreth)) {
    overlaps <- intersect(vijay_om[[i]]$gene, hildreth[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(hildreth[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("G_",j-1)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(sun)) {
    overlaps <- intersect(vijay_om[[i]]$gene, sun[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(sun[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("H_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(angueira)) {
    overlaps <- intersect(vijay_om[[i]]$gene, angueira[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(angueira[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("I_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(vijay_om[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(vijay_om[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(vijay_om[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(vijay_om[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(vijay_om[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(vijay_om[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(vijay_om[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(vijay_om[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(vijay_om[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("B_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#jaitin VS ALL (-2)
for (i in 2:length(jaitin)) {
  
  for (j in 2:length(merrick)) {
    overlaps <- intersect(jaitin[[i]]$gene, merrick[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(merrick[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("D_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(karunakaran)) {
    overlaps <- intersect(jaitin[[i]]$gene, karunakaran[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(karunakaran[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("E_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(acosta)) {
    overlaps <- intersect(jaitin[[i]]$gene, acosta[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(acosta[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("F_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 1:length(hildreth)) {
    overlaps <- intersect(jaitin[[i]]$gene, hildreth[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(hildreth[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("G_",j-1)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(sun)) {
    overlaps <- intersect(jaitin[[i]]$gene, sun[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(sun[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("H_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(angueira)) {
    overlaps <- intersect(jaitin[[i]]$gene, angueira[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(angueira[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("I_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(jaitin[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(jaitin[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(jaitin[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(jaitin[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(jaitin[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(jaitin[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(jaitin[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(jaitin[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(jaitin[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("C_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#merrick VS ALL (-3)
for (i in 2:length(merrick)) {
  
  for (j in 2:length(karunakaran)) {
    overlaps <- intersect(merrick[[i]]$gene, karunakaran[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(karunakaran[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("E_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(acosta)) {
    overlaps <- intersect(merrick[[i]]$gene, acosta[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(acosta[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("F_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 1:length(hildreth)) {
    overlaps <- intersect(merrick[[i]]$gene, hildreth[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(hildreth[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("G_",j-1)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(sun)) {
    overlaps <- intersect(merrick[[i]]$gene, sun[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(sun[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("H_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(angueira)) {
    overlaps <- intersect(merrick[[i]]$gene, angueira[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(angueira[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("I_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(merrick[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(merrick[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(merrick[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(merrick[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(merrick[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(merrick[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(merrick[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(merrick[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(merrick[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("D_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#karunakaran VS ALL (-4)
for (i in 2:length(karunakaran)) {
  
  for (j in 2:length(acosta)) {
    overlaps <- intersect(karunakaran[[i]]$gene, acosta[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(acosta[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("F_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 1:length(hildreth)) {
    overlaps <- intersect(karunakaran[[i]]$gene, hildreth[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(hildreth[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("G_",j-1)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(sun)) {
    overlaps <- intersect(karunakaran[[i]]$gene, sun[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(sun[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("H_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(angueira)) {
    overlaps <- intersect(karunakaran[[i]]$gene, angueira[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(angueira[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("I_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(karunakaran[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(karunakaran[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(karunakaran[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(karunakaran[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(karunakaran[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(karunakaran[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(karunakaran[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(karunakaran[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(karunakaran[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("E_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#acosta VS ALL (-5)
for (i in 2:length(acosta)) {
 
    for (j in 1:length(hildreth)) {
    overlaps <- intersect(acosta[[i]]$gene, hildreth[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(hildreth[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("G_",j-1)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(sun)) {
    overlaps <- intersect(acosta[[i]]$gene, sun[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(sun[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("H_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(angueira)) {
    overlaps <- intersect(acosta[[i]]$gene, angueira[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(angueira[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("I_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(acosta[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(acosta[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(acosta[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(acosta[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(acosta[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(acosta[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(acosta[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(acosta[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(acosta[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("F_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#hildreth VS ALL (-6)
for (i in 1:length(hildreth)) {
 
  for (j in 2:length(sun)) {
    overlaps <- intersect(hildreth[[i]]$gene, sun[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(sun[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("H_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(angueira)) {
    overlaps <- intersect(hildreth[[i]]$gene, angueira[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(angueira[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("I_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(hildreth[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(hildreth[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(hildreth[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(hildreth[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(hildreth[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(hildreth[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(hildreth[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("G_",i-1)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(hildreth[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(hildreth[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("G_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#sun VS ALL (-7)
for (i in 2:length(sun)) {

  for (j in 2:length(angueira)) {
    overlaps <- intersect(sun[[i]]$gene, angueira[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(angueira[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("I_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(sun[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(sun[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(sun[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(sun[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(sun[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(sun[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(sun[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(sun[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(sun[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("H_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#angueira VS ALL (-8)
for (i in 2:length(angueira)) {

  for (j in 2:length(lipidlab)) {
    overlaps <- intersect(angueira[[i]]$gene, lipidlab[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(angueira[[i]]) * 100
    perc2 <- numb/ nrow(lipidlab[[j]]) * 100
    output2[1,1] <- paste0("I_",i-2)
    output2[1,2] <- paste0("J_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(angueira[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(angueira[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("I_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(angueira[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(angueira[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("I_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(angueira[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(angueira[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("I_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(angueira[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(angueira[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("I_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(angueira[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(angueira[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("I_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(angueira[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(angueira[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("I_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(angueira[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(angueira[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("I_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#lipidlab VS ALL (-9)
for (i in 2:length(lipidlab)) {
  
  for (j in 2:length(leipzig_sc)) {
    overlaps <- intersect(lipidlab[[i]]$gene, leipzig_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(lipidlab[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_sc[[j]]) * 100
    output2[1,1] <- paste0("J_",i-2)
    output2[1,2] <- paste0("K_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(lipidlab[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(lipidlab[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("J_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(lipidlab[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(lipidlab[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("J_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(lipidlab[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(lipidlab[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("J_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(lipidlab[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(lipidlab[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("J_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(lipidlab[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(lipidlab[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("J_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(lipidlab[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(lipidlab[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("J_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#leipzig_sc VS ALL (-10)
for (i in 2:length(leipzig_sc)) {

  for (j in 2:length(leipzig_om)) {
    overlaps <- intersect(leipzig_sc[[i]]$gene, leipzig_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_sc[[i]]) * 100
    perc2 <- numb/ nrow(leipzig_om[[j]]) * 100
    output2[1,1] <- paste0("K_",i-2)
    output2[1,2] <- paste0("L_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(leipzig_sc[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_sc[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("K_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(leipzig_sc[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_sc[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("K_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(leipzig_sc[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_sc[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("K_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(leipzig_sc[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_sc[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("K_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(leipzig_sc[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_sc[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("K_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#leigpzig om VS ALL (-11)
for (i in 2:length(leipzig_om)) {
  
   for (j in 2:length(rosen_sn_om)) {
    overlaps <- intersect(leipzig_om[[i]]$gene, rosen_sn_om[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_om[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_om[[j]]) * 100
    output2[1,1] <- paste0("L_",i-2)
    output2[1,2] <- paste0("M_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(leipzig_om[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_om[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("L_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(leipzig_om[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_om[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("L_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(leipzig_om[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_om[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("L_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(leipzig_om[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(leipzig_om[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("L_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#rosen sn sc  VS ALL (-12)
for (i in 2:length(rosen_sn_om)) {
 
  for (j in 2:length(rosen_sn_sc)) {
    overlaps <- intersect(rosen_sn_om[[i]]$gene, rosen_sn_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sn_om[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sn_sc[[j]]) * 100
    output2[1,1] <- paste0("M_",i-2)
    output2[1,2] <- paste0("N_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(rosen_sn_om[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sn_om[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("M_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(rosen_sn_om[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sn_om[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("M_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(rosen_sn_om[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sn_om[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("M_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#rosen sn  om VS ALL (-13)
for (i in 2:length(rosen_sn_sc)) {
  
  for (j in 2:length(rosen_sc)) {
    overlaps <- intersect(rosen_sn_sc[[i]]$gene, rosen_sc[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sn_sc[[i]]) * 100
    perc2 <- numb/ nrow(rosen_sc[[j]]) * 100
    output2[1,1] <- paste0("N_",i-2)
    output2[1,2] <- paste0("O_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(rosen_sn_sc[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sn_sc[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("N_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(rosen_sn_sc[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sn_sc[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("N_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#rosen sc VS ALL (-14)
for (i in 2:length(rosen_sc)) {
 
  for (j in 2:length(lipid_1s)) {
    overlaps <- intersect(rosen_sc[[i]]$gene, lipid_1s[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sc[[i]]) * 100
    perc2 <- numb/ nrow(lipid_1s[[j]]) * 100
    output2[1,1] <- paste0("O_",i-2)
    output2[1,2] <- paste0("P_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
  for (j in 2:length(spatial)) {
    overlaps <- intersect(rosen_sc[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(rosen_sc[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("O_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}
#spatial vs ALL (-15)
for (i in 2:length(lipid_1s)) {
  for (j in 2:length(spatial)) {
    overlaps <- intersect(lipid_1s[[i]]$gene, spatial[[j]]$gene)
    numb <- length(overlaps)
    perc <- numb/ nrow(lipid_1s[[i]]) * 100
    perc2 <- numb/ nrow(spatial[[j]]) * 100
    output2[1,1] <- paste0("P_",i-2)
    output2[1,2] <- paste0("Q_",j-2)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output <- rbind(output, output2)
  }
}

output3 <- output[which(output$Cohort_1 != output$Cohort_2),]
write.table(output3, file.path(DIR_RES,"marker_overlaps_v3.txt"), sep="\t", row.names = F, col.names = T)
#read cluster tables for color & size of nodes
setwd(file.path(DIR_WD))
vijay_sc_c <- read.table("./scSEQ_SVF/grundberg/gr_sc_Counts.txt", header = T)
vijay_om_c <- read.table("./scSEQ_SVF/grundberg/gr_om_Counts.txt", header = T)
jaitin_c <- read.table("./jaitin_Counts.txt", header = T)
merrick_c <- read.table("./Merrick_Counts.txt", header = T)
karanukaran_c <- read.table("./karuna_all_Counts.txt", header = T)
acosta_c <- read.table("./acosta_Counts.txt", header = T)
hildreth_c <- read.table("./hildreth_Counts.txt", header = T)
sun_c <- read.table("./snSeq_AT/Sun_Wolfrum/Sun_Counts.txt", header = T)
angueira_c <- read.table("./snSeq_AT/Angueira_Seale/seale_Counts.txt", header = T)
lipidlab_c <- read.table("./snSeq_AT/lipidlab/lipidlab_Counts.txt", header = T)
pamela_sc_c <- read.table("./snSeq_AT/blueher/leipzig_sc_Counts.txt", header = T)
pamela_om_c <- read.table("./snSeq_AT/blueher/leipzig_Counts.txt", header = T)
rosen_sn_om_c <-read.table("./snSeq_AT/rosen/rosen_sn_omental_Counts.txt", header = T)
rosen_sn_sc_c <-read.table("./snSeq_AT/rosen/rosen_sn_other_Counts.txt", header = T)
rosen_sc_c <-read.table("./snSeq_AT/rosen/rosen_sc_Counts.txt", header = T)
lipid_1s_c <- read.table("./snSeq_AT/lipid_sn_sc/lipidlab_Counts.txt", header = T)
#add cluster annotation to match intercept table
vijay_sc_c$cluster_new <- paste0("A_", vijay_sc_c$seurat_clusters)
vijay_om_c$cluster_new <- paste0("B_", vijay_om_c$seurat_clusters)
jaitin_c$cluster_new <- paste0("C_", jaitin_c$seurat_clusters)
merrick_c$cluster_new <- paste0("D_", merrick_c$seurat_clusters)
karanukaran_c$cluster_new <- paste0("E_", karanukaran_c$seurat_clusters)
acosta_c$cluster_new <- paste0("F_", acosta_c$seurat_clusters)
hildreth_c$cluster_new <- paste0("G_", hildreth_c$seurat_clusters)
sun_c$cluster_new <- paste0("H_", sun_c$seurat_clusters)
angueira_c$cluster_new <- paste0("I_", angueira_c$seurat_clusters)
lipidlab_c$cluster_new <- paste0("J_", lipidlab_c$seurat_clusters)
pamela_sc_c$cluster_new <- paste0("K_", pamela_sc_c$seurat_clusters)
pamela_om_c$cluster_new <- paste0("L_", pamela_om_c$seurat_clusters)
rosen_sn_om_c$cluster_new <- paste0("M_", rosen_sn_om_c$seurat_clusters)
rosen_sn_sc_c$cluster_new <- paste0("N_", rosen_sn_sc_c$seurat_clusters)
rosen_sc_c$cluster_new <- paste0("O_", rosen_sc_c$seurat_clusters)
lipid_1s_c$cluster_new <- paste0("P_", lipid_1s_c$seurat_clusters)
count_color <- rbind(vijay_sc_c, vijay_om_c, jaitin_c, merrick_c, karanukaran_c, acosta_c, hildreth_c,
                     sun_c, angueira_c, lipidlab_c, pamela_sc_c, pamela_om_c, rosen_sn_om_c, rosen_sn_sc_c, rosen_sc_c, lipid_1s_c)
library(dplyr)
#calculate jaccard distance
output3$total <- output3$Number/output3$Percent_1 * 100
output3$total_2 <- output3$Number/output3$Percent_2 * 100
jaccard <- output3[NULL,]
for (i in unique(output3$Cohort_1)) {
  temp <- output3[output3$Cohort_1 == i,]
  temp$total <- max(temp$total, na.rm = T)
  jaccard <- rbind(jaccard,temp)
} 
jaccard2 <- jaccard[NULL,]
for (i in unique(jaccard$Cohort_2)) {
  temp <- jaccard[jaccard$Cohort_2 == i,]
  temp$total_2 <- max(temp$total_2, na.rm = T)
  jaccard2 <- rbind(jaccard2, temp)
}
jaccard2$jaccard <- jaccard2$Number / (jaccard2$total + jaccard2$total_2 - jaccard2$Number)
write.table(jaccard2, file.path(DIR_RES,"Jaccard_All.txt"), col.names = T, row.names = F)
#jaccard2 <- read.table("./Results_II/Jaccard_All.txt", header = T)
#plot jaccard
j_1 <- jaccard2[,c(1,8)]
colnames(j_1) <- c("cluster","jaccard")
j_2 <- jaccard2[,c(2,8)]
colnames(j_2) <- c("cluster","jaccard")
jac_plot <- rbind(j_1, j_2)
jac_plot <- jac_plot %>% group_by(cluster) %>% slice_max(jaccard, n=5)
jac_plot$cohort <- substr(jac_plot$cluster, 1, 1)
library(ggplot2)
jac_plot$cohort <- factor(jac_plot$cohort, levels = c("J","K","P","H","N","O","G","A","F","D","Q","L","M","C","B","E","I"))
ggplot(data = jac_plot, aes(x=cohort,y=jaccard)) +geom_boxplot(outlier.shape = NA)+ theme_classic() +coord_flip()+scale_x_discrete(limits=rev)
  geom_jitter(color="black", size=1, alpha=0.7, width = 0.2)  + theme_classic() +coord_flip()

output3$avgPERCENT <- (output3$Percent_1 + output3$Percent_2) / 2

#output_NUMBER <- output3[NULL,]
#for (i in 1:length(nodes2)) {
#  output_temp <- output3[which(output3$Cohort_1 == nodes2[[i]] | output3$Cohort_2 == nodes2[[i]]),]
#  output_temp <- output_temp %>% top_n(5, Number)
#  output_NUMBER <- rbind(output_NUMBER, output_temp)
#  }
### NUMBER not good selection
#remove Hildreth dataset as previous mapping showed very low overlap
output3 <- output3 %>% filter(!grepl("G_\\d", output3$Cohort_1))
output3 <- output3 %>% filter(!grepl("G_\\d", output3$Cohort_2))
nodes2 <- unique(c(output3$Cohort_1, output3$Cohort_2))

#Problem: if no cutoff, network is overloaded with very weak connections
#tested different cutoff based on percentage & top connections
#final minimum requirements for edge connections were >15% genes overlapped in one of the two nodes and >5% in both nodes
output_PERCENT <- output3[NULL,]
for (i in 1:length(nodes2)) {
  output_temp <- output3[which(output3$Cohort_1 == nodes2[[i]] | output3$Cohort_2 == nodes2[[i]]),]
  output_temp <- output_temp %>% top_n(5, avgPERCENT)
  output_PERCENT <- rbind(output_PERCENT, output_temp)
}

output_PERCENT <- output_PERCENT[which(output_PERCENT$Percent_1 > 15 | output_PERCENT$Percent_2 > 15 ),]
output_PERCENT <- output_PERCENT[which(output_PERCENT$Percent_1 > 5 & output_PERCENT$Percent_2 > 5),]
#remove duplicates
output_PERCENT <- output_PERCENT[!duplicated(output_PERCENT),]

#generate edge and node tables
nodes <- as.data.frame(unique(c(output3$Cohort_1, output3$Cohort_2)))
colnames(nodes) <- "label" 
nodes$id <- row.names(nodes)
nodes <- select(nodes, id, label)
nodes <- nodes %>%
  left_join(count_color, by=c("label"="cluster_new"))
nodes <- select(nodes, id, label, cluster_color, percent, Freq, cluster_anno, cluster_group, cohort)
colnames(nodes) <- c("id", "label", "color", "cells_percent", "cells_count", "annotation","annot_group", "cohort")

edges <- output_PERCENT%>% 
  left_join(nodes, by = c("Cohort_1" = "label")) %>% 
  dplyr::rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("Cohort_2" = "label")) %>% 
  dplyr::rename(to = id)
edges <- select(edges, from, to, Number, Percent_1, Percent_2)
colnames(edges) <- c("from", "to", "weight", "Percent_1", "Percent_2")
nodes <- nodes[!is.na(nodes$label),]
edges <- edges[!is.na(edges$to),]
edges <- edges[!is.na(edges$from),]

#Export to cytoscape using igraph 
library(igraph)
routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
plot(routes_igraph)
library(RCy3)
cytoscape <- createNetworkFromIgraph(
  routes_igraph,
  title = "TOP5 PERCENT NO HILDRETH FILTERED",
  collection = "Meta New")

#final visual adjustments and network clustering using cytoscape