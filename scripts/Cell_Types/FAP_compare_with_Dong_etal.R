#Comparison to Dong et al.
#Mandrup et al. not currently in the manuscript

#read packages and make functions
library(readxl)
library(tidyverse)
read_excel_netsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
library(nichenetr) #to tranfer mouse gene symbols to human 
dong <- read.delim("dong_marker.txt")
dong$human <- dong$gene %>% convert_mouse_to_human_symbols()
dong <- dong[!is.na(dong$human),]
fap_marker <- read_excel_netsheets("./MARKER/FIBRO_SC06.xlsx")
#fap_marker <- read_excel_netsheets("./MARKER/FIBRO_OM03.xlsx")

#filter gene lists with FC cutoff 0.5 to reduce, otherwise to many false positive results expected
fap_marker <- fap_marker[2:length(fap_marker)]
for (i in 1:length(fap_marker)) {
  fap_marker[[i]] <- fap_marker[[i]][which(fap_marker[[i]]$avg_log2FC > 0.5 & fap_marker[[i]]$p_val_adj < 0.05),]
}

#compare
dong2 <- split(dong, f=dong$cluster)
for (i in 1:length(dong2)) {
  dong2[[i]] <- dong2[[i]][which(dong2[[i]]$avg_log2FC > 0 & dong2[[i]]$p_val_adj < 0.05),]
}
output <-data.frame(Cohort_1=integer(),
                    Cohort_2=integer(),
                    Number=integer(),
                    Percent_1=integer(),
                    Percent_2=integer(),
                    Total_1 = integer(),
                    Total_2 = integer(),
                    stringsAsFactors=FALSE)
output2 <- output
for (i in 1:length(fap_marker)) {
  for (j in 1:length(dong2)) {
    overlaps <- intersect(fap_marker[[i]]$gene, dong2[[j]]$human)
    numb <- length(overlaps)
    perc <- numb/ nrow(fap_marker[[i]]) * 100
    perc2 <- numb/ nrow(dong2[[j]]) * 100
    output2[1,1] <- paste0("FAP_",i-1)
    output2[1,2] <- paste0("DONG_",j)
    output2[1,3] <- numb
    output2[1,4] <- perc
    output2[1,5] <- perc2
    output2[1,6] <- nrow(fap_marker[[i]])
    output2[1,7] <- nrow(dong2[[j]])
    output <- rbind(output, output2)
  }
}
View(output)
output$jaccard <- output$Number/(output$Total_1+output$Total_2-output$Number)
library(viridis)
output$Cohort_1 <- factor(output$Cohort_1, levels = unique(output$Cohort_1))
#ggplot(data = output, aes(x=Cohort_1, y=Cohort_2, fill=jaccard))+ geom_tile() +scale_fill_stepsn(colors=viridis(10), n.breaks=9)+
#  theme_classic()+ scale_x_discrete(expand=c(0,0))+ scale_y_discrete(expand=c(0,0))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggplot(data = output, aes(x=Cohort_1, y=Cohort_2, fill=Number))+ geom_tile() +scale_fill_stepsn(colors=viridis(10), n.breaks=10)+
#  theme_classic()+ scale_x_discrete(expand=c(0,0))+ scale_y_discrete(expand=c(0,0))

x <- split(output, output$Cohort_1)
y <- split(output, output$Cohort_2)

SCvsDong <- data.frame(row.names = y$DONG_1$Cohort_1)
for (i in 1:7) {
  SCvsDong[1:17,i] <- y[[i]]$jaccard
}
colnames(SCvsDong) <- x$FAP_0$Cohort_2
library(pheatmap)
color_func_purp <- colorRampPalette(colors = c("#f3f0ff", "#62376E"))

#Fig 4h
pheatmap(SCvsDong, color = viridis(8), border_color = NA)
pheatmap(SCvsDong, color = viridis(8), border_color = NA, cutree_rows = 2, cutree_cols = 2)
