#Harmony
#setwd("")
packages <- c("Seurat", "patchwork", "harmony", "DT", "dplyr", "writexl")
invisible(lapply(packages, library, character.only=TRUE))

seurat <- readRDS("./Seurat_Merged.rds")
var_int <- "orig.ident"
#tested different variables to integrate over, harmony can integrate over multiple variables
seurat_harmony <- RunHarmony(object = seurat, group.by.vars = var_int, assay.use="SCT", reduction = "ica", plot_convergence = TRUE)
red_use <- "harmony"
dims_use <- 1:30
nneigh <- 50
seurat_harmony <- RunUMAP(object = seurat_harmony, 
                        reduction = red_use, 
                        dims = dims_use,
                        n.neighbors = nneigh)
#run UMAP on unintegrated data to compare
seurat_harmony <- RunUMAP(object = seurat_harmony, 
                        reduction = "ica", 
                        dims = dims_use, 
                        n.neighbors = nneigh,
                        reduction.name = "umapICA",
                        seed.use = 42)
seurat_harmony <- AddMetaData(seurat_harmony, seurat_harmony@reductions$umap@cell.embeddings, col.name = c("UMAP_1", "UMAP_2"))
paste("k (spots):", dim(seurat_harmony@reductions$umap@cell.embeddings)[1])
seurat_harmony <- FindNeighbors(object = seurat_harmony, verbose = T, reduction = red_use, dims = dims_use)
seurat_harmony <- FindClusters(object = seurat_harmony, verbose = T, algorithm = 1, resolution = 1)
DimPlot(seurat_harmony)
DimPlot(seurat_harmony, reduction = "umap", group.by = c("orig.ident", "tissue", "method"))
saveRDS(seurat_harmony, file = "./Seurat_Harmony.rds")
