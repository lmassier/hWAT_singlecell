#rPCA
#Integrate seurat
#setwd("")
packages <- c("Seurat", "patchwork", "harmony", "DT", "dplyr", "writexl")
invisible(lapply(packages, library, character.only=TRUE))

seurat <- readRDS("./Seurat_Merged.rds")

seurat.list <- SplitObject(seurat, split.by = "orig.ident")
rm(seurat)
gc()
seurat.list <- lapply(X = seurat.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
seurat.list <- lapply(X = seurat.list, FUN = RunPCA, features = features)
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT",
                                       anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
seurat.combined.sct <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", dims = 1:30)
seurat.combined.sct <- RunPCA(seurat.combined.sct, verbose = FALSE)
seurat.combined.sct <- RunUMAP(seurat.combined.sct, reduction = "pca", dims = 1:30)
DimPlot(seurat.combined.sct, reduction = "umap", group.by = c("orig.ident", "tissue", "method"))
saveRDS(seurat.combined.sct, file = "./Seurat_rPCA.rds")
seurat.combined.sct <- FindNeighbors(object = seurat.combined.sct, verbose = T, reduction = red_use, dims = dims_use)
seurat.combined.sct <- FindClusters(object = seurat.combined.sct, verbose = T, algorithm = 1, resolution = 1)