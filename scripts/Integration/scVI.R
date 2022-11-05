#scVI
#following https://docs.scvi-tools.org/en/0.15.1/tutorials/notebooks/scvi_in_R.html
#restart R/ Rstudio for scVI using reticulate 
library(reticulate)
#use_python("define/python/environment", required=T)
library(sceasy)
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
#setwd("")
packages <- c("Seurat", "patchwork", "harmony", "DT", "dplyr", "writexl")
invisible(lapply(packages, library, character.only=TRUE))
seurat <- readRDS("Seurat_Merged.rds")
top2000 <- head(VariableFeatures(seurat), 2000) #other variable features possible, depending on expected differences between batches
seurat_top <- seurat[top2000]
adata <- convertFormat(seurat_top, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata) # Note generally in Python, dataset conventions are obs x var
# run setup_anndata
scvi$model$SCVI$setup_anndata(adata, batch_key = 'orig.ident') #batch key
#note: also possible to add further categorigal or continuous_covariate_keys, did not improve integration for this project
# create the model
model = scvi$model$SCVI(adata)
# train the model
model$train()
#  get the latent represenation
latent = model$get_latent_representation()

# put it back in the original Seurat object
latent <- as.matrix(latent)
rownames(latent) = colnames(seurat_top)
seurat[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_top))

seurat <- RunUMAP(seurat, dims = 1:10, reduction = "scvi", n.components = 2)

saveRDS(seurat, "./Seurat_scVI.rds")
rm(seurat)
