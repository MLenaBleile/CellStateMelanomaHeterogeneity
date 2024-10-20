library(reticulate)
library(sceasy)
library(anndata)
library(Seurat)

sc <- import("scanpy", convert = FALSE)


load("Github Data/PPD1_22")
load("Github Data/cleanreference")

adata <- convertFormat(spat, from="seurat", to="anndata", main_layer="counts", assay="Spatial",drop_single_values=FALSE)
adata

write_h5ad(
  anndata=adata,
  filename = "Github Data/PPD1_22.h5ad",
  compression = NULL,
  compression_opts = NULL,
  as_dense = list()
)
scdata <- convertFormat(cref, from="seurat", to="anndata", main_layer="counts", assay="RNA",drop_single_values=FALSE)


write_h5ad(
  anndata=scdata,
  filename = "Github Data/reference.h5ad",
  compression = NULL,
  compression_opts = NULL,
  as_dense = list()
)

desired.types = c("Neural_Crest_like","Melanocytic","Immune","Stroma")
cref.subset <- subset(cref, subset= Cluster %in% desired.types)

scdata2 <- convertFormat(cref.subset, from="seurat", to="anndata", main_layer="counts", assay="RNA",drop_single_values=FALSE)
write_h5ad(
  anndata=scdata2,
  filename = "Github Data/reference.subset.h5ad",
  compression = NULL,
  compression_opts = NULL,
  as_dense = list()
)
