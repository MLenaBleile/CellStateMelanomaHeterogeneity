library(reticulate)
library(sceasy)
library(anndata)
library(Seurat)

sc <- import("scanpy", convert = FALSE)

sampleid="PPD1_22"
load(paste("H:/Biostatistics/PICI/data/ST/sct/", sampleid, sep=""))
load("H:/Biostatistics/PICI/data/ST/SmithyMS/ChrisMarineData/data/cleanreference")

adata <- convertFormat(spat, from="seurat", to="anndata", main_layer="counts", assay="Spatial",drop_single_values=FALSE)
adata

write_h5ad(
  anndata=adata,
  filename = "H:/Biostatistics/PICI/data/ST/BleileMS/data/PPD1_22.h5ad",
  compression = NULL,
  compression_opts = NULL,
  as_dense = list()
)


desired.types = c("Neural_Crest_like","Melanocytic","Mesenchymal_like","Immune","Stroma")
cref.subset <- subset(cref, subset= Cluster %in% desired.types)

scdata2 <- convertFormat(cref.subset, from="seurat", to="anndata", main_layer="counts", assay="RNA",drop_single_values=FALSE)
write_h5ad(
  anndata=scdata2,
  filename = "H:/Biostatistics/PICI/data/ST/BleileMS/data/reference_6type.h5ad",
  compression = NULL,
  compression_opts = NULL,
  as_dense = list()
)
