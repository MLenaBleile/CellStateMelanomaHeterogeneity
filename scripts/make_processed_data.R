library(Seurat)
library(ggplot2)

sampleid="in_house_sample_id"

home.dir="in_house_home_dir"
one.file= "in_house_sample_location"

imgd=Read10X_Image(
  image.dir=paste(home.dir,one.file,"/outs/spatial",sep=""),
  image.name = "tissue_lowres_image.png"
)

spat= Load10X_Spatial(
  data.dir=paste(home.dir,one.file,"/outs",sep=""),
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "img",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = imgd,
)

##Compute percentage of mitochondrial reads
spat[["percent.mt"]] <- PercentageFeatureSet(spat, pattern = "^MT+")

##Filter out spots with too few reads
spat <- subset(spat, subset = nFeature_Spatial > 100)

##Correct for library size and adjust for percent mitochondrial reads
spat <- SCTransform(spat, assay = "Spatial", verbose = FALSE,method = "glmGamPoi",return.only.var.genes = F, vars.to.regress = c("percent.mt"))

##Compute principal components
spat <- RunPCA(spat, assay = "SCT", verbose = FALSE)

#Check how many PCs we seem to need
ElbowPlot(spat)

spat <- FindNeighbors(spat, reduction = "pca", dims = 1:15)
spat <- FindClusters(spat, verbose = FALSE)
spat <- RunUMAP(spat, reduction = "pca", dims = 1:15)

save(spat, file="VisiumData.rda")
