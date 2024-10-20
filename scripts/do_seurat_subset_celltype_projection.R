library(Seurat)

##requires: 1) cleanreference (CM data) and 2) sct data
setwd("H:/Biostatistics/PICI/data/ST/BleileMS")
library(ggplot2)
load("H:/Biostatistics/PICI/data/ST/SmithyMS/ChrisMarineData/data/cleanreference")
source("H:/Biostatistics/PICI/data/ST/PipelineMS/functions/ScSpatialFeaturePlot.R")

cref.subset=subset(cref, subset=(Cluster %in% c("Immune","Stroma","Neural_Crest_like","Melanocytic")))
load("data/PPD1_22")
# anchors <- FindTransferAnchors(reference = cref.subset, query = spat,
#                                normalization.method = "SCT")
 
# predictions.matrix <- TransferData(anchorset = anchors, refdata = cref.subset$Cluster, prediction.assay = F,
#                                    weight.reduction = spat[["pca"]], dims = 1:15)

# colnames(predictions.matrix) = stringr::str_replace(colnames(predictions.matrix), pattern="prediction.score.", "seurat.subset.")

# save(predictions.matrix,"data/seurat_subset_predictions_matrix")

load("data/seurat_subset_predictions_matrix")
mal.cols=setdiff(colnames(predictions.matrix),
                 c("seurat.subset.max", "predicted.id","seurat.subset.Immune",
                   "seurat.subset.Keratinocytes", "seurat.subset.Stroma")
)
spat@meta.data=cbind(spat@meta.data, predictions.matrix)
spat$projected.tumor.proportion = rowSums(predictions.matrix[,mal.cols])



ScSpatialFeaturePlot(spat, features=mal.cols, pt.size.factor=2.8, flip=3)
ScSpatialFeaturePlot(spat, features=c("seurat.subset.Neural_Crest_like", "seurat.subset.Melanocytic"), pt.size.factor=2.8, flip=3)


save(spat, file="data/PPD1_22")
