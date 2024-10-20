library(Seurat)

##requires: 1) cleanreference (CM data) and 2) sct data
setwd("H:/Biostatistics/PICI/data/ST/BleileMS")
library(ggplot2)
load("H:/Biostatistics/PICI/data/ST/SmithyMS/ChrisMarineData/data/cleanreference")
source("functions/ScSpatialFeaturePlot.R")


load("H:/Biostatistics/PICI/data/ST/sct/PPD1_22")
bothvarfeats=intersect(VariableFeatures(cref),VariableFeatures(spat))
# This part takes awhile to run
anchors <- FindTransferAnchors(reference = cref, query = spat,
                               normalization.method = "SCT", features=bothvarfeats)
# save(anchors, file="Github data/anchors_full")
load("Github data/anchors_full")

predictions.matrix <- TransferData(anchorset = anchors, refdata = cref$Cluster, prediction.assay = F,
                                   weight.reduction = spat[["pca"]], dims = 1:15)




colnames(predictions.matrix) = stringr::str_replace(colnames(predictions.matrix), pattern="prediction.score.", "")

#save(predictions.matrix,file="Github Data/seurat_predictions_matrix")

# load("Github Data/seurat_predictions_matrix")
mal.cols=setdiff(colnames(predictions.matrix),
                 c("max", "predicted.id","Immune",
                   "Keratinocytes", "Stroma")
)
spat@meta.data=cbind(spat@meta.data, predictions.matrix)
spat$projected.tumor.proportion = rowSums(predictions.matrix[,mal.cols])



ScSpatialFeaturePlot(spat, features=mal.cols, pt.size.factor=2.8, flip=3)

###there are only 13 Antigen presentation spots, so we eliminate those
table(predictions.matrix$predicted.id)

cell.states=c("Melanocytic","Neural_Crest_like","Immune", "Stroma")
spat$seurat.predicted.id = cell.states[apply(predictions.matrix[,cell.states], 1, which.max)]

col.vec=c("orange","orchid", "lightblue","darkblue")
names(col.vec)= cell.states
label.cols=c("black","black","black","white")
names(label.cols) = names(col.vec)

SpatialDimPlot(spat, group.by="seurat.predicted.id", pt.size.factor=2.8, repel=T,label=T, cols=col.vec, label.color = label.cols)+NoLegend()

