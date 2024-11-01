library(Seurat)

##requires: 1) cleanreference (CM data) and 2) sct data
setwd("H:/Biostatistics/PICI/data/ST/BleileMS")
library(ggplot2)
load("H:/Biostatistics/PICI/data/ST/SmithyMS/ChrisMarineData/data/cleanreference")
source("functions/ScSpatialFeaturePlot.R")

# save(cref, file="data/ChrisMarineReference")
#subset the data to only include exclusive cell types
cref.6type=subset(cref, subset=(Cluster %in% c("Immune","Stroma","Keratinocytes","Neural_Crest_like","Melanocytic", "Mesenchymal_like")))



load("H:/Biostatistics/PICI/data/ST/sct/PPD1_22")
bothvarfeats=intersect(rownames(cref),VariableFeatures(spat))
# This part takes awhile to run
anchors <- FindTransferAnchors(reference = cref.6type, query = spat,
                               normalization.method = "SCT")
# save(anchors, file="data/anchors")
# load("data/anchors")

predictions.matrix <- TransferData(anchorset = anchors, refdata = cref.6type$Cluster, prediction.assay = F,
                                   weight.reduction = spat[["pca"]], dims = 1:15)




colnames(predictions.matrix) = stringr::str_replace(colnames(predictions.matrix), pattern="prediction.score.", "")

#save(predictions.matrix,file="data/seurat_predictions_matrix")

# load("Github Data/seurat_predictions_matrix")
mal.cols=setdiff(colnames(predictions.matrix),
                 c("max", "predicted.id","Immune",
                   "Keratinocytes", "Stroma")
)

spat@meta.data=cbind(spat@meta.data[,1:8], predictions.matrix)
spat$projected.tumor.proportion = rowSums(predictions.matrix[,mal.cols])



ScSpatialFeaturePlot(spat, features=mal.cols, pt.size.factor=2.8, flip=3)

###No spots are dominated by Keratinocytes or Mesenchymal
table(predictions.matrix$predicted.id)



cell.states=c("Melanocytic","Neural_Crest_like","Immune", "Stroma")
spat$seurat.predicted.id = cell.states[apply(predictions.matrix[,cell.states], 1, which.max)]

col.vec=c("orange","orchid", "lightblue","darkblue")
names(col.vec)= cell.states

 
label.cols=c("black","black","white","black")
names(label.cols) = names(col.vec)

SpatialDimPlot(spat, group.by="seurat.predicted.id", pt.size.factor=2.8, repel=T,label=T, cols=col.vec, label.color = label.cols)+NoLegend()
