library(Seurat)
library(ScSpatialFeaturePlot)
load("data/VisiumData.rda")
load("results/deconProp_K8.rda")



predictions=colnames(deconProp)[apply(deconProp,1, which.max)]
names(predictions)=rownames(deconProp)
spat$LDA.predicted.id = predictions

SpatialDimPlot(spat, group.by="LDA.predicted.id",)

for(one.feature in colnames(deconProp)){
    spat[[one.feature]] = deconProp[,one.feature]
}

ScSpatialFeaturePlot(spat, uq=.99, lq=.00,LegLabel = "Proportion",
                     features=colnames(deconProp), flip=3)

