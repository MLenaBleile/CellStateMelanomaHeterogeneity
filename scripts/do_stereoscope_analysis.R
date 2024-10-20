library(Seurat)

setwd("H:/Biostatistics/PICI/data/ST/BleileMS")
sampleid="PPD1_22"
load(paste("H:/Biostatistics/PICI/data/ST/sct/", sampleid, sep=""))
source("functions/scSpatialFeaturePlot.R")
stereo = read.csv("data/stereo_projections_subset.csv", row.names=1)
old.meta.data=spat@meta.data
colnames(stereo) = paste("stereo",colnames(stereo), sep=".")
spat@meta.data=cbind(old.meta.data, stereo)
stromafeats=c(paste("stereo",c("Stroma","Immune","Keratinocytes"), sep="."),
              c("Stroma","Immune","Keratinocytes"))
spat$stereo.predicted.id = colnames(stereo)[apply(stereo,1, which.max)]
spat$stereo.prop.tumor=rowSums(stereo[,setdiff(colnames(stereo),c("stereo.Stroma","stereo.Immune","stereo.Keratinocytes"))])
p1=ScSpatialFeaturePlot(spat, uq=.95, lq=.05,LegLabel = "Proportion",
                     features="stereo.Neural_Crest_like", flip=3)

p2=ScSpatialFeaturePlot(spat, uq=.95, lq=.05,LegLabel = "Proportion",
                     features="stereo.Melanocytic", flip=3)
# p3 =ScSpatialFeaturePlot(spat, uq=.95, lq=.05,LegLabel = "Proportion",
#                          features="Stress..Hypoxia.Response.", flip=3)

gridExtra::grid.arrange(p1,p2, ncol=2)

ScSpatialFeaturePlot(spat, uq=.99, lq=.0,LegLabel = "Stereoscope\nPrediction Score",
                     features=c("stereo.Melanocytic"), flip=3)

ScSpatialFeaturePlot(spat, uq=.99, lq=.00,LegLabel = "Stereoscope\nPrediction Score",
                     features=c("stereo.Neural_Crest_like"), flip=3)


col.vec=c("orange","orchid", "lightblue","darkblue")
names(col.vec) = paste("stereo",c("Melanocytic","Neural_Crest_like","Immune", "Stroma"), sep=".")
names(col.vec)[6:8] = setdiff(unique(spat$stereo.predicted.id), names(col.vec))
label.cols=rep("black",length(col.vec))
names(label.cols) = names(col.vec)
label.cols[c("Neural_Crest_like","Stress..Hypoxia.Response.")]="white"
SpatialDimPlot(spat, group.by="stereo.predicted.id", pt.size.factor=2.8, repel=T,label=T, cols=col.vec, label.color = label.cols)+NoLegend()




