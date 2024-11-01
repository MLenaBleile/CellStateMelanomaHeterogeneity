library(Seurat)

setwd("H:/Biostatistics/PICI/data/ST/BleileMS")
sampleid="PPD1_22"
load(paste("H:/Biostatistics/PICI/data/ST/sct/", sampleid, sep=""))
old.meta.data=spat@meta.data
source("functions/scSpatialFeaturePlot.R")
stereo = read.csv("data/stereoscope_6type_400sc_2000st.csv", row.names=1)

colnames(stereo) = paste("stereo",colnames(stereo), sep=".")
spat@meta.data=cbind(old.meta.data[,colnames(old.meta.data)],stereo)
stromafeats=c(paste("stereo",c("Stroma","Immune","Keratinocytes"), sep="."),
              c("Stroma","Immune","Keratinocytes"))
spat$stereo.predicted.id = colnames(stereo)[apply(stereo,1, which.max)]
spat$stereo.prop.tumor=rowSums(stereo[,setdiff(colnames(stereo),c("stereo.Stroma","stereo.Immune","stereo.Keratinocytes"))])



#Some Mesenchymal-like cells are dominant in this setup but not many
table(spat$stereo.predicted.id)


stereo.major.states=names(table(spat$stereo.predicted.id))
all.col.vec=c("orange","orchid", "lightblue","darkblue","darkred")
col.vec=all.col.vec[1:length(stereo.major.states)]
names(col.vec)[1:5] = paste("stereo",c("Melanocytic","Neural_Crest_like","Immune", "Stroma","Keratinocytes"), sep=".")

label.cols=rep("black",length(col.vec))
names(label.cols) = names(col.vec)
#label.cols[c("Immune")]="white"
SpatialDimPlot(spat, group.by="stereo.predicted.id",
               pt.size.factor=2.8, repel=T,label=F, cols=col.vec, label.color = label.cols)




