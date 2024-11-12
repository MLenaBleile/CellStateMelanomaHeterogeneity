setwd("H:/Biostatistics/PICI/data/ST/BleileMS")
#load("data/PPD1_22")
sampleid="PPD1_22"
load(paste("H:/Biostatistics/PICI/data/ST/sct/", sampleid, sep=""))

load("data/ldas")
source("../PipelineMS/functions/scSpatialFeaturePlot.R")

K=8
optLDA = ldas$models[[K-4]]
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1)
deconProp <- results$theta
colnames(deconProp) = paste("LDA_Factor",colnames(deconProp), sep="")
deconGexp <- results$beta
#colnames(pos) = c("y","x")
predictions=colnames(deconProp)[apply(deconProp,1, which.max)]
names(predictions)=rownames(deconProp)
spat$LDA.predicted.id = predictions


spat@meta.data =cbind(spat@meta.data[rownames(deconProp),], deconProp)


ScSpatialFeaturePlot(spat, uq=.99, lq=.00,LegLabel = "Proportion",
                     features=colnames(deconProp), flip=3)
