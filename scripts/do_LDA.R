library(STdeconvolve)
# library(Matrix)
# library(Seurat)

# load("data/VisiumData.rda")
# cd <- GetAssayData(spat, slot="counts", assay="Spatial")
# pos <- GetTissueCoordinates(spat)

load("data/cd")
load("data/pos")

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
cleancorpus=t(as.matrix(corpus))
ldas <- fitLDA(cleancorpus, Ks = 5:12)

save(ldas, file="models/ldas")

K=8
optLDA = ldas$models[[K-4]]
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1)
deconProp <- results$theta
colnames(deconProp) = paste("LDA_Factor",colnames(deconProp), sep="")
save(deconProp, file="results/deconProp_K8.rda")

ldas <- fitLDA(cleancorpus, Ks = 2:4)

save(ldas, file="models/ldas_2to4")


K=4
optLDA = ldas_2to4$models[[K]]
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1)
deconProp <- results$theta
colnames(deconProp) = paste("LDA_Factor",colnames(deconProp), sep="")
save(deconProp, file="results/deconProp_K4.rda")
