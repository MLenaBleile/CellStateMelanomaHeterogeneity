library(STdeconvolve)
library(Matrix)
library(Seurat)
#setwd("H:/Biostatistics/PICI/data/ST/BleileMS")
load("data/VisiumData.rda")
cd <- GetAssayData(spat, slot="counts", assay="Spatial")
pos <- GetTissueCoordinates(spat)

counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)
cleancorpus=t(as(corpus, "sparseMatrix"))
ldas <- fitLDA(cleancorpus, Ks = 5:12)

save(ldas, file="Github Data/ldas")

ldas <- fitLDA(cleancorpus, Ks = 2:4)

save(ldas, file="Github Data/ldas_2to4")
