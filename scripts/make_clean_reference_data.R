library(Seurat)
library(ggplot2)

####only run once
acells=readRDS("data/GC_all_cells.rds")
mcells = readRDS("data/Malignant_cells.rds")

##only want the cells in both
shared=intersect(colnames(acells), colnames(mcells))
st=table(acells$Cluster[shared])
immnames=setdiff(colnames(acells)[acells$Cluster=="Immune"], shared)
kernames=setdiff(colnames(acells)[acells$Cluster=="Keratinocytes"], shared)
stronames=setdiff(colnames(acells)[acells$Cluster=="Stroma"], shared)

malnames=intersect(names(acells$Cluster[acells$Cluster=="Malignant"]), 
                   colnames(mcells)[!mcells$Malignant_clusters %in% c("Patient_specific_A", "Patient_specific_B","Mitochondrial(low_quality)")])
acells$cell.type=acells$Cluster
cref=acells[,c(immnames, kernames, stronames, malnames)]
cref$Cluster[malnames] = mcells$Malignant_clusters[malnames]


newlabels=read.csv("data/marine_seuratobj_NONMAL_annot.csv", skip=1)
rownames(newlabels) = newlabels$name.of.cell
cells.in.both=intersect(newlabels$name.of.cell, colnames(cref))
cref$cell.type.alli = cref$Cluster
cref$cell.type.alli[cells.in.both] = newlabels[cells.in.both,"Alli.s.SingleR.cell.annotation.as.more.general.categories"]
table(cref$cell.type.alli)
##remove Neurons: this is not brain tissue so we know these labels are incorrect
cref = cref[,cref$cell.type.alli!="Neurons"]
save(cref, file="data/cleanreference")
