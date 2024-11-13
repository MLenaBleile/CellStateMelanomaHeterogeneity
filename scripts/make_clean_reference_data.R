library(Seurat)
library(ggplot2)
setwd("data")

####only run once
acells=readRDS("GC_all_cells.rds")
mcells = readRDS("Malignant_cells.rds")

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
save(cref,file="Pozniak_reference_data.rda")

