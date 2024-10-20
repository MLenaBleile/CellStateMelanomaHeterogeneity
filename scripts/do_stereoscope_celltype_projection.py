import os
import tempfile
import matplotlib.pyplot as plt
import anndata
import muon
import numpy as np
#import pooch
import scanpy as sc
import scvi
import seaborn as sns
import torch
from scvi.model import CondSCVI, DestVI

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

#st_adata = sc.read_visium(path="H:/Biostatistics/PICI/data/ST/spaceranger/Sample_PPD1-22_D1_IGO_14577_C_18/outs")
st_adata = sc.read(filename="H:/Biostatistics/PICI/data/ST/BleileMS/data/PPD1_22.h5ad")
st_adata.var_names_make_unique()

sc_adata = sc.read(filename="H:/Biostatistics/PICI/data/ST/BleileMS/data/reference.mel.reduced.h5ad")
sc_adata.var_names_make_unique()


print("# cells, # genes before filtering:", sc_adata.shape)

# let us filter some genes
G = 1000
sc.pp.filter_genes(sc_adata, min_counts=10)

sc_adata.layers["counts"] = sc_adata.X.copy()

sc.pp.highly_variable_genes(
    sc_adata, n_top_genes=G, subset=True, layer="counts", flavor="seurat_v3"
)
print("# cells, # genes after filtering:", sc_adata.shape)


st_adata.layers["counts"] = st_adata.X.copy()

sc.pp.normalize_total(st_adata)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata

# filter genes to be the same on the spatial data
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)

#sc.pl.embedding(st_adata, basis="spatial", color="seurat_clusters", s=80)


CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="Cluster")

sc_model = CondSCVI(sc_adata, weight_obs=False)
sc_model.view_anndata_setup()

sc_model.train()

sc_model.history["elbo_train"].iloc[5:].plot()
#plt.show()

DestVI.setup_anndata(st_adata, layer="counts")

st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.view_anndata_setup()

st_model.train(max_epochs=2500)
st_model.history["elbo_train"].iloc[10:].plot()
plt.show()

st_adata.obsm["proportions"] = st_model.get_proportions()

st_adata.obsm["proportions"].to_csv("H:/Biostatistics/PICI/data/ST/BleileMS/data/stereo_projections_mel_reduced.csv")
