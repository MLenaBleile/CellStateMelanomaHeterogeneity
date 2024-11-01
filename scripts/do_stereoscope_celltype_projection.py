import os
import tempfile
import matplotlib.pyplot as plt
import anndata
import muon
import numpy as np
from ray import tune
import ray
from scvi import autotune
import scanpy as sc
import scvi
import seaborn as sns
import torch
from scvi.external import RNAStereoscope, SpatialStereoscope

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

#st_adata = sc.read_visium(path="H:/Biostatistics/PICI/data/ST/spaceranger/Sample_PPD1-22_D1_IGO_14577_C_18/outs")
st_adata = sc.read(filename="H:/Biostatistics/PICI/data/ST/BleileMS/data/PPD1_22.h5ad")
st_adata.var_names_make_unique()

sc_adata = sc.read(filename="H:/Biostatistics/PICI/data/ST/BleileMS/data/reference_6type.h5ad")
sc_adata.var_names_make_unique()


print("# cells, # genes before filtering:", sc_adata.shape)

# let us filter some genes
G = 2000
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


RNAStereoscope.setup_anndata(sc_adata, layer="counts", labels_key="Cluster")

train = False
if train:
    sc_model = RNAStereoscope(sc_adata)
    sc_model.train(max_epochs=400)
    sc_model.history["elbo_train"][10:].plot()
    sc_model.save("scmodel400ep", overwrite=True)
else:
    sc_model = RNAStereoscope.load("scmodel400ep", adata=sc_adata)
    print("Loaded RNA model from file!")

#plt.savefig("H:/Biostatistics/PICI/data/ST/BleileMS/data/stereoscope_sc_model_400ep.png")

SpatialStereoscope.setup_anndata(st_adata, layer="counts")
train = True
if train:
    #SpatialStereoscope.load("stmodel_400sc_1250ep", overwrite=True)
    spatial_model = SpatialStereoscope.from_rna_model(st_adata, sc_model)
    spatial_model.train(max_epochs=2000)
    spatial_model.history["elbo_train"][10:].plot()
    spatial_model.save("stmodel_6type_400sc_2000ep", overwrite=True)
else:
    spatial_model = SpatialStereoscope.load("stmodel", adata=st_adata)
    print("Loaded Spatial model from file!")


plt.savefig("H:/Biostatistics/PICI/data/ST/BleileMS/data/stereoscope_st_model_6type_st_2000_sc_400.png")

st_adata.obsm["proportions"] = spatial_model.get_proportions()

st_adata.obsm["proportions"].to_csv("H:/Biostatistics/PICI/data/ST/BleileMS/data/stereoscope_6type_400sc_2000st.csv")
