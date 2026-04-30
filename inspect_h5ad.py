import scanpy as sc
import scipy.sparse as sp
import scipy.io as sio
import numpy as np

file = r"C:\code\NEMO\GW18 CGE\2bf59d36.033aaed3-7a29-4dc1-b3d8-a454b5686d52.h5ad"

print("Loading file...")
adata = sc.read_h5ad(file)

print("Converting matrix...")

X = adata.X
if sp.issparse(X):
    X = X.toarray()

print("Building MATLAB struct...")

out = {
    "X": X,
    "umap": adata.obsm["X_umap"],
    "tsne": adata.obsm["X_tsne"],
    "pca": adata.obsm["X_pca"],
    "genes": np.array(adata.var["gene_symbol"]),
    "louvain": np.array(adata.obs["louvain"]),
    "orig_louvain": np.array(adata.obs["orig_louvain"])
}

out_file = r"C:\code\NEMO\nemo_matlab_ready.mat"

print("Saving to:", out_file)
sio.savemat(out_file, out)

print("DONE FILE CREATED")