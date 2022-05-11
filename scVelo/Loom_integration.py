# -*- coding: utf-8 -*-
"""

"""

from os.path import isfile, join

import scvelo as scv
import scanpy as sp
import pandas as pd
import numpy as np

data_file_path = '/home/hilmar/Work/Charite/Tacke_group/ABC/loom/'

# First step: combine individual files to single loom file
# import glob
# import loompy

# all_RV_files = [join(data_file_path, f) for f in glob.glob(data_file_path+"/A*.loom") if isfile(join(data_file_path, f))]

# loompy.combine(files = all_RV_files, output_file = join(data_file_path, "BCN_NAFLD_scVelo_combined.loom"), key="Accession")

# combined scVelo data from all samples
complete_dataset_file = "BCN_NAFLD_scVelo_combined.loom"
complete_data = scv.read(join(data_file_path, complete_dataset_file), cache=True)
complete_data.var_names_make_unique()

# Add UMAP from Seurat data integration
umap_from_seurat = pd.read_csv(join(data_file_path, "BCN_MoMF_DC_combined_annotated_umap.txt"), sep="\t")

cell_ids = [i.replace("_", ":").replace("-1", "x") for i in umap_from_seurat['Unnamed: 0'] ] 
umap_from_seurat['cell_ids_fixed'] = cell_ids

complete_data = complete_data[cell_ids,:]

umap_from_seurat.set_index("cell_ids_fixed", inplace=True)
umap_from_seurat_sorted = umap_from_seurat.loc[complete_data.obs_names,:]
tmp = np.stack([umap_from_seurat_sorted["UMAP_1"], umap_from_seurat_sorted["UMAP_2"]]).T
complete_data.obsm['X_umap'] = tmp


metadata_from_seurat = pd.read_csv(join(data_file_path, 'BCN_MoMF_DC_only_combined_annotated_metadata.txt'), sep="\t")
cell_ids = [i.replace("_", ":").replace("-1", "x") for i in metadata_from_seurat['Unnamed: 0'] ] 
metadata_from_seurat['cell_ids_fixed'] = cell_ids

metadata_from_seurat.set_index("cell_ids_fixed", inplace=True)
metadata_from_seurat_sorted = metadata_from_seurat.loc[complete_data.obs_names,:]

complete_data.obs['clusters'] = metadata_from_seurat_sorted['cluster_label']
complete_data.obs['louvain'] = metadata_from_seurat_sorted['cluster_label']

 
scv.pp.filter_and_normalize(complete_data, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(complete_data, n_pcs=30, n_neighbors=30)

#    scv.tl.velocity(complete_data, mode='stochastic')

scv.tl.recover_dynamics(complete_data, n_jobs=4)
scv.tl.velocity(complete_data, mode='dynamical', n_jobs=4)

scv.tl.velocity_graph(complete_data)
    
#sp.tl.louvain(complete_data)
#scv.tl.umap(complete_data)

scv.pl.velocity_embedding_stream(complete_data, basis='umap', save='RNAvelocity_stream_all_samples.png' )

#scv.pl.velocity_embedding(complete_data, basis='umap')
#scv.pl.velocity_embedding_grid(complete_data, basis='umap')
#scv.pl.velocity_embedding_stream(complete_data, basis='umap')
    