import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict
from sklearn.neighbors import NearestNeighbors
import umap.aligned_umap
import sklearn.datasets


class timeSeriesUMAP:
    def __init__(self, data_dir: str, files: List[str], n_pcs: int = 15):
        self.data_dir = data_dir
        self.files = files
        self.n_pcs = n_pcs
        self.meta_list = []
        self.pca_list = []
        self.anchor_dict = []
        self.umap_coords = pd.DataFrame()

    def read_files(self):
        for i in range(len(self.files)):
            df = pd.read_csv(self.data_dir + 'nmps_meta_' + self.files[i] + '.csv', index_col=0)
            self.meta_list.append(df)
            df = pd.read_csv(self.data_dir + 'nmps_global_' + self.files[i] + '.csv', index_col=0)
            X = df.values
            self.pca_list.append(X[:, 0:self.n_pcs])

    def find_neighbors(self, max_k: int, frac_k: float, max_dist: float, use_metric: str, cell_type_label: str):
        for i in range(len(self.pca_list) - 1):
            Y = self.pca_list[i]
            X = self.pca_list[i + 1]
            print("finding neighbors in " + str(i) + ":" + str(i + 1))
            nbrs = NearestNeighbors(n_neighbors=1, metric=use_metric).fit(Y)
            distances, indices = nbrs.kneighbors(X)
            neigh_distribution = np.concatenate(distances, axis=0)
            neigh_indexes = np.concatenate(indices, axis=0)

            pairs = pd.DataFrame({'neighbor': neigh_indexes, 'dist': neigh_distribution})
            pairs.reset_index(inplace=True)
            pairs.rename(columns={'index': 'cell_target'}, inplace=True)
            pairs['cell_type'] = self.meta_list[i + 1][cell_type_label].values

            df1 = pairs.groupby(['cell_type'])
            df2 = df1.apply(lambda x: x.sort_values(["dist"]))
            df3 = df2.reset_index(drop=True)

            pairs_rank = df3.groupby('neighbor').head(1)
            pairs_rank = pairs_rank.groupby('cell_type').head(max_k)
            pairs_rank = pairs_rank[pairs_rank['dist'] < max_dist]
            pairs_dict = {pairs_rank['neighbor'].values[j]: pairs_rank['cell_target'].values[j] for j in
                          range(pairs_rank.shape[0])}
            self.anchor_dict.append(pairs_dict)

        print("Finished finding all anchors")

    def run_aligned_umap(self, k_neighbors : int = 20, align_regularisation : float = 0.01, align_window  : float = 2, epochs : int = 200, rand_state : int = 42, use_metric : str = 'cosine'):
        print("Running aligned UMAP")
        aligned_mapper = umap.AlignedUMAP(metric=use_metric,
                                                        n_neighbors=k_neighbors,
                                                        alignment_regularisation=align_regularisation,
                                                        alignment_window_size=align_window,
                                                        n_epochs=epochs,
                                                        random_state=rand_state, ).fit(self.pca_list,
                                                                                relations=self.anchor_dict)
        self.aligned_mapper = aligned_mapper

    def get_umap_coords(self, cell_type_label = str):
        print("Building UMAP coordinates")
        all_timepoints = []
        for i in range(0, len(self.files)):
            aligned_umap_coord = pd.DataFrame({
                'UMAP_1': self.aligned_mapper.embeddings_[i].T[0],
                'UMAP_2': self.aligned_mapper.embeddings_[i].T[1], 
                'timepoint' : self.files[i],
                'cell_type' : self.meta_list[i][cell_type_label].values, 
                'cell_id' : self.meta_list[i].cell_id.values
            })
            all_timepoints.append(aligned_umap_coord)
        
        self.umap_coords = pd.concat(all_timepoints)
        return(self.umap_coords)