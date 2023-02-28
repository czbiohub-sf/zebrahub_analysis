# zebrahub_analysis

Contains the notebooks used for the following analyses:
* scRNAseq processing
* embryo gene expression variability
* RNA velocity
* State transition analysis from RNA velocity fields

## Organization
The structure of this repo is illustrated below. 
```
├── embryo_gene_expression_variability (contains notebooks for analysis)
│   ├── Inter-Embryo_Divergence_DataFrames.ipynb/ (notebook that creates mpKLD dataframes, and null model mpKLD dataframe)
│   ├── Inter-Embryo_Divergence_PathwayAnalysis_Plots.ipynb/ (notebook for cleaning, pathway analysis, plots found in Fig2a-2d)
│   ├── zf2mouse.txt/ (file with analogous genes for pathway analysis)
│   ├── README.md/
├── pre-processing (contains pre-processing, timepoint integration and clustering notebooks)
│   ├── Aligned_UMAPs_EarlyTimepoints.ipynb/ (notebook creates aligned UMAP for early timepoints, using k-NN for UMAP clusters in adjacent timepoints)
│   ├── Integrated_Embedding_3D_UMAPs.ipynb/ (notebook for intergration of 3D UMAPs)
│   ├── Integrated_Embedding_EarlyTimepoints.ipynb/ (notebook for integrated UMAP global embedding using Seurat)
│   ├── Sequencing_QualityControl.ipynb/ (notebook loads CellRanger files, generated adata, conducts QC/clustering/UMAP)
└── README.md
```


## Set up the Conda Environment

### 1. Pre Processing


### 2. Embryo Gene Expression Variability

![image](https://user-images.githubusercontent.com/35573897/221049586-b0bc8f05-b035-4279-9116-62ebf9b97c53.png)

We set up the conda environment to run the inter-embryo analysis notebooks by running the following commands in the terminal:

```
# Create the environment with Python version 3.8.13 to be compatible with other packages
conda create -c conda-forge python=3.8.13 -n zebrahub_inter-embryo_analysis

#Activate the conda environment
conda activate zebrahub_inter-embryo_analysis

# Install Jupyter Lab
conda install -c conda-forge jupyterlab

# Install necessary packages
conda install ipykernel

# we need matplotlib version 3.6 to be compatible with scanpy

conda install seaborn
conda install scanpy
conda install matplotlib=3.6
conda install -c conda-forge re2
conda install -c jmcmurray json
conda install -c anaconda requests

```
