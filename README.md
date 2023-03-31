# Zebrahub: RNA-sequencing analysis 


![image](https://user-images.githubusercontent.com/35573897/221049586-b0bc8f05-b035-4279-9116-62ebf9b97c53.png)
*Schematic of the analysis workflow to understand differences between genetically identical embryos. For each gene, we compare the distribution of counts across single-cells using the Kullback–Leibler divergence metric.*


<!--
<figure>
  <img src="https://user-images.githubusercontent.com/35573897/221049586-b0bc8f05-b035-4279-9116-62ebf9b97c53.png" alt="image">
  <figcaption>Schematic of the analysis workflow to understand differences between genetically identical embryos. For each gene, we compare the distribution of counts across single-cells using the Kullback–Leibler divergence metric. </figcaption>
</figure>
-->

Elucidating the developmental process of an organism will require the complete cartography of cellular lineages in the spatial, temporal, and molecular domains. We present Zebrahub, a comprehensive dynamic atlas of zebrafish embryonic development that combines single-cell sequencing time course data with light-sheet microscopy-based lineage reconstructions. Zebrahub is a foundational resource to study developmental processes at both transcriptional and spatiotemporal levels. It is publicly accessible as a web-based resource, providing an open-access collection of datasets and tools.<br>

This repository contains a collection of scripts and notebooks to analyze Zebrahub's single-cell RNA seq data and generate the corresponding figures as published in the manuscript. <br> For analysis of imaging data, please refere to https://github.com/royerlab/in-silico-fate-mapping <br>

The manuscript is  published here:  https://www.biorxiv.org/content/10.1101/2023.03.06.531398v1
Data is accesible through the Zebrahub portal: https://zebrahub.ds.czbiohub.org/
<br>
<br>

## Table of Contents
- [Repo organization](#organization)
- [Setting up an environment](#conda)
- [Credits](#credits)


## Organization

Each folder contains the notebooks used for the following analyses:
* scRNAseq processing (pre_processing + clustering)
* RNA velocity (scVelo + single-cell velocity quantification)
* State transition analysis from RNA velocity methods (state_transitions)
* embryo gene expression variability (Kullback–Leibler divergence)


The structure of this repo is illustrated below. 
```
├── pre-processing (contains pre-processing, timepoint integration and clustering notebooks)
│   ├── Aligned_UMAPs_EarlyTimepoints.ipynb/ (notebook creates aligned UMAP for early timepoints, using k-NN for UMAP clusters in adjacent timepoints)
│   ├── Integrated_Embedding_3D_UMAPs.ipynb/ (notebook for intergration of 3D UMAPs)
│   ├── Integrated_Embedding_EarlyTimepoints.ipynb/ (notebook for integrated UMAP global embedding using Seurat)
│   ├── Sequencing_QualityControl.ipynb/ (notebook loads CellRanger files, generated adata, conducts QC/clustering/UMAP)
├── RNA_velocity (contains RNA velocity analyses)
│   ├── RNA_Velocity_FullAtlas.ipynb/ (notebook creates RNA velocity graph and visualization for all early timepoints)
│   ├── RNA_velocity_NMP_lineages.ipynb/ (notebook creates RNA velocity graph and visualization for NMP lineages, individual timepoints)
├── state_transitions (contains state transitions analysis from RNA velocity notebook)
│   ├── RNA_velocity_transition_graph.ipynb/ (notebook creates state transition graph from single-cell RNA velocity)
├── embryo_gene_expression_variability (contains mpKLD analysis notebooks)
│   ├── Inter-Embryo_Divergence_DataFrames.ipynb/ (notebook that creates mpKLD dataframes, and null model mpKLD dataframe)
│   ├── Inter-Embryo_Divergence_PathwayAnalysis_Plots.ipynb/ (notebook for cleaning, pathway analysis, plots found in Fig2a-2d)
│   ├── zf2mouse.txt/ (file with analogous genes for pathway analysis)
│   ├── README.md/
└── README.md
```


## Setting up a conda environment <a name="conda"></a>


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
