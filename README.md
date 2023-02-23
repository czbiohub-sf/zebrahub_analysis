# zebrahub_analysis

Contains the notebooks used for the following analyses:
* scRNAseq processing
* embryo gene expression variability
* RNA velocity
* State transition analysis from RNA velocity fields


## Set up the Conda Environment


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

## Organization
The structure of this repo is illustrated below. 
```
├── embryo_gene_expression_variability (contains notebooks for analysis, example dataframes and example adata)
│   ├── Inter-Embryo_Divergence_DataFrames.ipynb/ (notebook that creates mpKLD dataframes, and null model mpKLD dataframe)
│   ├── Inter-Embryo_Divergence_PathwayAnalysis_Plots.ipynb/ (notebook for cleaning, pathway analysis, plots found in Fig2a-2d)
│   ├── example_adata/ (h5ad, example anndata)
│   ├── example_interembryo_divergence_df.csv/ (example csv for 12hpf, neural anterior)
│   ├── example_nullmodel_df.csv/ (example csv for 12hpf, neural anterior)
│   ├── zf2mouse.txt/ (file with analogous genes for pathway analysis)
│   ├── README.md/
└── README.md
```

