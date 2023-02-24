# Embryo Gene Expression Variability

This folder contains the notebooks for **inter-embryo divergence analysis**. 
We developed a framework based on the Kullback-Leibler divergence (KLD) to quantify the differences across siblings for each gene.
![image](https://user-images.githubusercontent.com/35573897/221049586-b0bc8f05-b035-4279-9116-62ebf9b97c53.png)

## 1. The notebook *Inter-Embryo_Divergence_DataFrames.ipynb* contains the code to create the mpKLD dataframe and null model dataframe, executing the framework seen in the panels above. 
* We upload the anndata object for a specific timepoint, which contains the cellxgene matrix with log-normalized RNA counts for all embryos. 
* For a specific gene and cell cluster, we use the Kullback-Leibler divergence to compare the distributions, which have been smoothed using a Gaussian kernel. 
After running the notebook, the output is a csv with columns including gene, timepoint, cluster (cell type), 
and median (the median pairwise Kullback-Leibler divergence score, mpKLD).

<img width="626" alt="image" src="https://user-images.githubusercontent.com/35573897/221313426-93f3d0aa-4c30-42e3-a226-5a49eea35264.png">

* The second part of the notebook is the null model for the mpKLD scores for a **single timepoint**. 
To establish a negative control for the difference between embryos, we randomly reassign the embryo-ID for each cell, 
and then re-run the analysis of comparing the gene distributions. We repeat this process n=20 times.
The null model is output in a csv containing the same columns as the mpKLD dataframe, 
but includes a column 'Random_State' stating which randomized iteration the analysis comes from.

<img width="643" alt="image" src="https://user-images.githubusercontent.com/35573897/221313827-01ac9429-8468-4849-b81f-9f78db6b7c7c.png">



## 2. The notebook *Inter-Embryo_Divergence_PathwayAnalysis_Plots.ipynb* contains the code to create many of the plots seen in the paper (Fig. and Supp.), as well as filtering of the mpKLD dataframe and pathway analysis.
* Depending on the desired analysis, the notebook requires uploading the mpKLD dataframe and null model dataframe from part (1) and the anndata object.
* We calculate the gene metrics from the anndata object, and then filter. Genes with mean counts less than 0.1 are removed, and genes with fraction of 0s higher than 0.8 are also removed.
We remove ribosomal and mitochondrial genes, only keeping the RNA genes.
* If pathway analysis is desired, we can identify the top n genes with high (most divergent between embryos) and low (most synchronous between embryos) mpKLD scores,
and then perform pathway analysis using Enrichr on the top n high and low genes. The zf2mouse.txt file is a list of analogous genes used for pathway analysis.
This yields a dataframe with pathway information, including columns for the pathway, p-value, list of overlapping genes, database, timepoint, and gene type (high or low inter-fish/mpKLD).
* The last part of the notebook includes code to make the following figures found in the figure and supplementary:
  * The first plot is a demonstration of a few specific genes, showing their gene count distributions per embryo and heatmap of pairwise Kullback-Leibler scores.

  <img width="700" alt="image" src="https://user-images.githubusercontent.com/35573897/221312139-321f3746-6d85-49a3-9431-311a33e50f3e.png">
 
  * The second plot is a multiplot that compares the null model to the original data.

  <img width="771" alt="image" src="https://user-images.githubusercontent.com/35573897/221311925-6c75fc8b-b741-47ea-b45d-773cee3f60fd.png">

  * This notebook also contains the code for the Kolmogorov-Smirnov test between adjacent timepoints.
  
  <img width="410" alt="image" src="https://user-images.githubusercontent.com/35573897/221312393-6da00134-1d9b-474d-becb-06def49c7d04.png">

