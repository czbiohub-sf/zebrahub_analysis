library(data.table)
library(dplyr)
library(Seurat)
library(viridis)
library(ggExtra)
library(ggplot2)



ATLAS_DIR = '/mnt/ibm_lg/alejandro/danio-atlas/atlas_objects/seurat_objects/'        # Read QC'd h5ad from here 
SCT_ATLAS = '/mnt/ibm_lg/alejandro/danio-atlas/atlas_objects/timepoint_integration/' # Write integrated Atlas (after SCT pipeline) here
ATLAS_FILE  = 'danio_atlas_annotated_38k.rds'

options(future.globals.maxSize = 6000 * 1024^2) # to avoid memory issues

# output files 
OUT_FILE = 'integration_2023.rdata'

setwd(ATLAS_DIR)

print(paste("Loading object from ", ATLAS_DIR ))
load(ATLAS_FILE)

# change the default name 
danio.combined <- seurat_obj
rm(seurat_obj) 

# Rename columns with typos - optional 
danio.combined@meta.data <- danio.combined@meta.data %>%
  mutate(
    timepoint = if_else(timepoint == "5somite", "05somite", timepoint),
    global_annotation = case_when(
      global_annotation == "Differentating_Neurons" ~ "Differentiating_Neurons",
      global_annotation == "Muscles" ~ "Muscle",
      global_annotation == "Somite" ~ "Somites",
      global_annotation == "Epiderm" ~ "Epidermal",
      TRUE ~ global_annotation
    )
  )

# Integration 
print(paste("Split object by timepoint"))
# 1. Split object and process each time point independently 
danio.list = SplitObject(danio.combined, split.by='timepoint')

danio.list <- lapply(X = danio.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})



features <- SelectIntegrationFeatures(object.list = danio.list)

danio.list <- lapply(X = danio.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

print(paste("Integration ... "))
# 2. Find integration anchors and integrate data 
danio.anchors <- FindIntegrationAnchors(object.list = danio.list, anchor.features = features,
                                       normalization.method = 'LogNormalize', #c("LogNormalize", "SCT"),
                                       dims = 1:50, # default 1:30
                                       k.anchor = 5, #default 5
                                       k.filter = 200, #default 200 for a query cell, If the anchor reference cell is found within the first k.filter (200) neighbors, then we retain this anchor.
                                       k.score = 30, # default 30: For each reference anchor cell, we determine its k.score (30) nearest within-dataset neighbors and its k.score nearest neighbors in the query dataset
                                       reduction = "rpca", # default cca, rpca should be faster 
                                       reference = c(1,6) 
                                       )
danio.combined <- IntegrateData(anchorset = danio.anchors)



# 3. Generate an integrated embedding: run PCA on integrated (corrected) counts 
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(danio.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
danio.combined <- ScaleData(danio.combined, verbose = FALSE)
danio.combined <- RunPCA(danio.combined, npcs = 100, verbose = FALSE)
print(paste("Runnning UMAP on the integrated PCA embedding..."))
# 4. UMAP on integrated embbedding 
danio.combined <- RunUMAP(danio.combined, reduction = "pca", dims = 1:30,
                         metric='euclidean',
                         n.neighbors = 30,
                         local.connectivity  =1, # 1 default
                         repulsion.strength = 1, # 1 default
                         )
danio.combined <- FindNeighbors(danio.combined, reduction = "pca", dims = 1:30)
danio.combined <- FindClusters(danio.combined, resolution = 0.5)


# 5. Export R object 
save(danio.combined, file = paste0(SCT_ATLAS, OUT_FILE)) 