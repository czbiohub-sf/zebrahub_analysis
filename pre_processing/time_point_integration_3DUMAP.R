library(stringr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(viridis)

options(future.globals.maxSize = 6000 * 1024^2)


ROOT_DIR = '/mnt/ibm_lg/alejandro/danio-atlas/zebrahub/'
file_path = paste0(ROOT_DIR, 'review_objects/v1/')
export_dir = paste0(ROOT_DIR, '/r_objects/') # Seurat objects (exported from Scanpy)

file_list = list.files(file_path)
file_list = as.array(file_list)
timepoint_list = str_split(file_list, "_", simplify = T)[,3]


# Split the timepoints into early and late --- for sorting 
# This is not necessary if the timepoint have numeric names 
late_points = grep(timepoint_list, pattern = 'dpf', value = T) 
early_points = timepoint_list[ !timepoint_list %in% late_points]

sorted_early = paste(sort(as.numeric(str_split(early_points, "s", simplify= T)[,1]) ) , "s", sep = "") 

sorted_late = paste(sort(as.numeric(str_split(late_points, "d", simplify= T)[,1]) ) , "dpf", sep = "") 

# all timepoints in order 
sorted_timepoints = c(sorted_early, sorted_late) 

# Make file names for loading
sorted_files = c() 
for(i in 1:length(sorted_timepoints))
    sorted_files[i] = grep(file_list , pattern = paste0("_", sorted_timepoints[i],"_"), value = T) 


# Load Seurat Objects 
seurat_list <- list() 

for(i in 1:length(sorted_timepoints)){
    print(paste("Loading ", sorted_timepoints[i], "from", export_dir))
    load( paste0(export_dir, sorted_timepoints[i], "/seurat_obj.rds") )    
    seurat_list[[i]] <- seurat_obj   
}

print(paste("Finished loading all objects..."))


# Integration functions 

normalize_list <-function(danio.list = list() , n_features = 2000){
    danio.list <- lapply(X = danio.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n_features)
    })


    features <- SelectIntegrationFeatures(object.list = danio.list)

    danio.list <- lapply(X = danio.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
    })
 
    return( list(danio.list, features) ) 
}

integration_rpca <- function(danio.list = list() , 
                             features = c() , 
                             use_dims = 1:30, 
                             anchor_n = 3, 
                             use_reference = NULL ){ 
    danio.anchors <- FindIntegrationAnchors(object.list = danio.list, anchor.features = features,
                                           normalization.method = 'LogNormalize', #c("LogNormalize", "SCT"),
                                           dims = use_dims, # default 1:30
                                           k.anchor = anchor_n, #default 5
                                           k.filter = 200, #default 200
                                           reduction = "rpca", # default cca 
                                           reference = use_reference
                                           )
    danio.combined <- IntegrateData(anchorset = danio.anchors)
    
    return(danio.combined)
    
}


# Visualization functions 
# Assumes the Seurat object contains a 3D UMAP embedding 
plot3DUMAP <- function( adata_obj = c(), 
                       plot_vars =   c("UMAP_1", "UMAP_2", "UMAP_3", 'cell_type_annotation_V3', 'seurat_clusters') ,
                       color_list = c("lightseagreen",
                              "gray50",
                              "darkgreen",
                              "red4",
                              "red",
                              "turquoise4",
                              "black",
                              "yellow4", 'peachpuff3'
                         ), point_size = 2 , 
                          sub_sample = F , n_cells = 10000){
    plot.data <- FetchData(object = adata_obj, vars = plot_vars)
    # label each cell by their identify 
    plot.data$label <- paste(rownames(plot.data), plot.data[, plot_vars[4] ], plot.data$seurat_clusters)

    plot.data$ident <- plot.data[, plot_vars[4]]

    # for visualization purposes, we can subsample the dataset. The structure of the UMAP should be still visible
    if(sub_sample){
       plot.data <- plot.data[sample(row.names(plot.data), n_cells) , ] 
    }

    fig <- plot_ly(data = plot.data, 
                   x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
                   color = ~ident, 
                   colors = color_list,
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = point_size, width=2), # controls size of points
                   text=~label, #This is that extra column we made earlier for which we will use for cell ID
                   hoverinfo="text", 
                   width = 800, height = 600) #When you visu



    # xaxis
    axx <- list(
        nticks = 4,
        range = c(-20,20) #select range of xaxis
    )

    # yaxis
    axy <- list(
        nticks = 4,
        range = c(-20,20) #select range of yaxis
    )

    #zaxis
    axz <- list(
        nticks = 4,
        range = c(-20,20) #select range of zaxis
    )

    fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
    fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
    
    
   # fig_cube <- fig_cube %>% layout(width = 800, height = 600)
    fig_cube

    
}


# specify the column with the cell type annotations 
# not necessary for newer versions 
cell_type_col <- c('cell_ontology_class','cell_ontology_class','cell_ontology_class','ontology','Ontology_annotation_xz',
                         'Ontology.annotation.xz', "Ontology_annotation_xz","new_cell_annotation","cell_annotation","annotations_keir")


# 1. Normalize 
res_list = normalize_list( seurat_list )
seurat_list_norm = res_list[[1]]
features = res_list[[2]]


# 2. Integration 
# We will use the middle timepoints as anchors for integration 
seurat.combined <- integration_rpca(danio.list = seurat_list_norm, 
                                    features = features, 
                                    use_reference = c(4,5,6))


# 3. Standard workflow (PCA + UMAP + clustering)
print(paste("PCA on integrated data..."))
# specify that we will perform downstream analysis on the corrected counts. NOTE that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
seurat.combined <- RunPCA(seurat.combined, npcs = 100, verbose = FALSE)


# We are making a 3D UMAP
print(paste("UMAP ... "))
seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:30,
                         metric='cosine',
                         n.components = 3, 
                         n.neighbors = 10,
                         local.connectivity  =1, # 1 default
                         repulsion.strength = 1, # 1 default
                         )

seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:30, annoy.metric = 'cosine')
seurat.combined <- FindClusters(seurat.combined, resolution = 1)


# 4. Visualization 
make_plot <- function(seurat.combined = c() ){
    viridis_cols <- mako(10)
    names(viridis_cols) <- seurat.combined@meta.data$timepoint %>% unique( )


    plot3DUMAP(seurat.combined, 
            plot_vars =   c("UMAP_1", "UMAP_2", "UMAP_3", 'timepoint') ,
            color_list = viridis_cols ) 

}

# 5. Save object
print("Exported object")
save(seurat.combined, file = paste(ROOT_DIR, 'r_objects/oct_manuscript_version/seurat_combined_ALL_UMAP3D_Apr2023.rds'))


# Export functions 

export_pca <- function(seurat.combined = c() , 
                       PCA_DIR = '/mnt/ibm_lg/alejandro/danio-atlas/zebrahub/embeddings/PCA_projection/'){

    pca_timepoint = Embeddings(seurat.combined, reduction = 'pca')
    write.csv(pca_timepoint, file = paste0(PCA_DIR,
                                         'pca_FullAtlas_all_cells.csv'), quote = F, row.names =T)    
}

export_h5ad <- function(seurat.combined = c(), 
    H5_DIR = '/mnt/ibm_lg/alejandro/danio-atlas/zebrahub/integrated/', 
    file_name ='zf_atlas_integratedv2' ){
    
    export_h5ad(seurat_object = seurat.combined, 
           EXP_DIR = H5_DIR, 
           h5ad_name = file_name)
}

export_pca_timepoint <- function(seurat.combined = c() ,
    PCA_DIR = '/mnt/ibm_lg/alejandro/danio-atlas/zebrahub/embeddings/PCA_projection/'){ 

    # We want to use the integrated assay
    DefaultAssay(seurat.combined) <- 'integrated'
    # Split by timepoint
    nmp.list_global = SplitObject(seurat.combined, split.by='timepoint')
    # sort in biological order 
    # save as independent PCA objects -- used later in Aligned UMAP
    for(f in 1:length(nmp.list_global)){
        pca_timepoint = Embeddings(nmp.list_global[[f]], reduction = 'pca')
        write.csv(pca_timepoint, file = paste0(PCA_DIR,
                                            'pca_',
                                            sorted_timepoints[f],
                                            '.csv'), quote = F, row.names =T)    
    }

    #Export metadata 
    for(f in 1:length(nmp.list_global)){
        meta_data = nmp.list_global[[f]]@meta.data
        write.csv(meta_data, file = paste0(PCA_DIR,
                                            'metadata/meta_',
                                            sorted_timepoints[f],
                                            '.csv'), quote = F, row.names =T)    
    }

}