###---------------------------------------------------------------------###
###                         INTRODUCTION 								###
###---------------------------------------------------------------------###

# Implementation of silhoutte scores, Calinski-Harabasz Index and Davies-Bouldin Index
# For now the output is always splitted based on the resolution for clustering
# The next thing to do is to evaluate which hyperparameters need to be fine tuned


###---------------------------------------------------------------------###
###                       Function for the Metrics 						###
###---------------------------------------------------------------------###

# Function to evaluate different hyperparameters
test_clust <- function(seurat_object, number_neighbors, distance_metrics, algorithm_clustering, min_distance_nn,
                       n_pcs, cluster_resolution) {

  # Packages
  require(Seurat)
  library(ggplot2)
  library(patchwork)
  library(cluster)
  library(fpc)
  library(clusterSim)
  library(tidyverse)


  # At this point onw should have already done QC, Normalization with SCT and PCA

  # Load in a list of params
  params_list <- list(num_neighbors        = number_neighbors,
                      dist_metrics         = distance_metrics,
                      algorithm_clustering = algorithm_clustering,
                      min_distance         = min_distance_nn,
                      cluster_resolution   = cluster_resolution)

  # cluster resolution is omitted in the loops cause it'll  be the primary output
  # of each UMAP

  for (i in names(params_list[-length(params_list)])) {
    for(j in length(params_list$min_distance)) {


  # RUN UMAP, Neighbors finding and Clustering
  seurat <- RunUMAP(
    seurat,
    min.dist = paramms_list$min_distance[[j]],
    n.neighbors = params_list$num_neighbors[[i]],
    reduction.name = "UMAP",
    reduction.key = "UMAP_",
    dims = 1:n_pcs,
    n.components = 2,
    seed.use = 100
  )

  seurat <- FindNeighbors(seurat, dims = 1:n_pcs, k.param = params_list$num_neighbors[1])
  seurat <- FindClusters(seurat, resolution = params_list$cluster_resolution[1])

  # Perform UMAP with the first parameters


  # Define count matrix and resolutions to evaluate

  count_matrix <- seurat@reductions$pca@cell.embeddings
  resolutions_to_evaluate <- names(seurat@meta.data) %>% str_subset(pattern = "^SCT")

  # Silhouette Score

Silhouette_scores_V2 <- map(resolutions_to_evaluate_v2[2:length(resolutions_to_evaluate_v2)], 
                        ~ mean(cluster::silhouette(as.numeric(seurat_V2@meta.data[[.x]]), dist(count_matrix_V2))[,3])) %>%
  set_names(resolutions_to_evaluate_v2[2:length(resolutions_to_evaluate_v2)])

# Calinski Harabasz Index

CH_index_V2 <- map(resolutions_to_evaluate_v2[2:length(resolutions_to_evaluate_v2)],
                ~ fpc::calinhara(x = count_matrix_V2, clustering = as.numeric(seurat_V2@meta.data[[.x]]))) %>% 
  set_names(resolutions_to_evaluate_v2[2:length(resolutions_to_evaluate_v2)])

# Davies Bouldin Index

DB_index_V2 <- map(resolutions_to_evaluate_v2[2:length(resolutions_to_evaluate_v2)],
                ~ index.DB(x = count_matrix_V2, cl = (as.numeric(seurat_V2@meta.data[[.x]])))) %>% 
  map(pluck, 1) %>% 
  set_names(resolutions_to_evaluate_v2[2:length(resolutions_to_evaluate_v2)])
    }
  }
}



###---------------------------------------------------------------------###
###                      Best Candidates Function 						###
###---------------------------------------------------------------------###


top_candidates_metrics <- function(df_silhouette, df_CH, df_DB){

  require(tidyr)
  require(magrittr)
  require(dplyr)
  require(purrr)

  #create a list of the DFs
  df_list <- list(silhouette_score = df_silhouette,
                  CH_index        = df_CH,
                  DB_index        = df_DB)

  # Join them together
  temp <- df_list %>%
    map(as_tibble) %>%
    map(~ .x %>%
          pivot_longer(everything())) %>%
    imap(~ rename(.x, "{.y}" := value)) %>%
    reduce(left_join)

  #Extract the best resolutions by metric
  silh_res <- temp %>% arrange(desc(silhouette_score)) %>%
    select(name) %>% rename(silhouette_score = name)

  CH_res <- temp %>% arrange(desc(CH_index)) %>%
    select(name) %>% rename(CH_index = name)

  DB_res <- temp %>% arrange(DB_index) %>%
    select(name) %>% rename(DB_index = name)

  # joining ordered resolutions
  result <- data.frame(silh_res, CH_res, DB_res)

  return(result)
}



###---------------------------------------------------------------------###
###                         Plotting Function 						    ###
###---------------------------------------------------------------------###

hyperparameter <- "PCA Components"

pdf(here::here(OUTDIR, DOCNAME, "test_clustering.pdf"), width = 20, height = 8)
# Silhouette scores
Silhouette_scores %>% 
  as_tibble %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = name, y = value, group = 1)) +
  geom_line(col = "darkviolet", linewidth = 2) +
  theme_custom +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        legend.position = "none") +
  labs(x = "Cluster Resolution", y = "Silhouette Score", title = glue("Silhouette score for {hyperparameter} hyperparameter")) |
  

# CH Index

CH_index %>% 
  as_tibble() %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = name, y = value, group = 1)) +
  geom_line(col = "darkviolet", linewidth = 2) +
  theme_custom +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        legend.position = "none") +
  labs(x = "Cluster Resolution", y = "Calinski-Harabasz Index", title = glue("CH index for {hyperparameter} hyperparameter")) |
  
  
# Db Index
  
DB_index %>% 
  as_tibble() %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x = name, y = value, group = 1)) +
  geom_line(col = "darkviolet", linewidth = 2) +
  theme_custom +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        legend.position = "none") +
  labs(x = "Cluster Resolution", y = "Davies-Bouldin Index", title = glue("DB index for {hyperparameter} hyperparameter"))

dev.off()