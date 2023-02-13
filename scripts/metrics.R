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
  require(ggplot2)
  require(patchwork)
  require(cluster)
  require(fpc)
  require(clusterSim)
  require(tidyverse)
  require(glue)
  
  
  # At this point one should have already done QC, Normalization with SCT and PCA
  
  # Create and empty list that stores the results. Each Hyperparameter will be an element
  # of the list and each element will have another list as element
  
  res <- list()
  
  # Load in a list of params
  params_list <- list(num_neighbors        = number_neighbors,
                      dist_metrics         = distance_metrics,
                      algorithm_clustering = algorithm_clustering,
                      min_distance         = min_distance_nn,
                      cluster_resolution   = cluster_resolution)
  
  # cluster resolution is omitted in the loops cause it'll  be the primary output
  # of each UMAP
  
  #for (i in names(params_list[-length(params_list)])) {
  for(j in length(params_list$algorithm_clustering)) {
    for(k in length(params_list$dist_metrics)) {
      for(l in length(params_list$num_neighbors)) {
        for(m in length(params_list$min_distance)) {
          
          # RUN UMAP, Neighbors finding and Clustering
          seurat_tmp <- RunUMAP(
            seurat_object,
            min.dist = params_list$min_distance[[m]],
            n.neighbors = params_list$num_neighbors[[l]],
            reduction.name = "UMAP",
            reduction.key = "UMAP_",
            dims = 1:n_pcs,
            n.components = 2,
            seed.use = 100
          )
          
          # Run Neighbors and Clustering
          seurat_tmp <- FindNeighbors(seurat_tmp, 
                                      dims = 1:n_pcs, 
                                      k.param = params_list$num_neighbors[[l]],
                                      annoy.metric = params_list$dist_metrics[[k]])
          
          seurat_tmp <- FindClusters(seurat_tmp, 
                                     resolution = params_list$cluster_resolution,
                                     algorithm = algorithm_clustering[[j]])
          
          
          # Define count matrix and resolutions to evaluate
          
          count_matrix <- seurat_tmp@reductions$pca@cell.embeddings
          resolutions_to_evaluate <- names(seurat_tmp@meta.data) %>% str_subset(pattern = "^SCT")
          
          # Silhouette Score
          
          Silhouette_scores <- map(resolutions_to_evaluate[2:length(resolutions_to_evaluate)], 
                                   ~ mean(cluster::silhouette(as.numeric(seurat_tmp@meta.data[[.x]]), dist(count_matrix))[,3])) %>%
            set_names(resolutions_to_evaluate[2:length(resolutions_to_evaluate)])
          
          # Calinski Harabasz Index
          
          CH_index <- map(resolutions_to_evaluate[2:length(resolutions_to_evaluate)],
                          ~ fpc::calinhara(x = count_matrix, clustering = as.numeric(seurat_tmp@meta.data[[.x]]))) %>% 
            set_names(resolutions_to_evaluate[2:length(resolutions_to_evaluate)])
          
          # Davies Bouldin Index
          
          DB_index <- map(resolutions_to_evaluate[2:length(resolutions_to_evaluate)],
                          ~ index.DB(x = count_matrix, cl = (as.numeric(seurat_tmp@meta.data[[.x]])))) %>% 
            map(pluck, 1) %>% 
            set_names(resolutions_to_evaluate[2:length(resolutions_to_evaluate)])
          
          # Reset Seurat Object
          
          seurat_tmp <- NULL
          
          # Create a list that containes the three indexes
          
          tmp_indexes_df <- list(Silhouette_scores = Silhouette_scores,
                                 CH_index = CH_index,
                                 DB_index = DB_index)
          
          # Append the three metrics to a list
          #tmp name
          tmp_name <- glue::glue("algo", "_", j, "_",
                                 "distMet", "_", k, "_",
                                 "nn", "_", l, "_",
                                 "minDist", "_", m)
          
          res[[tmp_name]] <- tmp_indexes_df
        }
      }
    }
  }
  
  return(res)
  
}
#}




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