###---------------------------------------------------------------------###
###                     hyperparameters evaluation 						###
###---------------------------------------------------------------------###


test_clust <- function(seurat_object, algorithm_clustering, distance_metrics, number_neighbors, min_distance_nn,
                       n_pcs, cluster_resolution) {
  
  # Packages
  require(Seurat)
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
    for(algorithm in params_list$algorithm_clustering) {
      for(distance_metric in params_list$dist_metrics) {
        for(numb_neighbors in params_list$num_neighbors) {
          for(minimum_dist in params_list$min_distance) {
            
            
            # RUN UMAP, Neighbors finding and Clustering
            seurat_tmp <- RunUMAP(
              seurat_object,
              min.dist = minimum_dist,
              n.neighbors = numb_neighbors,
              reduction.name = "UMAP",
              reduction.key = "UMAP_",
              dims = 1:n_pcs,
              n.components = 2,
              seed.use = 100
            )
            
            # Run Neighbors and Clustering
            seurat_tmp <- FindNeighbors(seurat_tmp, 
                                    dims = 1:n_pcs, 
                                    k.param = numb_neighbors,
                                    annoy.metric = distance_metric)
            
            seurat_tmp <- FindClusters(seurat_tmp, 
                                   resolution = params_list$cluster_resolution,
                                   algorithm = algorithm)
            
            
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
            tmp_name <- glue::glue("algo", "_", algorithm, "_",
                                   "distMet", "_", distance_metric, "_",
                                   "nn", "_", numb_neighbors, "_",
                                   "minDist", "_", minimum_dist)
            
            res[[tmp_name]] <- tmp_indexes_df
          }
        }
      }
    }
  
  return(res)
  
  }


###---------------------------------------------------------------------###
###                         Plotting Function 						    ###
###---------------------------------------------------------------------###

# This function takes in the output of test_clust and plots the metrics per combination

plot_clusters <- function(metrics, outdir){
  
  # Packages
  require(tidyverse)
  require(glue)
  
  
  #Set theme for later
  theme_custom <- theme_bw(base_size = 16) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black"),
          legend.key.size  = unit(0.4, units = "cm"),
    ) +
    theme(plot.title = element_text(hjust = 0.5)
    ) +
    theme(
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    theme(
      strip.text = element_text(colour = 'black'),
      strip.background = element_rect(fill="white")
    )
  
  ## Plotting function
  
  for (i in names(metrics)){
    
    plot <-   metrics[i] %>%
      as.data.frame() %>% 
      pivot_longer(everything()) %>% 
      mutate(split_col = case_when(str_detect(name, "Silhouette") ~ "Silhouette",
                                   str_detect(name, "CH") ~ "Calinski-Harabasz",
                                   TRUE ~ "Davies-Bouldin")) %>% 
      split(~ as_factor(.$split_col)) %>% 
      #separate(name, c("metric", "index1", "index2"), sep = "_", remove = FALSE) %>% 
      imap(~ ggplot(.x, 
                    aes(x = name, y = value, group = 1)) +
             geom_line(col = "darkviolet", linewidth = 2) +
             geom_point() +
             theme_custom +
             theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                   legend.position = "none") +
             labs(x = "Cluster Resolution", y = .x$split_col, title = glue("Score/Index for   {.y}")))
    
    pdf(glue::glue(outdir,"/", i, "_", "plot.pdf"))
    print(plot)
    dev.off()
    
  }
}

###---------------------------------------------------------------------###
###                      Best Candidates Function 						###
###---------------------------------------------------------------------###


# top_candidates_metrics <- function(df_silhouette, df_CH, df_DB){
  
#   require(tidyr)
#   require(magrittr)
#   require(dplyr)
#   require(purrr)
  
#   #create a list of the DFs
#   df_list <- list(silhouette_score = df_silhouette,
#                   CH_index        = df_CH,
#                   DB_index        = df_DB)
  
#   # Join them together
#   temp <- df_list %>%
#     map(as_tibble) %>%
#     map(~ .x %>%
#           pivot_longer(everything())) %>%
#     imap(~ rename(.x, "{.y}" := value)) %>%
#     reduce(left_join)
  
#   #Extract the best resolutions by metric
#   silh_res <- temp %>% arrange(desc(silhouette_score)) %>%
#     select(name) %>% rename(silhouette_score = name)
  
#   CH_res <- temp %>% arrange(desc(CH_index)) %>%
#     select(name) %>% rename(CH_index = name)
  
#   DB_res <- temp %>% arrange(DB_index) %>%
#     select(name) %>% rename(DB_index = name)
  
#   # joining ordered resolutions
#   result <- data.frame(silh_res, CH_res, DB_res)
  
#   return(result)
# }
