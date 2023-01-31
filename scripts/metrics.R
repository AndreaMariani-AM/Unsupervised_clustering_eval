###---------------------------------------------------------------------###
###                         INTRODUCTION 								###
###---------------------------------------------------------------------###

# Implementation of silhoutte scores, Calinski-Harabasz Index and Davies-Bouldin Index
# For now the output is always splitted based on the resolution for clustering
# The next thing to do is to evaluate which hyperparameters need to be fine tuned


###---------------------------------------------------------------------###
###                       Function for the Metrics 						###
###---------------------------------------------------------------------###

# Retrieve PCA cell embeddings used later to as a space for clustering
count_matrix <- seurat@reductions$pca@cell.embeddings
# Retrieve resolutions, mine starts from 0, hence the 2:lenght(resolutions) when computing
# the scores
resolutions_to_evaluate <- names(seurat@meta.data) %>% str_subset(pattern = "^SCT")


# Silhouette Score

Silhouette_score <- map(resolutions_to_evaluate[2:length(resolutions_to_evaluate)], 
                        ~ mean(cluster::silhouette(as.numeric(seurat@meta.data[[.x]]), dist(count_matrix))[,3])) %>%
  set_names(resolutions_to_evaluate[2:length(resolutions_to_evaluate)])

# Calinski Harabasz Index

CH_index <- map(resolutions_to_evaluate[2:length(resolutions_to_evaluate)],
                ~ fpc::calinhara(x = count_matrix, clustering = as.numeric(seurat@meta.data[[.x]]))) %>% 
  set_names(resolutions_to_evaluate[2:length(resolutions_to_evaluate)])

# Davies Bouldin Index

DB_index <- map(resolutions_to_evaluate[2:length(resolutions_to_evaluate)],
                ~ index.DB(x = count_matrix, cl = (as.numeric(seurat@meta.data[[.x]])))) %>% 
  map(pluck, 1) %>% 
  set_names(resolutions_to_evaluate[2:length(resolutions_to_evaluate)])

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