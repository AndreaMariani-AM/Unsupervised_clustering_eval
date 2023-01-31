# Workflow

This R markdown file contains a general workflow with some explanations behind the set of functions created in the [metrics.R](https://github.com/AndreaMariani-AM/Unsupervised_clustering_eval/blob/main/scripts/metrics.R) file.  
The main idea was to make the steps before clustering more rigorous, and assign some values to the clustering params that could help me out
in deciding which visual representation is the one that can retain the highest information of my dataset.  

## Steps

1. Perform classical QC for single cell  
2. Perform PCA  						
3. Decide which hyperparameter needs tuning  
4. Compute the Score/Indices for that hyperparameter  
5. Inspect the plots and decides a value for the hyperparameter  

## Decide which hyperparameter need tuning

1. N. of PCs  
2. N. of neighbors when constructing the graph  
3. Distance Metric  
4. Community detection algorithm  

## Pseudo code (for now)

```{r}
read_in_10X_data <- function(x) {}

perform_QC_steps

hyperameter <- "value/char"

compute_metrics <- function(x) {}

plot_metrics <- function(x) {}
```