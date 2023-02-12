# Workflow

This R markdown file contains a general workflow with some explanations behind the set of functions created in the [metrics.R](https://github.com/AndreaMariani-AM/Unsupervised_clustering_eval/blob/main/scripts/metrics.R) file.  
The main idea was to make the steps before clustering a little bit more rigorous, and assign some values to the clustering params that could help me out
in deciding which visual representation is the one that can retain the highest information of my dataset.  
Evaluating all the possible combinations is going to scale exponentially and probably it's not even worth it. In this case i'll approach the problem by diving the hyperparameter tuning into three main steps, where one o two parameters are being evaluated and the held fixed when an optimal has been found.

## Steps
Steps that needs to be performed:  
1. Perform QC, normalization (SCT) and PCA 
2. Choose a certain Number of PCA (elbow method or others), then keep it fixed  
3. Decide which hyperparameters need tuning

## Hyperparameters

1. Number of Neighbors when constructing the graph (k.param)  
2. Min dist neighbors  
3. Distance Metrics for annoy(euclidean, cosine, manhattan, hamming)  
4. Algorithm for finding clusters( OG Louvain, Louvain refinement, SLM, Leiden)  

Hyperparametr `1` and `2` are going to be evaluated at the same time when performing UMAP. Param `3` is evaluated on the neighbors finding and lastly, param `4` is tested when performing clustering.

## Pseudo code (for now)

```{r}
read_in_10X_data <- function(x) {}

perform_QC_steps

hyperameter <- "value/char"

test_clust <- function(x) {}

top_candidates_metrics

plot_metrics <- function(x) {}
```
