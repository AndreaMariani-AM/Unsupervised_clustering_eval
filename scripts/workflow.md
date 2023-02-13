# Workflow

This R markdown file contains a general workflow with some explanations behind the set of functions created in the [metrics.R](https://github.com/AndreaMariani-AM/Unsupervised_clustering_eval/blob/main/scripts/metrics.R) file.  
The main idea was to make the steps before clustering a little bit more rigorous, and assign some values to the clustering params that could help me out
in deciding which visual representation is the one that can retain the highest information of my dataset.  
The chosen hyperparameters have to be evaluated at the same to guarantee some sort of "best" possible combinations.

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

Hyperparameters `1`,`2`, `3`, `4` are evaluated at the same time

## Pseudo code (for now)

```{r}
read_in_10X_data <- function(x) {}

perform_QC_steps

hyperameter <- "value/char"

test_clust <- function(x) {}

top_candidates_metrics

plot_metrics <- function(x) {}
```
