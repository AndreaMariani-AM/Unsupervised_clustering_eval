# Workflow

This R markdown file contains a general workflow with some explanations behind the set of functions created in the [metrics.R](https://github.com/AndreaMariani-AM/Unsupervised_clustering_eval/blob/main/scripts/metrics.R) file.  
The main idea was to make the steps before clustering a little bit more rigorous, and assign some values to the clustering params that could help me out
in deciding which visual representation is the one that can retain the highest information of my dataset.  
The chosen hyperparameters have to be evaluated at the same to guarantee some sort of "best" possible combinations.

## Steps
Steps that needs to be performed:  
1. Perform QC, normalization and PCA 
2. Choose a certain Number of PCA (elbow method or others), then keep it fixed  
3. Decide which hyperparameters need tuning

## Hyperparameters

1. Number of Neighbors when constructing the graph (k.param)  
2. Min dist neighbors  
3. Distance Metrics for annoy(euclidean, cosine, manhattan, hamming)  
4. Algorithm for finding clusters( OG Louvain, Louvain refinement, SLM, Leiden)  

Hyperparameters `1`,`2`, `3`, `4` are evaluated at the same time

## Pseudocode for basic analysis

#### Read in expression matrix
```{r}
read_in_10X_data <- Read10X()
```

#### QC/Normalization/PCA

This steps require a litte bit of work. Ideally you shoul find a QC routine that fits you needs. For normalization influence on clustering metrics
look at the end of the document.

```{r}
perform_QC_steps
perform_normalization
perform_PCA
```

#### Tune Hyperparameters

This section employes the first function of [metrics.R](https://github.com/AndreaMariani-AM/Unsupervised_clustering_eval/blob/main/scripts/metrics.R).  
This function iteratively performs different combinations of hyperparameters specified by the user (look above on the ones implemented as of now).  
To assess running time based on the number of cells, see the [README file](https://github.com/AndreaMariani-AM/Unsupervised_clustering_eval). I haven't tested this yet, mainly due to a lack of time, but i expect an exponentiall growth in running time as the hyperparameters grow.

```{r}
metrics <- test_clust()
```

#### Retrieve Metrics Plots

This section uses the second function to actully plot the metrics. The output ig going to be a PDF, which name is derived by the combination of hyperparameters tested, that contains the three metrics.

```{r}
plot_clusters
```

## TO DO

1) Find a way to return to the user the `TOP` resolutions for each combination. This could be useful to avoid too much visual inspection of metrics plots(tho always recommended).
2) As of now, i don't know how much normalization/QC steps can influence clustering metrics. I guess that Qc isn't gonna be a big of a deal since i assume the user is going to carefully assess quality metrics `BEFORE` even thinking about clustering and other downstream analysis. Normalization, on the other end, i think is gonna cause some troubles. 
