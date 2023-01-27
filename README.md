# INTRODUCTION

This repo contains a script implemented in R that can be used to evaluate the performance of clustering when supervised labels are not available (unsupervised clustering).  
A classical example in biology is clustering of cells in a scRNA-seq experiments.  


## METRICS

The aim is to use some evaluation metrics for clustering, and guide the decision on the number of clusters that can best explain the data.  
When clustering data, we should strive for clusters that are:
  - Dense  
  - Well-separated  
  - Non-overlapping  

Three key metrics are:
  - Silhoutte Score  
  - Calinski-Harabasz Index  
  - Davies-Bouldin Index
