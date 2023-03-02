# INTRODUCTION (UNDER DEVELOPMENT)

This repo contains some functions implemented in R that can be used to evaluate the performance of clustering when supervised labels are not available (unsupervised clustering).  
A classical example in biology is clustering of cells in a scRNA-seq experiments.  
A set of Hyperparameters is decided a priori and then is evaluated through the metrics reported here.  
This will hopefully retrieve optimal `values` for our clustering problem.  

A brief workflow is reported here in the same repo

## HYPERPARAMETERS

Here a brief discussion on which hyperparameters i think needs evaluation and tuning.  
I'll eventually add more, if i have reasons to believe others need to be included.  

1. The Number of Neighbors that are used when constructing the graph
2. Minimum distance between cells
3. Distance Metrics used for computing neighbors
4. Algorithm for community detection

## METRICS

The aim is to use some evaluation metrics for clustering, and guide the decision on the number of clusters that can best explain the data.  
When clustering data, we should strive for clusters that are:
  - Dense  
  - Well-separated  
  - Non-overlapping  

Three key metrics are:
  - Silhouette Score  $$\Large s(i) = \frac{b(i) - a(i)}{max(a(i),b(i))}$$
  - Calinski-Harabasz Index  $$\Large CH = \frac{Tr(B)}{Tr(W)} \cdot \frac{n - k}{k - 1}$$
  - Davies-Bouldin Index $$\Large DB = \frac{1}{k} \sum_{i =1}^{k} max_{i \not = j} Sim_{ij}$$ where $$\large Sim_{ij} = \frac{s_i + s_j}{d_{ij}}$$
  
  ## RUNNING TIME
  
As a small note, i've tested if running time of that shitty (i mean *cute*) function can scale up at least decently, when considering an increasingly higher number of cells. Turns out that my biggest fear, aka scaling exponentially, seems like not true. Though slow af, computing time doesn't seem to grow exponentially, at least for the number of cells.  
A whole different story are the hyperparameters, for which i kinda know for sure it's gonna be exponentially if, let's say, one decide to use 2 insted of 1 Algorithm for community detection and so on. As a general rule of thumb, this should increase for each additional `value` an hyperparameter can have.


<p align="center">
  <img width="900" height="750" src="https://github.com/AndreaMariani-AM/Unsupervised_clustering_eval/blob/main/images/running_time.jpeg">
</p>
