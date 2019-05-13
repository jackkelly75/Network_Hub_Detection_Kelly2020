# Network_Hub_Iter
Permutation test to detect important hub genes in networks

Weighted Gene Co-expression Network Analysis (WGCNA) (Langfelder and Horvath, 2008) and other network analysis tools are becoming more and more prevalent in network biology.

This function takes a user created adjacency matrix, and a vector of genes with assigned modules, and calculated various hub scores to identify important genes.
To calculate module membership the expression data is required 

**Scores**
* Betweenness Centrality
* Closeness Centrality
* Kleinberg's Hubscore
* PageRank
* Module Membership
* Edge betweenness

**Running**

An adjacency matrix is required to analyse. An adjacency matrix can be created using any preffered method. One popular way when dealing with a gene expression matrix is the adjacency() function in WGCNA.

adjacency.Rdata contains an example of a subset of an adjacency gene matrix produced from healthy controls.

The moduleColors variable is a vector containing all nodes in the network (genes) with there names being the module they have been assigned.

moduleColors.Rdata contains a vector of the genes in the adjacency matrix provided above, with there assigned module after hierachical clustering of the full dataset.


The provided Rdata files should give an idea of how data should be layed out, and can be used with the function to show the hub genes within these modules.


**Output**



**Ref.**

Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008 Dec 29; 9:559.
