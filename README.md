<img width="468" height="56" alt="image" src="https://github.com/user-attachments/assets/dc518332-410f-4b58-ba17-7e7e1ec567ca" /># Sediment Metagenomic Functional Network Analysis

This repository contains the analysis scripts used for the metagenomic functional network analyses described in our manuscript. The workflow focuses on the construction of MAG–KO bipartite networks, calculation of network topological properties, temporal functional changes, and the relationships between environmental variables, redox conditions, and microbial functional organization.

## Overview

The analyses in this repository include:

1. Construction of MAG–KO bipartite networks
2. Calculation of network topological properties
3. Generation and analysis of GraphML network files
4. Non-metric multidimensional scaling analysis
5. Functional changes across different temporal phases
6. Correlations between environmental variables and functional genes
7. Generalized additive models for temporal trends
8. Piecewise linear regression analysis
9. Variance inflation factor analysis
10. Relative importance analysis
11. Causal machine learning analysis

## Repository structure

```text
.
├── Calculate topological propertiess/
├── Causual machine learning/
├── Construct Bipartite Network/
├── Correlations between env and functions/
├── Folder changes in phases/
├── Generalized additive models/
├── Graphml networks/
├── Non-metric multidimensional scaling/
├── Piecewise linear regressions/
├── Relative importance/
└── Variance inflation factors/
````

## Folder description

### `Construct Bipartite Network/`

These networks describe the associations between microbial bins and functional genes.

Typical outputs include edge lists, adjacency matrices, and network files used for downstream analyses.

### `Graphml networks/`

This folder contains GraphML files generated from the MAG–KO bipartite networks. These files can be used for network visualization and for calculating network-level properties.

### `Calculate topological parameters/`

This folder contains scripts used to calculate network topological properties, including:

* number of nodes
* number of edges
* network size
* average degree
* average path length
* clustering coefficient
* network density
* modularity
* betweenness centrality

The output files summarize the structural properties of each network.

### `Folder changes in phases/`

This folder contains scripts used to analyze functional changes across different temporal phases. These analyses were used to compare the temporal dynamics of functional genes and microbial functional potential among sedimentary periods.

### `Non-metric multidimensional scaling/`

This folder contains scripts for non-metric multidimensional scaling analysis. NMDS was used to evaluate differences in functional gene composition among samples, lakes, or temporal phases.

### `Correlations between env and functions/`

This folder contains scripts for analyzing relationships between environmental variables and functional genes. These analyses were used to identify potential environmental drivers of functional gene variation.

### `Generalized additive models/`

This folder contains scripts for generalized additive models. GAMs were used to evaluate nonlinear temporal trends in microbial functional profiles and network properties.

### `Piecewise linear regressions/`

This folder contains scripts for piecewise linear regression analysis. These analyses were used to detect potential thresholds or breakpoint responses between environmental variables and network properties.

### `Variance inflation factors/`

This folder contains scripts for calculating variance inflation factors. VIF analysis was used to evaluate multicollinearity among environmental variables before downstream statistical analyses.

### `Relative importance/`

This folder contains scripts for estimating the relative importance of different environmental or ecological drivers. These analyses help identify which variables contribute most strongly to the observed variation in network or functional properties.

### `Causual machine learning/`

This folder contains scripts for causal machine learning analysis. These analyses were used to evaluate potential causal relationships among nutrient conditions, climate variables, redox conditions, and microbial network properties.

## Input data

Large raw sequencing files are not included in this repository. Raw sequencing data should be accessed through the public database accession number provided in the manuscript.

## Output files

The scripts generate outputs including:

* MAG–KO bipartite network files
* GraphML network files
* network topology summary tables
* functional gene abundance summaries
* NMDS ordination results
* correlation analysis results
* GAM fitting results
* piecewise regression results
* VIF summary tables
* relative importance results
* causal machine learning outputs
* figures and supplementary results

## Software

The analyses were mainly performed in R. The major R packages used include:
tidyverse
igraph
ggraph
tidygraph
vegan
mgcv
segmented
car
relaimpo
ggplot2
...

All scripts are organized according to the main analytical steps in the manuscript. File paths may need to be modified according to the local directory structure before running the scripts.

For reproducibility, users are encouraged to run the scripts in the order listed above and keep the input and output file names consistent with those used in the scripts.

## Citation

If you use the scripts or workflow from this repository, please cite the associated manuscript. The full citation and DOI will be provided after the manuscript is accepted:

Siwei, Y., Xiaofeng, C., Yushuai, W., Deen, F., Tong, C., Pengfei, C., Shaopeng, G., Hongwei, Y., Weixiao, Q., Jianfeng, P., Huijuan, L., Jiuhui, Q., 2026. Long-term deoxygenation reorganizes genome-resolved microbial functional networks in anoxic lake sediments. (Submitted)

## Contact

If you have any questions, please contact:

Dr. Siwei Yu  
Email: yusiwei_07@163.com
