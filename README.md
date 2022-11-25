# SODA
## Overview
SODA is an R package designed to facilitate differentially expressed gene and heterogeneity analysis of tumor molecular landscape for spatial transcriptomics data 
by providing interactive visualizations that display the synthesize high-resolution and continuous gene expression landscapes of a given tumor sample
and integration of such landscapes to identify and map the enriched regions of the pathways of interest with identification of genes
that have a spatial differential expression at locations representing specific phenotypic contexts.
An important feature of SODA is the computation of spatial entropy measures for quantification and objective characterization of intratumor heterogeneity. 
The package is built around a shiny "gadget" to allow the exploration of the data with multiple plots in parallel and an interactive UI.

## Installation
```
# install.packages("devtools")
devtools::install_github("mortezahaji/SODA")
require(SODA)

#Functions
?SODA_Kriging
?SODA_Plot
?SODA_Wass
?SODA_Pathway
?SODA_Data

```
