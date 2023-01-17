# SODA
## Overview
Spatial Omics Data Analysis (SODA) is an R package designed to facilitate differentially expressed gene and heterogeneity analysis of tumor molecular landscape for spatial transcriptomics data 
by providing interactive visualizations that display the synthesize high-resolution and continuous gene expression landscapes of a given tumor sample
and integration of such landscapes to identify and map the enriched regions of the pathways of interest with identification of genes
that have a spatial differential expression at locations representing specific phenotypic contexts.
An important feature of SODA is the computation of spatial entropy measures for quantification and objective characterization of intratumor heterogeneity. 
The package is built around a shiny "gadget" to allow the exploration of the data with multiple plots in parallel and an interactive UI.

## Installation
```
install.packages("devtools")
devtools::install_github("mortezahaji/SODA")
require(SODA)

#Functions & Data
?SODA_Data
?SODA_Kriging
?SODA_Plot
?SODA_Wass
?SODA_Pathway
?SODA_Entrp

```
### SODA_Data
#### This dataset contains a subset of the Human Breast Cancer data we used in our recent publications [listed at the end of this page] downloaded from 10x Genomics online resource.
```
?SODA_Data
data(SODA_Data)
View(SODA_Data)

```
### SODA_Kriging (Data,X,Y,Expand,Margin)
#### This function generates the interpolated dataset using Ordinary Kriging method.
```
?SODA_Kriging
data(SODA_Data)
View(SODA_Data)
Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)
```

### SODA_Plot (Data,BM,X,Y,Type=c("3D","2D","Interactive"),Res=500)
#### This function generates 2D, 3D, and interactive plots based on the interpolated data.
```
?SODA_Plot
#Data
data(SODA_Data)
View(SODA_Data)

Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)

#For 2D Plot
SODA_Plot(Data=Intrp_Data,BM="S100A10",X="Longitude",Y="Latitude",Type="2D")

#For 3D Plot
SODA_Plot(Data=Intrp_Data,BM="S100A10",X="Longitude",Y="Latitude",Type="3D")

#For Interactive Plot
SODA_Plot(Data=Intrp_Data,BM="S100A10",X="Longitude",Y="Latitude",Type="Interactive")
```

### SODA_Wass(Data1,Data2,X,Y,BM,cutoff,Rad,Method=c("TS","OS"))
#### This function performes spatial differentially expressed gene (DEG) analysis between the cells located in the area with upper than the expressed threshold for a specific biomarker (local cells) and cells located in other locations of the cell (Non-local cells).
```
?SODA_Wass
#Data
data(SODA_Data)
View(SODA_Data)
Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)	

#Rad is one of the percentiles of the distance matrix between different cells
pts <- SpatialPoints(data.frame(Longitude=SODA_Data$Longitude,Latitude=SODA_Data$Latitude))
# longlat = FALSE returns Euclidean distance
euclidDist <- sp::spDists(pts,longlat = FALSE)
Distance=ecdf(euclidDist)
X=data.frame(Dist=euclidDist[upper.tri(euclidDist, diag = FALSE)])
Q=quantile(X[,1],names = TRUE,seq(0,1,0.05))
Q
#Q[2] is 55.33
	
#Based on the “TS” for the two-stage method
Result1=SODA_Wass(Data1=Intrp_Data,Data2=SODA_Data,X="Longitude",Y="Latitude",
			BM="S100A10",cutoff=1,Rad=55.33,Method="TS")

#Based on the “OS” for the one-stage method
Result1=SODA_Wass(Data1=Intrp_Data,Data2=SODA_Data,X="Longitude",Y="Latitude",
			BM="S100A10",cutoff=1,Rad=55.33,Method="OS")
```

### SODA_Pathway (Data,Pathway,Type=c("3D","2D","Interactive"),species,category)
#### This function finds the list of the genes in human hallmarks from the Molecular Signatures Database (MSigDB) and displays the distribution of human hallmarks over the interpolated coordinates.
```
?SODA_Pathway
#Data
data(SODA_Data)
View(SODA_Data)
Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)	

#For 2D Plot
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="2D",species = "Human", category = "H")

#For 3D Plot
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="3D",species = "Human", category = "H")

#For Interactive Plot
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="Interactive",species = "Human", category = "H")
```

### SODA_Entrp (Data,X,Y,BM=NULL,Pathway=NULL,cutoff,Part,species, category)
#### This function calculates Batty's Entropy index for various number of partitions for a biomarker or a geneset.
```
?SODA_Entrp
#Data
data(SODA_Data)
View(SODA_Data)
Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)	

#For a Biomarker
SODA_Entrp(Data=Intrp_Data,X="Longitude",Y="Latitude",
				BM="S100A10",Pathway=NULL,cutoff=1.5,Part=2:10)

#For a Gene Set
SODA_Entrp(Data=Intrp_Data,X="Longitude",Y="Latitude",
				BM=NULL,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",
				species = "Human", category = "H",cutoff=1.5,Part=2:10)
```
###SODA_GWLCT(Data,X,Y,ovrp,BM,kernel=c("Bisquare","Tricube"),method=c("adaptive","fixed"),bw,pthres,qthres,nbPermutations)
####


## Useful Links
[Geostatistical Modeling and Heterogeneity Analysis of Tumor Molecular Landscape](https://www.mdpi.com/2072-6694/14/21/5235)\
[Geographically Weighted Linear Combination Test for Gene Set Analysis of a Continuous Spatial Phenotype as applied to Intratumor Heterogeneity](https://www.biorxiv.org/content/10.1101/2022.10.09.511477v1.abstract)
