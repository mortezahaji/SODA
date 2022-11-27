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

#Functions & Data
?SODA_Data
?SODA_Kriging
?SODA_Plot
?SODA_Wass
?SODA_Pathway

```
### SODA_Data
```
?SODA_Data
data(SODA_Data)
View(SODA_Data)

```
### SODA_Kriging (Data,X,Y,Expand,Margin)
```
?SODA_Kriging
data(SODA_Data)
View(SODA_Data)
Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)
```

### SODA_Plot (Data,Var,X,Y,Type=c("3D","2D","Interactive"),Res=500)
```
?SODA_Plot
#Data
data(SODA_Data)
View(SODA_Data)

Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)

#For 2D Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="2D")

#For 3D Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="3D")

#For Interactive Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="Interactive")
```

### SODA_Wass(Data1,Data2,X,Y,BM,cutoff,Rad,Method=c("TS","OS"))
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
Rad=Q[2]
	
#Based on the “TS” for the two-stage method
Result1=SODA_Wass(Data1=Intrp_Data,Data2=SODA_Data,X="Longitude",Y="Latitude",
			BM="S100A10",cutoff=1,Rad=55.33,Method="TS")

#Based on the “OS” for the one-stage method
Result1=SODA_Wass(Data1=Intrp_Data,Data2=SODA_Data,X="Longitude",Y="Latitude",
			BM="S100A10",cutoff=1,Rad=55.33,Method="OS")
```

### SODA_Pathway (Data,Pathway,Type=c("3D","2D","Interactive"))
```
?SODA_Pathway
#Data
data(SODA_Data)
View(SODA_Data)
Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)	

#For 2D Plot
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="2D")

#For 3D Plot
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="3D")

#For Interactive Plot
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="Interactive")
```
