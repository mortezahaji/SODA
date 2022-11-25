\name{SODA_Wass}
\alias{SODA_Wass}

\title{
Spatial Differentially Expressed Gene Analysis Using 2-Wasserstein Distance Test
}
\description{
This function performes spatial differentially expressed gene test between the cells located in the area with upper than the expressed threshold for a specific biomarker (local cells) and cells located in other locations of the cell (Non-local cells).
}
\usage{
SODA_Wass(Data1,Data2,X=,Y=,BM,cutoff,Rad,Method=c("TS","OS"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data1}{  A data frame contains interpolated data, output of \code{\link{SODA_Kriging}}}
  \item{Data2}{  A data frame contains gene expression data with samples in the rows and genes in the columns.}
  \item{X}{Name of the variable that indicate the longitude of the cells in quotation. For example X="Longitude".}
  \item{Y}{Names of the variable that indicates the latitude of the cells in quotation. For example X="Latitude".}
  \item{BM}{ The name of a specific biomarker}
  \item{cutoff}{  The threshold for the biomarker z-score that defines local and non-local cells.}
  \item{Rad}{ The radious of the distance from the biomarker that defines the local and non-local areas.}
  \item{Method}{ The 2-Wasserstein distance test method, more details: \code{\link{wasserstein.sc}}}
}

\value{
  A data frame with four column: \cr
  N_High is the number of local cells (area with the biomarker z-score over the threshold) \cr
  N_Low  is the number of non-local cells (area with the biomarker z-score under the threshold) \cr
  Pvalue is the p_value of the test, examining if each gene was differentially expressed between two areas \cr
  Adjusted_Pvalue is an ordered adjusted p_value according to the method of Benjamini-Hochberg.
}
\references{
Hajihosseini, M., Amini, P., Voicu, D., Dinu, I., & Pyne, S. (2022). Geostatistical Modeling and Heterogeneity Analysis of Tumor Molecular Landscape. Cancers, 14(21), 5235. \cr
Amini, P., Hajihosseini, M., Pyne, S., & Dinu, I. (2022). Geographically Weighted Linear Combination Test for Gene Set Analysis of a Continuous Spatial Phenotype as applied to Intratumor Heterogeneity. bioRxiv. \cr
Ramdas, A., García Trillos, N., & Cuturi, M. (2017). On wasserstein two-sample testing and related families of nonparametric tests. Entropy, 19(2), 47.
}
\author{
Morteza Hajihosseini \cr
Payam Amini \cr
Dan Voicu \cr
Irina Dinu \cr
Saumyadipta Pyne
}

\examples{
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

}

