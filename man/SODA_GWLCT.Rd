\name{SODA_GWLCT}
\alias{SODA_GWLCT}

\title{
Geographically Weighted Linear Combination Test for Gene Set Analysis of a Continuous Spatial Phenotype
}
\description{
This function performs the geographically weighted linear combination test at the pathway level
}
\usage{
SODA_GWLCT(Data,X,Y,ovrp,BM,kernel=c("Bisquare","Tricube"),method=c("fixed","adaptive"),bw,pthres,qthres,nbPermutations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{ A data frame contains interpolated data, output of \code{\link{SODA_Kriging}}}
  \item{X} { Name of the variable that indicate the longitude of the points in quotation. For example X="Longitude".}
  \item{Y} { Names of the variable that indicates the latitude of the points in quotation. For example X="Latitude".}
  \item{ovrp} { A value between 0 and 1 that shows the required overlap between MSigDB human hallmark gene sets and genes in the original data set. for example ovrp=0.20}
  \item{BM} { The name of a specific biomarker.}
  \item{kernel} { Indicates the Kernel function including Bisquare and Tricube , more details are available in reference section at Amini P, et al.}
  \item{method} {Indicates the Kernel function including fixed and adaptive, more details are available in reference section at Amini P, et al.}
  \item{bw} { A value indicating the bandwidth, more details are available in reference section at Amini P, et al. }
  \item{pthres} { The threshold of interest for the p_value.}
  \item{qthres} { The threshold of interest for the q_value.}
  \item{nbPermutations} { The number of permutations.}
}
\details{
This function tests all human hallmark gene sets available in MSigDB that pass the overlap threshold.
}

\value{
A list that contains global LCT and GWLCT results.
An interactive plot showing the number of significant pathways at each location.
}
\references{
Amini, P., Hajihosseini, M., Pyne, S., & Dinu, I. (2022). Geographically Weighted Linear Combination Test for Gene Set Analysis of a Continuous Spatial Phenotype as applied to Intratumor Heterogeneity. bioRxiv. \cr
Hajihosseini, M., Amini, P., Voicu, D., Dinu, I., & Pyne, S. (2022). Geostatistical Modeling and Heterogeneity Analysis of Tumor Molecular Landscape. Cancers, 14(21), 5235.
}
\author{
Payam Amini \cr
Morteza Hajihosseini \cr
Irina Dinu \cr
Saumyadipta Pyne
}

\examples{
#Data
data(SODA_Data)
View(SODA_Data)
Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=35,Margin=10)

Result=SODA_GWLCT(Data=Intrp_Data,X="Longitude",Y="Latitude",Ovrp=0.20,
		BM="COL10A1",kernel="Bisquare",method="adaptive",
		bw=7,pthres=0.05,qthres=0.001,nbPermutations=1000)
}

