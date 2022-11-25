\name{SODA_Plot}
\alias{SODA_Plot}

\title{
2D, 3D, and Interactive Gene Expression or Pathway Plots
}
\description{
This function generates 2D, 3D, and interactive plots based on the interpolated data
}
\usage{
SODA_Plot (Data,Var,X,Y,Type=c("3D","2D","Interactive"),Res=500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{  A data frame contains interpolated data, output of \code{\link{SODA_Kriging}}}
  \item{X}{  Name of the variable that indicate the longitude of the points in the interpolated data in quotation. For example X="Longitude".}
  \item{Y}{  Names of the variable that indicates the latitude of the points in the interpolated data in quotation. For example X="Latitude".}
  \item{Type}{ Type of the plot could be 2D, 3D, or interactive using plotly package.}
  \item{Res}{  Indicating resolution for 3D plots.}
}

\value{
  A 2D, 3D, or interactive plot
}
\references{
Hajihosseini, M., Amini, P., Voicu, D., Dinu, I., & Pyne, S. (2022). Geostatistical Modeling and Heterogeneity Analysis of Tumor Molecular Landscape. Cancers, 14(21), 5235. \cr
Amini, P., Hajihosseini, M., Pyne, S., & Dinu, I. (2022). Geographically Weighted Linear Combination Test for Gene Set Analysis of a Continuous Spatial Phenotype as applied to Intratumor Heterogeneity. bioRxiv.

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

#For 2D Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="2D")

#For 3D Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="3D")

#For Interactive Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="Interactive")

}

