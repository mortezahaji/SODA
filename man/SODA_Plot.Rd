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
  \item{Data}{  A data frame contains interpolated data, output of \code{SODA_Kriging}}
  \item{X}{  Name of the variable that indicate the longitude of the points in the interpolated data in quotation. For example X="Longitude".}
  \item{Y}{  Names of the variable that indicates the latitude of the points in the interpolated data in quotation. For example X="Latitude".}
  \item{Type}{ Type of the plot could be 2D, 3D, or interactive using plotly package.}
  \item{Res}{  Indicating resolution for 3D plots.}
}

\value{
  A 2D, 3D, or interactive plot
}
\references{
Hajihosseini M, Amini P, Voicu D, Dinu I, Pyne S. Geostatistical Modeling and Heterogeneity Analysis of Tumor Molecular Landscape. Cancers. 2022; 14(21):5235. https://doi.org/10.3390/cancers14215235

}
\author{
Morteza Hajihosseini
}

\examples{
#For 2D Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="2D")

#For 3D Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="3D")

#For Interactive Plot
SODA_Plot(Data=Intrp_Data,Var="S100A10",X="Longitude",Y="Latitude",Type="Interactive")

}

