\name{SODA_Kriging}
\alias{SODA_Kriging}
\title{
Interpolation Using Ordinary Kriging
}
\description{
This function makes the interpolated dataset using Ordinary Kriging method
}
\usage{
SODA_Kriging (Data,X,Y,Expand,Margin)
}

\arguments{
  \item{Data}{ A dataframe conaints samples in rows, genes in coloums with longitude and latitude of cells as extra columns.}
  \item{X}{Name of the variable that indicate the longitude of the cells in quotation. For example X="Longitude".}
  \item{Y}{Names of the variable that indicates the latitude of the cells in quotation. For example X="Latitude".}
  \item{Expand}{An integer indicating how much expansion is needed from the original coordinates, see example below.}
  \item{Margin}{An integer indicating how much margin need to be added to the X and Y, see example below.}
}
\details{
Expand and Margin will be used in the "expand.grid" function. For example: \cr
  expand.grid(seq(min(X-Margin),max(X+Margin),l=Expand),seq(Y-Margin),max(Y+Margin),l=Expand))
}
\value{
A dataframe contains interpolated data with coordinates of the points.
}
\references{
Hajihosseini, M., Amini, P., Voicu, D., Dinu, I., & Pyne, S. (2022). Geostatistical Modeling and Heterogeneity Analysis of Tumor Molecular Landscape. Cancers, 14(21), 5235. \cr
Amini, P., Hajihosseini, M., Pyne, S., & Dinu, I. (2022). Geographically Weighted Linear Combination Test for Gene Set Analysis of a Continuous Spatial Phenotype as applied to Intratumor Heterogeneity. bioRxiv. \cr
Goovaerts, P. (2000). Geostatistical approaches for incorporating elevation into the spatial interpolation of rainfall. Journal of hydrology, 228(1-2), 113-129.
}
\author{
Morteza Hajihosseini \cr
Payam Amini \cr
Dan Voicu \cr
Irina Dinu \cr
Saumyadipta Pyne
}

\examples{
data(SODA_Data)
Intrp_Data=SODA_Kriging (Data=SODA_Data,X="Longitude",Y="Latitude",Expand=30,Margin=10)
}

