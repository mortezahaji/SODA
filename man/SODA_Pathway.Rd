\name{SODA_Pathway}
\alias{SODA_Pathway}

\title{
Ditribution of Human Hallmarks Over the Interpolated Coordinates
}
\description{
This function finds the list of the genes in human hallmarks from the Molecular Signatures Database (MSigDB) and displays the distribution of human hallmarks over the interpolated coordinates
}
\usage{
SODA_Pathway (Data,Pathway,Type=c("3D","2D","Interactive"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{  A data frame contains interpolated data, output of \code{\link{SODA_Kriging}}}
  \item{Pathway}{  The names of the hallmark in MSigDB in quotation. For example Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY".}
  \item{Type}{  Type of the plot could be 2D, 3D, or interactive using plotly package.}
}

\value{
  A 2D, 3D, or interactive plot
}
\references{
Hajihosseini, M., Amini, P., Voicu, D., Dinu, I., & Pyne, S. (2022). Geostatistical Modeling and Heterogeneity Analysis of Tumor Molecular Landscape. Cancers, 14(21), 5235. \cr
Amini, P., Hajihosseini, M., Pyne, S., & Dinu, I. (2022). Geographically Weighted Linear Combination Test for Gene Set Analysis of a Continuous Spatial Phenotype as applied to Intratumor Heterogeneity. bioRxiv. \cr
Liberzon, A., Subramanian, A., Pinchback, R., Thorvaldsdóttir, H., Tamayo, P., & Mesirov, J. P. (2011). Molecular signatures database (MSigDB) 3.0. Bioinformatics, 27(12), 1739-1740.

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
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="2D")

#For 3D Plot
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="3D")

#For Interactive Plot
SODA_Pathway(Data=Intrp_Data,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",Type="Interactive")

}

