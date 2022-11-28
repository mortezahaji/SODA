\name{SODA_Entrp}
\alias{SODA_Entrp}

\title{
This function calculates Batty's Entropy index for various number of partitions for a biomarker or a geneset.
}
\description{
This function calculates Batty's Entropy index and generates a table and two figures
}
\usage{
SODA_Entrp (Data,X,Y,BM=NULL,Pathway=NULL,cutoff,Part)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{  A data frame contains interpolated data, output of \code{\link{SODA_Kriging}}}
  \item{X}{Name of the variable that indicate the longitude of the cells in quotation. For example X="Longitude".}
  \item{Y}{Names of the variable that indicates the latitude of the cells in quotation. For example X="Latitude".}
  \item{BM}{ The name of a specific biomarker. It is NULL when "Pathway"" specified.}
  \item{Pathway}{  The names of the hallmark in MSigDB in quotation. It is NULL when "BM"" specified. For example Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY".}
  \item{cutoff}{  The threshold for the biomarker or geneset z-score that defines over-expressed and under-expressed areas.}
  \item{Part}{ A number or a sequence of numbers indicating the number of partitions for Batty's Entropy index.}
}

\value{
  A table with four columns: \cr
  1. The number of partitions\cr
  2. The average of entropy index over permutations\cr
  3. The standard deviation of entropy index over permutations\cr
  4. The median of entropy index over permutations \cr
  Two graphs:\cr
  1. Showing the distribution of the over-expressed and under-expressed points\cr
  2. Showing the results of Batty's Entropy index for various number of partitions over 1000 permutaions
}
\references{
Hajihosseini, M., Amini, P., Voicu, D., Dinu, I., & Pyne, S. (2022). Geostatistical Modeling and Heterogeneity Analysis of Tumor Molecular Landscape. Cancers, 14(21), 5235. \cr
Amini, P., Hajihosseini, M., Pyne, S., & Dinu, I. (2022). Geographically Weighted Linear Combination Test for Gene Set Analysis of a Continuous Spatial Phenotype as applied to Intratumor Heterogeneity. bioRxiv. \cr
Batty, M. (1976). Entropy in spatial aggregation. Geographical Analysis, 8(1), 1-21.

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

#For a Biomarker
SODA_Entrp(Data=Intrp_Data,X="Longitude",Y="Latitude",
				BM="S100A10",Pathway=NULL,cutoff=1.5,Part=2:10)

#For a Gene Set
SODA_Entrp(Data=Intrp_Data,X="Longitude",Y="Latitude",
				BM=NULL,Pathway="HALLMARK_ESTROGEN_RESPONSE_EARLY",cutoff=1.5,Part=2:10)

}

