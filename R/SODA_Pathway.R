

SODA_Pathway=function(Data,Pathway,Type=c("3D","2D","Interactive"),species,category)
{

suppressPackageStartupMessages({
require(msigdbr)
require(ExperimentHub)
require(GSEABase)
})


Gnames=colnames(Data)

#Retrieve human H (hallmark) gene set
msigdbr_df <- msigdbr(species = species, category = category)
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

set=unique(unlist(pathwaysH[[Pathway]]))
Nset=length(which(set %in% Gnames))

Pathway_Score=data.frame(rowSums(Data[colnames(Data) %in% set]/sqrt(length(Nset))
,na.rm=TRUE))
names(Pathway_Score)=Pathway

Data_Plot=data.frame(Data,Pathway_Score)

if(Type=="3D"){
SODA_Plot(Data=Data_Plot,BM=Pathway,X="Longitude",Y="Latitude",Type="3D")

}else if(Type=="2D"){
SODA_Plot(Data=Data_Plot,BM=Pathway,X="Longitude",Y="Latitude",Type="2D")

}else if(Type=="Interactive"){
SODA_Plot(Data=Data_Plot,BM=Pathway,X="Longitude",Y="Latitude",Type="Interactive")
}

}
