SODA_Pathwayovrp<-function(Data,X,Y,ovrp,species,category)
{


suppressPackageStartupMessages({
require(msigdbr)
require(ExperimentHub)
require(GSEABase)
})


Lo=Data[,X]
La=Data[,Y]
DATA=Data[,!(names(Data) %in% c(X,Y))]
DATA=t(Data)

Gnames=rownames(DATA)
#Retrieve H (hallmark) gene set
msigdbr_df <- msigdbr(species = species, category = category)
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

N=c()
GS=matrix(NA,length(Gnames),length(pathwaysH)+1)
GS[,1]=Gnames
for(j in 2:(length(pathwaysH)+1))
{
set=unique(unlist(pathwaysH[[j-1]]))
N=c(N,length(unique(set)))
	for(i in 1:length(Gnames))
	{
		GS[i,j]=ifelse(GS[i,1] %in% set,1,0)
	}
}
colnames(GS)=c("Genes",names(pathwaysH))
rownames(GS)=Gnames
GS2=data.frame(GS[,-1])
GS2=GS2 %>% mutate_if(is.character, as.numeric)
CS=apply(GS2,2,sum)
Pro=CS/N
OvrpNam=names(which(Pro>=ovrp))

Names=colnames(GS2)[colnames(GS2) %in% OvrpNam]

return(Names)
}