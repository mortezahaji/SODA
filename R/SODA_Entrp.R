
SODA_Entrp=function(Data,X,Y,BM=NULL,Pathway=NULL,cutoff,Part)
{

suppressPackageStartupMessages({
require(SpatEntropy)
require(spatstat)
require(ggplot2)
require(doBy)
require(msigdbr)
require(ExperimentHub)
require(GSEABase)
require(ggpubr)
})


if(!is.null(BM) & is.null(Pathway))
{
Data$DG=ifelse(Data[,BM]>=cutoff,1,0)
Data$DG=factor(Data$DG,levels=c(0,1),
			labels=c(paste0("z_",BM,"<",cutoff),paste0("z_",BM,">=",cutoff)))


P1=ggplot(Data, aes(Data[,X],Data[,Y],color=DG)) +
    		geom_point(size=2) +
    		xlab("Longitude")+ ylab("Latitude")+
		scale_color_manual(values=c("gray","red"))+
		ggtitle(paste0("Distribution of Over-Expressed Points for"," ",BM))+
		labs(color="Biomarker Z-score")+theme_classic()


Data2=ppp(x=Data[,X],y=Data[,Y],
		window=owin(c(min(Data[,X]),max(Data[,X])),
					c(min(Data[,Y]),max(Data[,Y]))),
		marks=as.numeric(Data$DG))

BE=list()
i=min(Part)
while(i<=max(Part))
{
mat=matrix(NA,1000,2)
colnames(mat)=c("iter","Batty_Entropy")
	for(j in 1:1000)
		{
			mat[j,1]=j
			mat[j,2]=round(batty(Data2,partition=i)$rel.batty,3)
		}
BE[[i-1]]=mat
i=i+1
}

Data3=as.data.frame(rbind(BE[[1]],BE[[2]],BE[[3]],BE[[4]],BE[[5]],BE[[6]],
			BE[[7]],BE[[8]],BE[[9]]))
Data3$Partitions=rep(Part,each=1000)
Data3$Partitions=as.numeric(Data3$Partitions)

P2=ggplot(Data3,aes(x=factor(Partitions), y=Batty_Entropy,col=factor(Partitions),
			group=factor(Partitions))) + 
	geom_boxplot() + stat_summary(inherit.aes = FALSE,
	aes(x=factor(Partitions), y=Batty_Entropy,
		group=1),fun.y=median,geom="line",size=1)+
		theme_classic()+
	ggtitle(paste0("Batty's Entropy for"," ",BM))+
	labs(x ="Number of Partitions",color="Partitions", y ="Batty's Entropy")+
	theme(legend.position="bottom")

Entrp=round(summaryBy(Batty_Entropy ~ Partitions, data = Data3, 
          FUN = list(mean,sd,median)),3)

print(ggarrange(P1,P2,labels = c("A", "B")))
return(Entrp)

}else if (!is.null(Pathway) & is.null(BM)){

Gnames=colnames(Data)
#Retrieve human H (hallmark) gene set
msigdbr_df <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

set=unique(unlist(pathwaysH[[Pathway]]))
Nset=length(which(set %in% Gnames))

Pathway_Score=data.frame(rowSums(Data[colnames(Data) %in% set]/sqrt(length(Nset))
,na.rm=TRUE))
names(Pathway_Score)=Pathway
Data=data.frame(Data,Pathway_Score)

Data$DP=ifelse(Data[,Pathway]>=cutoff,1,0)
Data$DP=factor(Data$DP,levels=c(0,1),
			labels=c(paste0("z_","Pathway","<",cutoff),paste0("z_","Pathway",">=",cutoff)))

P3=ggplot(Data, aes(Data[,X],Data[,Y],color=DP)) +
    		geom_point(size=2) +
    		xlab("Longitude")+ ylab("Latitude")+
		scale_color_manual(values=c("gray","red"))+
		ggtitle(paste0("Distribution of Over-Expressed Points for"," ",Pathway))+
		labs(color="Pathway Z-score")+theme_classic()


Data2=ppp(x=Data[,X],y=Data[,Y],
		window=owin(c(min(Data[,X]),max(Data[,X])),
					c(min(Data[,Y]),max(Data[,Y]))),
		marks=as.numeric(Data$DP))

BE=list()
i=min(Part)
while(i<=max(Part))
{
mat=matrix(NA,1000,2)
colnames(mat)=c("iter","Batty_Entropy")
	for(j in 1:1000)
		{
			mat[j,1]=j
			mat[j,2]=round(batty(Data2,partition=i)$rel.batty,3)
		}
BE[[i-1]]=mat
i=i+1
}

Data3=as.data.frame(rbind(BE[[1]],BE[[2]],BE[[3]],BE[[4]],BE[[5]],BE[[6]],
			BE[[7]],BE[[8]],BE[[9]]))
Data3$Partitions=rep(Part,each=1000)
Data3$Partitions=as.numeric(Data3$Partitions)

P4=ggplot(Data3,aes(x=factor(Partitions), y=Batty_Entropy,col=factor(Partitions),
			group=factor(Partitions))) + 
	geom_boxplot() + stat_summary(inherit.aes = FALSE,
	aes(x=factor(Partitions), y=Batty_Entropy,
		group=1),fun.y=median,geom="line",size=1)+
		theme_classic()+
	ggtitle(paste0("Batty's Entropy for"," ",Pathway))+
	labs(x ="Number of Partitions",color="Partitions", y ="Batty's Entropy")+
	theme(legend.position="bottom")

Entrp=round(summaryBy(Batty_Entropy ~ Partitions, data = Data3, 
          FUN = list(mean,sd,median)),3)

print(ggarrange(P3,P4,labels = c("A", "B")))
return(Entrp)
}

}