
SODA_Wass=function(Data1,Data2,X,Y,BM,cutoff,Rad,Method=c("TS","OS"))
{

suppressPackageStartupMessages({
    require("waddR")
    require("SingleCellExperiment")
    require("BiocFileCache")
})

Coord=data.frame(Data2[,X],Data2[,Y])
names(Coord)=c(X,Y)
Data2=Data2[,!(names(Data2) %in% c(X,Y))]

ColN=dim(Data2)[2]
NamCol=names(Data2)

P_matrix=matrix(NA,1,ColN+1)
colnames(P_matrix)=c("CAFhi","CAFlo",
NamCol[-which(NamCol %in% BM)])
rownames(P_matrix)=BM

SCell=which(Data1[,BM]>=cutoff)
LCel=c()
Cell=c()
	for(i in 1:length(SCell))
	{
		
		for(j in 1:dim(Coord)[1])
			{
				dx = abs(Data1[SCell[i],1]-Coord[j,1])
				dy = abs(Data1[SCell[i],2]-Coord[j,2])
				R =Rad
					if (dx^2 + dy^2 <= R^2)
					{
						LCel=c(LCel,j)
						Cell=c(Cell,rownames(Coord)[j])
					}
			}
	}

		Data2$LNL=NA
		Data2$LNL[rownames(Data2) %in% unique(Cell)]<-1 #High
		Data2$LNL[!(rownames(Data2) %in% unique(Cell))]<-0 #Low
		T=table(Data2$LNL)

		D2=Data2[,c(dim(Data2)[2])]
		D3=Data2[,-which(names(Data2) %in% c(BM,"LNL"))]
		
		#Wasserstein Distance Test
if(Method=="TS")
{
		res <- wasserstein.sc(t(D3),D2, 
			method="TS",permnum=1000,seed=24)
		P_Value=res[,9]
		P_Adjusted=res[,18]
P=data.frame(N_High=rep(T[2],length(P_Value)),
		N_Low=rep(T[1],length(P_Value)),
		Pvalue=P_Value,
		Adjusted_Pvalue=P_Adjusted)
P=P[order(P$Adjusted_Pvalue),]

}else if(Method=="OS"){
		res <- wasserstein.sc(t(D3),D2, 
			method="OS",permnum=1000,seed=24)
		P_Combined=res[,9]
		P_Adjusted_Combined=res[,16]
P=data.frame(N_High=rep(T[2],length(P_Combined)),
		N_Low=rep(T[1],length(P_Combined)),
		Pvalue=P_Combined,
		Adjusted_Pvalue=P_Adjusted_Combined)
P=P[order(P$Adjusted_Pvalue),]
}

return(P)
}
