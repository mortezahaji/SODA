

SODA_GWLCT <- function(Data=Intrp_Data,X="Longitude",Y="Latitude",Ovrp=0.25,
		BM="S100A10",kernel=c("Bisquare","Tricube"),method=c("fixed","adaptive"),
		bw=6,pthres=0.05,qthres=0.001, nbPermutations=500, silent=FALSE)
{

suppressPackageStartupMessages({
require(DT)
require(dplyr) 
require(corpcor)
require(qvcalc)
require(stringr)
require(MAVTgsa)
require(msigdbr)
require(ExperimentHub)
require(GSEABase)
})

cl=Data[,BM]

Lo=Data[,X]
La=Data[,Y]
DATA=Data[,!(names(Data) %in% c(X,Y))]


Gnames=colnames(Data)
#Retrieve human H (hallmark) gene set
msigdbr_df <- msigdbr(species = "human", category = "H")
pathwaysH = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

set=unique(unlist(pathwaysH[[Pathway]]))

GS=matrix(NA,length(Gnames),2)
GS[,1]=Gnames
GS[,2]=NA
for(i in 1:length(Gnames))
{
	GS[i,2]=ifelse(GS[i,1] %in% set,1,0)
}
colnames(GS)=c("Genes",Pathway)


  genes <- rownames(DATA)       # gene names of the microarray data 
  nb.Samples  <- ncol(DATA)     # nb of samples
  nb.GeneSets <- dim(GS)[2]     # nb of gene sets
  ngenes <- dim(GS)[1]
  Zp=matrix(0,nb.Samples,3)
  Zq=matrix(0,nb.Samples,3)
  ZZp=vector()
  ZZq=vector()
  
  Z1=array(0,c(nb.GeneSets,5,nb.Samples))
  Coordination=cbind(Lo[1:nb.Samples],La[1:nb.Samples])
  DistanceT <- dist(Coordination)
  dmatrix <- as.matrix(DistanceT)  
  
  if(kernel=="Bisquare" & method == 'adaptive' )
  { 
    Ne <- bw
    cat("\nBisquare Kernel: Adaptive\nNeightbours:", Ne)
  } 
  else if(kernel=="Bisquare" & method == 'fixed')
  {
    cat("\nBisquare Kernel: Fixed\nBandwidth:", bw)
  }
  else if(kernel=="Tricube" & method == 'adaptive'  )
  { 
    Ne <- bw
    cat("\nTricube Kernel: Adaptive\nNeightbours:", Ne)
  }
  else if(kernel=="Tricube" & method == 'fixed')
  {
    cat("\nTricube Kernel: Fixed\nBandwidth:", bw)
  }
  
  
  LCT_Global<-WLCT(GS,DATA ,cl,weight=rep(1/nb.Samples,nb.Samples), nbPermutations=1000, silent=FALSE)
  
  
  for(m in 1:nb.Samples){
    DNeighbour <- dmatrix[m,]
    DataSet <- data.frame(t(rbind(DATA,DNeighbour=DNeighbour)))
    DataSetSorted<- data.frame(t(DataSet[order(DataSet$DNeighbour),]))
    if(method == 'adaptive')
    { 
      SubSet1 <- DataSetSorted[,1:Ne]
      DATA1 <- DATA[,1:Ne]
      cl1 <- cl[1:Ne]
      Kernel_H <- max(SubSet1[ngenes+1,])
    } 
    else 
    { 
      if(method == 'fixed')
      {
        SubSet1 <- DataSetSorted[,1:bw]
        DATA1 <- DATA[,1:bw]
        cl1 <- cl[1:bw]
        Kernel_H <- bw
      }
    }
    
    if (kernel=="Bisquare"){
      Wts<-(1-(as.vector(unlist(SubSet1[ngenes+1,]))/Kernel_H)^2)^2
    }
    else if(kernel=="Tricube"){
      Wts<-(1-(as.vector(unlist(SubSet1[ngenes+1,]))/Kernel_H)^3)^3
    }
    
    LCTmodel=WLCT(GS,DATA1 ,cl1,weight=Wts, nbPermutations=nbPermutations, silent=FALSE)
    
    #Store in table
    Z1[,1,m]=rep(m,nb.GeneSets)
    Z1[,2,m]=LCTmodel$`GS stats`$`GS name`
    Z1[,3,m]=LCTmodel$`GS stats`$`GS size`
    Z1[,4,m]=LCTmodel$`GS stats`$`GS p-value`
    Z1[,5,m]=LCTmodel$`GS stats`$`GS q-value`
    ZZp[m]=str_c(Z1[,2,m][Z1[,4,m]<pthres], collapse = ",")
    ZZq[m]=str_c(Z1[,2,m][Z1[,5,m]<qthres], collapse = ",")
    Zp[m,]=c(m,sum(Z1[,4,m]<pthres),ZZp[m])
    Zq[m,]=c(m,sum(Z1[,5,m]<qthres),ZZq[m])
    
  }
  
  colnames(Z1)=c("Point","Pathway","No. of genes","P-value","Q-value")
  colnames(Zp)=c("Points","Number of Significant GS","Significant GS")
  colnames(Zq)=c("Points","Number of Significant GS","Significant GS")
  Zq1=data.frame(cbind(Zq[,2:3]))
  Zp1=data.frame(cbind(Lo,La,Zp[,2:3],Zq1))
  colnames(Zp1)=c("Longitude","Latitute","SGsF-p","SGs-p","SGsF-q","SGs-q")
  
  #figp <- plot_ly(as.data.frame(Zp1), x = ~Zp1[,1], y = ~Zp1[,2], z = ~Zp1[,3],
  #                text = ~Zp1[,4], marker = list(color = ~Zp1[,3], colorscale = "Rainbow",
  #                                              size=3, showscale = TRUE))
  #figp <- figp %>% add_markers()
  #figp <- figp %>% layout(title=paste("Geneset"),
  #                        scene = list(xaxis = list(title = 'Longitude'),
  #                                     yaxis = list(title = 'Latitude'),
  #                                     zaxis = list(title ="No.GeneSets")))
  #htmlwidgets::saveWidget(as_widget(figp),paste0("D:\\3DPlots\\pGS.html"))
  # 
  # 
  #   figq <- plot_ly(as.data.frame(Zp1), x = ~Zp1[,1], y = ~Zp1[,2], z = ~Zp1[,5],
  #                 text = ~Zp1[,6], marker = list(color = ~Zp1[,5], colorscale = "Rainbow",
  #                                              size=3, showscale = TRUE))
  # figq <- figq %>% add_markers()
  # figq <- figq %>% layout(title=paste("Geneset"),
  #                        scene = list(xaxis = list(title = 'Longitude'),
  #                                      yaxis = list(title = 'Latitude'),
  #                                      zaxis = list(title ="No.GeneSets")))
  # htmlwidgets::saveWidget(as_widget(figq),paste0("D:\\3DPlots\\qGS.html"))
  
  
  return(list("LCT"=LCT_Global, "GWLCT"=Zp1))  
}


