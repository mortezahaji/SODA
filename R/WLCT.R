
GS.format.dataframe.to.list <- function(GS){
  if(is.data.frame(GS)){
    genes <- rownames(GS)
    L <- NULL
    for(ags in names(GS)){
      w <- which(GS[,ags]==1)
      if(length(w)>0)  {
        L <- c(L,list(genes[w]))
        names(L)[length(L)] <- ags
      }
    }
    L
  }else{
    GS
  }
}



T2.like.SAMGS <- function(DATA, cl){
  cl<-as.matrix(c(as.numeric(cl)))
  DATA<-as.matrix(DATA)
  cl.DATA<-DATA%*%cl
  sum(cl.DATA^2)
}

   WLCT <- function(GS, DATA, cl,weight, nbPermutations){

  genes <- rownames(DATA)       # gene names of the microarray data 
  nb.Samples  <- ncol(DATA)     # nb of samples
  nb.GeneSets <- dim(GS)[2]     # nb of gene sets

    GS <-  GS.format.dataframe.to.list(GS);
    # change format of GS from data.frame to                         
    # list
    GS <-  lapply(GS,function(z) as.numeric(which(genes %in% z)));
    GS.sizes <- sapply(GS,length) # size of each gene set
    # numericalized index of each GS
    GS.data <- lapply(GS, function(z) as.matrix(DATA[z, ],ncol=nb.Samples));
    # creat data of each GS (rows=genes,  
    # columns=samples)
    GS.data <- lapply(GS.data,function(z) scale(t(z))); 
    # standardized genes in each GS 
    # (columns=gene, rows=samples)
    
    cl=scale(cl) #standardized response
    
    # (2) Eigen-decomposition of shrinkage pooled covariance matrix for each GS
    
    Cov.Pooled<-lapply(GS.data, function(z) cov.shrink(z,verbose=FALSE, lambda.var=0,w=ww));
    # pooled covariance of genes in each GS
    for (i in 1:nb.GeneSets){
      EIGEN.decom<-eigen(Cov.Pooled[[i]]);
      # eigen decomposition of pooled covariance for each GS
      D<-EIGEN.decom$values;      # shrinkage by adding a positive constant s0
      U<-EIGEN.decom$vectors;
      
      GS.data[[i]]<-t(GS.data[[i]]%*%U)/sqrt(D)
      # adjust data of each GS (rows=genes, columns=samples)
      
    }

    # (3) T-like stats obtained on 'true' data
    sam.sumsquareT.obs  <- sapply(GS.data, function(z) T2.like.SAMGS(z,cl))
    # the T-like statistics obtained on 'true' data
    
    # (4) stats obtained on 'permuted' data
    sam.sumsquareT.permut <- matrix(NA,nbPermutations,nb.GeneSets)
    for(i in 1:nbPermutations) {
      ind <- sample(nb.Samples)
      sam.sumsquareT.permut[i,] <- sapply(GS.data, function(z) T2.like.SAMGS(z[,ind],cl))
    }
    
    # (5) p-value and q-value
    GeneSets.pval <- apply(t(sam.sumsquareT.permut) >= sam.sumsquareT.obs,1,sum)/nbPermutations
    
    if(nb.GeneSets>=2){
      GeneSets.qval=qvalue(GeneSets.pval,lambda=0)$qvalues
      
      res <- as.data.frame(cbind("GS size"              = GS.sizes,
                                 "GS p-value" 	       = GeneSets.pval,
                                 "GS q-value"           = GeneSets.qval ))
      res <- cbind(res,"GS name"= names(GS))[c(4,1:3)]
    }
    
    if(nb.GeneSets==1){
      #if there is only one set, no need to calculate q-value.
      res <- as.data.frame(cbind("GS size"              = GS.sizes, ##GeneSets.sizes,
                                 "GS p-value" 	       = GeneSets.pval))
      res <- cbind(res,"GS name"= names(GS))[c(3,1:2)]
    }
    rownames(res)<-NULL
    return(list("GS stats"=res))
  
}
