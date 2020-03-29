#' @title AStest
#'
#' @description
#' \code{Astest} is a function to implement GWAS using covariates.
#'
#' @details
#' This function use GLM modeling snp, covariates and phenotype.
#'
#' @param y phenotype data: n by 1
#' @param X genotype data: n by m
#' @param C covariated defined by user: n by t
#'
#' @return p values list for each SNP
#'
#' @examples
#' P_values <- AStest(y=myY, X=myGD, C=myC)
#'
#' @references
#' Lecture 13_GLM in Crops545 Statistical Genomics
#'
#' @export
AStest=function(y,X,C){
  y=y[,2]
  C0=C[,-1] #remove taxa name of covariates
  G=X[,-1]

  #these code come from lecture of GLM
  n=nrow(G)
  m=ncol(G)
  P=matrix(NA,1,m)
  for (i in 1:m){
    x=G[,i]
    if(max(x)==min(x)){
      p=1}else{
        X1=cbind(mean(y), C0, x)
        X1=as.matrix(X1)
        LHS=t(X1)%*%as.matrix(X1)
        C1=solve(LHS)
        RHS=t(X1)%*%y
        b=C1%*%RHS
        yb=as.matrix(X1)%*%b
        e=y-yb
        n=length(y)
        ve=sum(e^2)/(n-1)
        vt=C1*ve
        t=b/sqrt(diag(vt))
        p=2*(1-pt(abs(t),n-2))
      } #end of testing variation
    P[i]=p[length(p)]
  } #end of looping for markers
  return(P)
}


#' @title GWAS
#'
#' @description
#' \code{GWAS} is a function to implement GWAS using covariates and PCs from PCA.
#'
#' @details
#' This function use GLM modeling snp, covariates, PCs and phenotype.
#'
#' @param y phenotype data: n by 1
#' @param X genotype data: n by m
#' @param C covariated defined by user: n by t
#' @param pcas the number of pcas used in GWAS. Defalut is 3. We can use 'False' to cancel the usage of pcas
#'
#' @return p values list for each SNP
#'
#' @examples
#' # to include 3 PCs in calculation
#' P_values01 <- GWAS(y=myY, X=myGD, C=myC, pcas=3)
#' # Dosen't include PCs in calculation
#' P_values02 <- GWAS(y=myY, X=myGD, C=myC, pcas=False)
#'
#' @references
#' Lecture 13_GLM in Crops545 Statistical Genomics
#'
#' @export
GWAS=function(y,X,C, pcas=3){
  C0=C[,-1] #remove taxa name of covariates
  #GWAS related content
  if(pcas == FALSE){
    C_all = C0
  }else{
    G=X[,-1] #remove taxa name of genotype
    ###remove PC that are in linear dependent to the covariates
    PCA=prcomp(G)
    r<-cor(PCA$x,C0)
    index1=r > 0.9
    r[index1]=NA
    r_remain=na.omit(r)
    keep_PC<-PCA$x[,colnames(PCA$x) %in% rownames(r_remain)]
    C_all = cbind(C0, keep_PC[,1:pcas])
  }
  P = AStest(y, X, C=C_all)
  return(P)
}
