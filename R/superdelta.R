## turn a vector of t-statistics into two sided pvalues
.t2p <- function(x, df, side="abs"){
  if (side=="abs"){
    return(2-2*pt(abs(x),df))
  } else if (side=="lower"){
    return(pt(x, df))
  } else if (side=="upper"){
    return(1-pt(x, df))
  } else {
    stop(paste("side must be one of the following three options: abs (default), lower, or upper."))
  }
}

## Median; fold; trim; median
mftm <- function(tvec, trim=0.2){
  med1 <- median(tvec, na.rm=TRUE)
  tvec2 <- abs(tvec-med1)
  t.est <- median(tvec[tvec2<=quantile(tvec2, 1-trim, na.rm=TRUE)])
  return(t.est)
}

## cluster tvec by mclust
.tclust1 <- function(tvec, ...) {
  mm.i <- summary(mclustBIC(tvec, G=3, modelName="V", ...), tvec)
  NC1 <- sum(mm.i$classification==1); NC2 <- sum(mm.i$classification==2); NC3 <- sum(mm.i$classification==3)
  largest.cluster <- which.max(c(NC1, NC2, NC3))
  t.est <- as.real(mm.i$parameters$mean[largest.cluster])
  return(t.est)
}

.tclust2 <- function(tvec, trim=0.2, ...){
  ## mclust is not that robust, so use the following wrapper as a backup.
  .substitute <- function(x) {
    warning("Mclust does not converge, use the robust method instead.")
    mftm(tvec, trim=trim)
  }
  t.est <- tryCatch(.tclust1(tvec, ...), warning=.substitute, error=.substitute)
  return(t.est)
}

## this function takes tmat, which is a matrix of
## t-statistics. tmat[i,j] is the tstat associated with the ith gene
## normalized by the jth gene.  It estimates one best t-stat according
## to different methods.
.est.t <- function(tmat, method=c("robust", "mean", "median", "mclust"), trim=0.2, ...){
  method=match.arg(method)
  if (method=="robust"){
    t.est <- apply(tmat, 1, mftm, trim=trim)
  } else if (method=="median"){
    t.est <- rowMedians(tmat, na.rm=TRUE)
  } else if (method=="mean"){
    t.est <- rowMeans(tmat, na.rm=TRUE)
  } else if (method=="mclust"){
    t.est <- apply(tmat, 1, .tclust2)
  } else {
    stop(paste("Method", method, "is not implemented!"))
  }
  ## 7/30/2011.  Try t-stat with discounted variance.
  return(sqrt(2)*t.est)
}

superdelta <- function(X, classlabel, test="t", side="abs", na=.mt.naNUM,
                       methods="robust", trim=0.2, baseline="auto", ...)
{
  ngenes <- dim(X)[1]; nslides <- dim(X)[2]
  teststat <- matrix(0, nrow=ngenes, ncol=length(methods))
  colnames(teststat) <- methods; rownames(teststat) <- rownames(X)

  if(is.factor(classlabel)) classlabel<-unclass(classlabel)-1
  extra<-max(classlabel)+1
  mt.checkothers(na=na,nonpara="n")
  tmp<-mt.transformX(X,classlabel,test,na,nonpara="n")
  options<-c(test,"abs","y"); #"abs"  and "y" has no meaning here
  ## decide whether a heuristic should be used or not
  if (baseline=="auto") {
    baseline <- ifelse(ngenes>1000, 1000, 0)
  } else if (baseline=="all"){
    baseline <- 0
  }

  if (length(baseline)>1){                #assume this is a vector of pre-selected baseline genes (heuristic)
    baseline.genes <- baseline
  } else if (baseline>0){                      #use heuristic
    baseline.genes <- sample(ngenes, baseline)
  } else {                              #use all genes
    baseline.genes <- 1:ngenes
  }
  baseline <- length(baseline.genes)

  res<-.C("superdelta_stats",as.double(tmp$X),as.integer(tmp$m),
          as.integer(tmp$n),as.integer(tmp$classlabel),as.double(na),
          as.integer(baseline.genes-1), as.integer(baseline),
          teststat=double(tmp$m*(baseline-1)),as.character(options),
          as.integer(extra), PACKAGE="xmulttest")$teststat
  res[abs(res)>=0.9*1e20]<-NA
  tmat <- matrix(res, nrow=ngenes)

  ## Now estimate the centers
  for (mm in methods){
    teststat[,mm] <- .est.t(tmat, method=mm, trim=trim, ...)
  }

  ## two sided or one sided test; different test statistics, etc
  if (test=="t"){
    rawp=.t2p(teststat, df=nslides-2, side=side)
  } else {
    stop(paste("Test statistic", test, "is not implemented yet!"))
  }

  ## two representations
  if (ncol(teststat)==1){
    rr <- data.frame(teststat=teststat, rawp=rawp)
    colnames(rr) <- c("teststat", "rawp")
    return(rr)
  } else {
    return(list(teststat=teststat, rawp=rawp))
  }
}

