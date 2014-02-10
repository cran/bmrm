


setClass("rrmSolver",contains="VIRTUAL",slots=c(LAMBDA="numeric"))
setMethod("initialize",signature(.Object="rrmSolver"),function(.Object,LAMBDA) {
  .Object@LAMBDA <- LAMBDA
  .Object
})
setGeneric("destroy",function(.Object) standardGeneric("destroy"))



#
# -- Define Solver for L1 regularization
# 
setClass("rrmSolverL1",contains="rrmSolver",slots=c(lp="ANY"))
setMethod("initialize",signature(.Object="rrmSolverL1"),function(.Object,...) {
  .Object <- callNextMethod()
  .Object@lp <- initProbCLP()
  setLogLevelCLP(.Object@lp,0)
  setObjDirCLP(.Object@lp,-1)
  .Object
})
setMethod("destroy","rrmSolverL1",function(.Object) {
  delProbCLP(.Object@lp)
})
setMethod("rbind2",signature(x="rrmSolverL1",y="ANY"),function(x,y,b,...) {
  nc <- 2L*length(y)+1L
  if (getNumRowsCLP(x@lp)<1L) {
    resizeCLP(x@lp,0L,nc)
    chgObjCoefsCLP(x@lp,c(-1,rep_len(-x@LAMBDA,nc-1L)))
    chgColLowerCLP(x@lp,rep_len(0,nc))
    chgColUpperCLP(x@lp,rep_len(.Machine$double.xmax,nc))
  }
  addRowsCLP(x@lp,1L,-.Machine$double.xmax,-b,c(0L,nc),seq_len(nc)-1L,c(-1,y,-y))
})
setMethod("solve",signature(a="rrmSolverL1"),function(a,...) {
  nc <- getNumColsCLP(a@lp)
  solveInitialCLP(a@lp)
  if (getSolStatusCLP(a@lp)!=0) warning("issue in the LP solver:",status_codeCLP(getSolStatusCLP(a@lp)))
  W <- getColPrimCLP(a@lp)
  u <- W[-1][seq_len(nc%/%2)]
  v <- W[-1][-seq_len(nc%/%2)]
  w <- u - v
  return(list(w = w, obj = -getObjValCLP(a@lp), regval = sum(abs(w))))  
})






#
# -- Define Solver for L2 regularization
# 
setClass("rrmSolverL2",contains="rrmSolver",slots=c(env="environment"))
setMethod("initialize",signature(.Object="rrmSolverL2"),function(.Object,...) {
  .Object <- callNextMethod()
  .Object@env <- new.env()
  with(.Object@env,{A <- matrix(,0,0);b <- numeric(0)})
  .Object
})
setMethod("destroy","rrmSolverL2",function(.Object) {
})
setMethod("rbind2",signature(x="rrmSolverL2",y="ANY"),function(x,y,b,...) {
  if (ncol(x@env$A)!=length(y)) dim(x@env$A)[2] <- length(y)
  x@env$A <- rbind(x@env$A,y)
  x@env$b <- c(x@env$b,b)
})
setMethod("solve",signature(a="rrmSolverL2"),function(a,...) {
  Ale <- matrix(1,1,nrow(a@env$A))
  H <- tcrossprod(a@env$A)
  opt <- ipop(c=-a@LAMBDA*a@env$b,H,Ale,0,rep_len(0,ncol(Ale)),rep_len(1,ncol(Ale)),1,sigf=5)
  alpha <- primal(opt)
  w <- as.vector(-crossprod(a@env$A,alpha) / a@LAMBDA)
  regval <- 0.5*crossprod(w)
  R <- max(0,a@env$A %*% w + a@env$b)
  return(list(w = as.vector(w), obj = a@LAMBDA*regval+R, regval = regval))
})


#' Bundle Methods for Regularized Risk Minimization
#' 
#' Implement Bundle Methods for Regularized Risk Minimization as described in Teo et. al 2007.
#' 
#' @param lossfun the loss function to use in the optimization (e.g.: hingeLoss, softMarginVectorLoss). 
#'   The function must evaluate the loss value and its gradient for a given point vector (w).
#'   The function must be of the form lossfun(w,...,cache=NULL), i.e. accept as first parameter the vector of weight w, and unused arguments to bmrm().
#'   The return value must be a list(value,gardient,cache), where value is the numeric value of the loss for w, and gradient is the gradient vector of the function at point w.
#'   The "cache" parameter and the "cache" element in the return value can be used to store variable from one call to the next call. 
#'   The "cache" parameter is set to NULL at the first call, and is set to the previous returned "cache" value at next calls.
#' @param LAMBDA control the regularization strength in the optimization process. This is the value used as coefficient of the regularization term.
#' @param MAX_ITER the maximum number of iteration to perform. The function stop with a warning message if the number of iteration exceed this value
#' @param EPSILON_TOL control optimization stoping criteria: the optimization end when the optimization gap is below this threshold
#' @param regfun type of regularization to consider in the optimization. It can either be the character string "l1" for L1-norm regularization, 
#'   or "l2" (default) for L2-norm regularization.
#' @param verbose a length one logical. Show progression of the convergence on stdout
#' @param ... additional argument passed to the loss function
#' @return a list of 2 fileds: "w" the optimized weight vector; "log" a data.frame showing the trace of important values in the optimization process.
#' @export
#' @import clpAPI
#' @import kernlab
#' @import methods
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @author Julien Prados
#' @seealso \code{\link{hingeLoss}} \code{\link{softMarginVectorLoss}}
#' @examples
#'   # -- Create a 2D dataset with the first 2 features of iris, with binary labels
#'   x <- data.matrix(iris[1:2])
#'   y <- c(-1,1,1)[iris$Species]
#'   
#'   # -- Add a constant dimension to the dataset to learn the intercept
#'   x <- cbind(x,1)
#'   
#'   # -- train scalar prediction models with maxMarginLoss and fbetaLoss 
#'   models <- list(
#'     svm_L1 = bmrm(x,y,lossfun=hingeLoss,LAMBDA=0.1,regfun='l1',verbose=TRUE),
#'     svm_L2 = bmrm(x,y,lossfun=hingeLoss,LAMBDA=0.1,regfun='l2',verbose=TRUE),
#'     f1_L1 = bmrm(x,y,lossfun=fbetaLoss,LAMBDA=0.01,regfun='l1',verbose=TRUE)
#'   )
#'   
#'   # -- Plot the dataset and the predictions
#'   layout(matrix(1:2,1,2))
#'   plot(x,pch=20+y,main="dataset & hyperplanes")
#'   legend('bottomright',legend=names(models),col=seq_along(models),lty=1,cex=0.75,lwd=3)
#'   for(i in seq_along(models)) {
#'     m <- models[[i]]
#'     if (m$w[2]!=0) abline(-m$w[3]/m$w[2],-m$w[1]/m$w[2],col=i,lwd=3)
#'   }
#'   
#'   rx <- range(na.rm=TRUE,1,unlist(lapply(models,function(e) nrow(e$log))))
#'   ry <- range(na.rm=TRUE,0,unlist(lapply(models,function(e) e$log$epsilon)))
#'   plot(rx,ry,type="n",ylab="epsilon gap",xlab="iteration",main="evolution of the epsilon gap")
#'   for(i in seq_along(models)) {
#'     m <- models[[i]]
#'     lines(m$log$epsilon,type="o",col=i,lwd=3)
#'   }
#'   
bmrm <- function(...,LAMBDA=1,MAX_ITER=100,EPSILON_TOL=0.01,lossfun=hingeLoss,regfun=c('l1','l2'),verbose=FALSE) {		
	regfun <- match.arg(regfun)
  rrm <- switch(regfun,l1=new("rrmSolverL1",LAMBDA=LAMBDA),l2=new("rrmSolverL2",LAMBDA=LAMBDA))
  on.exit(destroy(rrm))
  loss <- list(cache=NULL)
  opt <- list(w=0,regval=0)
  ub <- +Inf
	log <- list(loss=numeric(),regVal=numeric(),lb=numeric(),ub=numeric(),epsilon=numeric(),nnz=integer())
	for (i in 1:MAX_ITER) {
	  loss <- lossfun(opt$w,cache=loss$cache,...)
    loss$gradient <- as.vector(loss$gradient)
    opt$w <- rep_len(opt$w,length(loss$gradient))
    rbind2(rrm,loss$gradient,loss$value - crossprod(opt$w,loss$gradient))
	  ub <- min(ub,loss$value + LAMBDA*opt$regval)
	  opt <- solve(rrm)
    lb <- opt$obj
    log$loss[i]<-loss$value;log$regVal[i]<-opt$regval;log$lb[i]<-lb;log$ub[i]<-ub;log$epsilon[i]<-ub-lb;log$nnz[i]<-sum(opt$w!=0)
	  if (verbose) {cat(sprintf("i=%d,eps=%g (=%g-%g),nnz=%d,loss=%g,reg=%g\n",i,log$epsilon[i],log$ub[i],log$lb[i],log$nnz[i],log$loss[i],log$regVal[i]))}
    if (ub-lb < EPSILON_TOL) break
	}
	if (i >= MAX_ITER) warning('max # of itertion exceeded')
	return(list(w=opt$w,log=as.data.frame(log)))
}





