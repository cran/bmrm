

# L2 regularizer
# TODO: Enable hot starting if possible (do not forgot to check/rescale the starting point to ensure it is feasible)
l2 <- function(A,b,LAMBDA) {
  D <- tcrossprod(A,A) / LAMBDA
  Ale <- matrix(1,1,nrow(D))
  lb <- matrix(0,nrow(D),1)
  ub <- matrix(1,nrow(D),1)
  opt <- ipop(-b,D,Ale,0,lb,ub,1,margin=0.01)
  alpha <- primal(opt)
  value <- t(alpha) %*% D %*% alpha / 2 - b %*% alpha
  w <- -crossprod(A,alpha) / LAMBDA
  return(list(
    w = as.vector(w),
    objectiveValue = -value,
    regularizationValue = LAMBDA*(t(w) %*% w)/2
  ))
}


# L1 regularizer
l1 <- function(A,b,LAMBDA) {
  lp <- initProbCLP()
  on.exit(delProbCLP(lp))
  setLogLevelCLP(lp,0)
  
  setObjDirCLP(lp,-1)
  loadProblemCLP(lp,ncols=1,nrows=nrow(A),
                 ia=seq_len(nrow(A))-1L,ja=c(0L,nrow(A)),ra=rep(-1,nrow(A)),
                 lb=0,obj_coef=-1,rub=-b)
  lb <- rep(0,ncol(A))
  ub <- rep(.Machine$double.xmax,ncol(A))
  obj <- rep(-LAMBDA,ncol(A))
  addColsCLP(lp,ncol(A),lb=lb,ub=ub,obj=obj,colst=seq(0L,length(A),by=nrow(A)),rows=row(A)-1L,val=A)
  addColsCLP(lp,ncol(A),lb=lb,ub=ub,obj=obj,colst=seq(0L,length(A),by=nrow(A)),rows=row(A)-1L,val=-A)
  solveInitialCLP(lp)

  if (getSolStatusCLP(lp)!=0) warning("issue in the LP solver:",status_codeCLP(getSolStatusCLP(lp)))

  W <- getColPrimCLP(lp)
  u <- W[-1][1:ncol(A)]
  v <- W[-1][-1:-ncol(A)]
  w <- u - v
  return(list(
    w = w,
    objectiveValue = -getObjValCLP(lp),
    regularizationValue = LAMBDA*sum(abs(w))
  ))
}





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
#' @param w0 a numeric vector used to initialize the minimization process
#' @param verbose a length one logical. Show progression of the convergence on stdout
#' @param ... additional argument passed to the loss function
#' @return a list of 2 fileds: "w" the optimized weight vector; "log" a data.frame showing the trace of important values in the optimization process.
#' @export
#' @import clpAPI
#' @import kernlab
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
#'   train.prediction.model <- function(x,y,lossfun=hingeLoss,...) {
#'     m <- bmrm(x,y,lossfun=lossfun,...)
#'     m$f <- x %*% m$w
#'     m$y <- sign(m$f)
#'     m$contingencyTable <- table(y,m$y)
#'     print(m$contingencyTable)
#'     return(m)
#'   }
#'   
#'   # -- train scalar prediction models with maxMarginLoss and fbetaLoss 
#'   models <- list(
#'     svm_L1 = train.prediction.model(x,y,lossfun=hingeLoss,LAMBDA=0.01,regfun='l1'),
#'     svm_L2 = train.prediction.model(x,y,lossfun=hingeLoss,LAMBDA=0.1,regfun='l2'),
#'     f1_L1 = train.prediction.model(x,y,lossfun=fbetaLoss,LAMBDA=0.01,regfun='l1')
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
bmrm <- function(...,LAMBDA=1,MAX_ITER=100,EPSILON_TOL=0.01,lossfun=hingeLoss,regfun=c('l2','l1'),w0=0,verbose=FALSE) {
# Inner Variables:
#       A: current gradient matrix
#       b: current optimization vector
#    loss: current loss value
#     reg: current optimization status
		
	regfun <- match.arg(regfun)
	regfun <- get(regfun, mode = "function",parent.env(environment()))
	loss <- lossfun(w0,cache=NULL,...)
	A <- matrix(as.vector(loss$gradient),1,)
	b <- loss$value - crossprod(rep_len(w0,length(loss$gradient)),as.vector(loss$gradient))

	loop <- list(loss=loss$value,regVal=NA,lb=NA,ub=NA,epsilon=NA,nnz=NA)
	for (i in 2:MAX_ITER) {
		reg <- regfun(A,b,LAMBDA) # Find minimizer w_i = argmin(J_i(w)), using [A,b]
		
		# -- compute gradient ai and offset bi, and update [A,b]
		loss <- lossfun(reg$w,cache=loss$cache,...)
    loss$gradient <- as.vector(loss$gradient)
		A <- rbind(A,loss$gradient)
		b <- c(b,loss$value - crossprod(reg$w,loss$gradient))

		# -- check convergeance (duality gap epsilon close to zero)
		loop$loss[i] <- loss$value
		loop$regVal[i] <- reg$regularizationValue
		loop$lb[i] <- reg$objectiveValue
    loop$nnz[i] <- sum(reg$w != 0)
    loop$numCst[i] <- length(b)
		if (i>3) {
      Ri <- max(0,(A %*% reg$w) + b)
			loop$ub[i] <- min(loop$ub[i-1],Ri + loop$regVal[i],na.rm=TRUE)
			loop$epsilon[i] <- loop$ub[i] - loop$lb[i]
			if (loop$epsilon[i] < EPSILON_TOL) break
		}
    if (verbose) {cat(paste("i=",i,sep=''),paste(names(loop),sapply(loop,'[',i),sep='='),sep=',');cat("\n")}
	}
	if (i >= MAX_ITER) warning('max # of itertion exceeded')
	return(list(
    w=reg$w,
    log=as.data.frame(loop)
  ))
}

