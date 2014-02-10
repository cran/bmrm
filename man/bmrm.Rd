\name{bmrm}
\alias{bmrm}
\title{Bundle Methods for Regularized Risk Minimization}
\usage{
bmrm(..., LAMBDA = 1, MAX_ITER = 100, EPSILON_TOL = 0.01,
  lossfun = hingeLoss, regfun = c("l1", "l2"), verbose = FALSE)
}
\arguments{
  \item{lossfun}{the loss function to use in the
  optimization (e.g.: hingeLoss, softMarginVectorLoss).
  The function must evaluate the loss value and its
  gradient for a given point vector (w).  The function must
  be of the form lossfun(w,...,cache=NULL), i.e. accept as
  first parameter the vector of weight w, and unused
  arguments to bmrm().  The return value must be a
  list(value,gardient,cache), where value is the numeric
  value of the loss for w, and gradient is the gradient
  vector of the function at point w.  The "cache" parameter
  and the "cache" element in the return value can be used
  to store variable from one call to the next call.  The
  "cache" parameter is set to NULL at the first call, and
  is set to the previous returned "cache" value at next
  calls.}

  \item{LAMBDA}{control the regularization strength in the
  optimization process. This is the value used as
  coefficient of the regularization term.}

  \item{MAX_ITER}{the maximum number of iteration to
  perform. The function stop with a warning message if the
  number of iteration exceed this value}

  \item{EPSILON_TOL}{control optimization stoping criteria:
  the optimization end when the optimization gap is below
  this threshold}

  \item{regfun}{type of regularization to consider in the
  optimization. It can either be the character string "l1"
  for L1-norm regularization, or "l2" (default) for L2-norm
  regularization.}

  \item{verbose}{a length one logical. Show progression of
  the convergence on stdout}

  \item{...}{additional argument passed to the loss
  function}
}
\value{
a list of 2 fileds: "w" the optimized weight vector; "log"
a data.frame showing the trace of important values in the
optimization process.
}
\description{
Implement Bundle Methods for Regularized Risk Minimization
as described in Teo et. al 2007.
}
\examples{
# -- Create a 2D dataset with the first 2 features of iris, with binary labels
  x <- data.matrix(iris[1:2])
  y <- c(-1,1,1)[iris$Species]

  # -- Add a constant dimension to the dataset to learn the intercept
  x <- cbind(x,1)

  # -- train scalar prediction models with maxMarginLoss and fbetaLoss
  models <- list(
    svm_L1 = bmrm(x,y,lossfun=hingeLoss,LAMBDA=0.1,regfun='l1',verbose=TRUE),
    svm_L2 = bmrm(x,y,lossfun=hingeLoss,LAMBDA=0.1,regfun='l2',verbose=TRUE),
    f1_L1 = bmrm(x,y,lossfun=fbetaLoss,LAMBDA=0.01,regfun='l1',verbose=TRUE)
  )

  # -- Plot the dataset and the predictions
  layout(matrix(1:2,1,2))
  plot(x,pch=20+y,main="dataset & hyperplanes")
  legend('bottomright',legend=names(models),col=seq_along(models),lty=1,cex=0.75,lwd=3)
  for(i in seq_along(models)) {
    m <- models[[i]]
    if (m$w[2]!=0) abline(-m$w[3]/m$w[2],-m$w[1]/m$w[2],col=i,lwd=3)
  }

  rx <- range(na.rm=TRUE,1,unlist(lapply(models,function(e) nrow(e$log))))
  ry <- range(na.rm=TRUE,0,unlist(lapply(models,function(e) e$log$epsilon)))
  plot(rx,ry,type="n",ylab="epsilon gap",xlab="iteration",main="evolution of the epsilon gap")
  for(i in seq_along(models)) {
    m <- models[[i]]
    lines(m$log$epsilon,type="o",col=i,lwd=3)
  }
}
\author{
Julien Prados
}
\references{
Teo et al.  A Scalable Modular Convex Solver for
Regularized Risk Minimization.  KDD 2007
}
\seealso{
\code{\link{hingeLoss}} \code{\link{softMarginVectorLoss}}
}

