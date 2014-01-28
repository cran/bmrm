\name{softMarginVectorLoss}
\alias{softMarginVectorLoss}
\title{Soft Margin Vector Loss function for multiclass SVM}
\usage{
softMarginVectorLoss(w, x, y, l = "0/1", cache = NULL)
}
\arguments{
  \item{w}{weight vector where the function have to be
  evaluated}

  \item{x}{instance matrix, where x(t,) defines the
  features of instance t}

  \item{y}{target vector where y(t) is an integer encoding
  target of x(t,)}

  \item{l}{loss matrix. l(t,p(t)) must be the loss for
  predicting target p(t) instead of y(t) for instance t. By
  default, the parameter is set to character value "0/1" so
  that the loss is set to a 0/1 loss matrix.}

  \item{cache}{if NULL, the parameters are checked. If set
  to a list with element "l", the this element will be used
  to replace the parameter l without any checking}
}
\value{
a 2 element list (value,gradient) where "value" is the
value of the function at point w, and "gradient" is the
gradient of the loss function at w
}
\description{
Soft Margin Vector Loss function for multiclass SVM
}
\examples{
# -- Load the data
  x <- data.matrix(iris[1:2])
  y <- as.integer(iris$Species)

  # -- Add a constant dimension to the dataset to learn the intercept
  cst <- sqrt(max(rowSums(x*x)))
  x <- cbind(x,cst)

  # -- compute the loss for predicting a given category instead of the true category
  l <- matrix(1,length(y),max(y))
  l[cbind(seq_along(y),y)] <- 0

  # -- train a multi-class SVM & compute the predictions
  train.multiclass <- function(x,y,...) {
    m <- bmrm(x,y,lossfun=softMarginVectorLoss,...)
    m$w <- matrix(m$w,ncol(x))
    m$f <- x \%*\% m$w
    m$y <- max.col(m$f)
    m$contingencyTable <- table(y,m$y)
    print(m$contingencyTable)
    return(m)
  }
  # train a L1-regularized multi-class SVM
  m <- train.multiclass(x,y,l=l,MAX_ITER=50,regfun='l1',LAMBDA=1)

  # -- Plot the dataset and the decision boundaries
  gx <- seq(min(x[,1]),max(x[,1]),length=200) # positions of the probes on x-axis
  gy <- seq(min(x[,2]),max(x[,2]),length=200) # positions of the probes on y-axis
  Y <- outer(gx,gy,function(a,b){
     max.col(cbind(a,b,cst) \%*\% m$w)
  }) # matrix of predictions for all probes
  layout(matrix(1:2,1))
  image(gx,gy,Y,asp=1,main="dataset & decision boundaries")
  points(x,pch=19+y)
  plot(m$log$epsilon,type="o",ylab="epsilon gap",xlab="iteration")
}
\references{
Teo et al.  A Scalable Modular Convex Solver for
Regularized Risk Minimization.  KDD 2007
}

