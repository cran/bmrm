\name{ordinalRegressionLoss}
\alias{ordinalRegressionLoss}
\title{The loss function for ordinal regression}
\usage{
ordinalRegressionLoss(w, x, y, C = "0/1", impl = c("loglin", "quadratic"),
  cache = NULL)
}
\arguments{
  \item{w}{weight vector where the function have to be
  evaluated}

  \item{x}{matrix of training instances (one instance by
  row)}

  \item{y}{integer vector of positive values (>=1)
  representing the training labels for each instance in x}

  \item{C}{the cost matrix to use, C[i,j] being the cost
  for predicting label i instead of label j.}

  \item{impl}{either the string "loglin" or "quadratic",
  that define the implementation to use for the computation
  of the loss.}

  \item{cache}{if NULL (which is the case at the first
  call) parameters values are checked}
}
\value{
a 2 element list (value,gradient) where "value" is the
value of the function at point w, and "gradient" is the
gradient of the loss function at w
}
\description{
The loss function for ordinal regression
}
\examples{
# -- Load the data
x <- data.matrix(iris[1:4])
y <- as.integer(iris$Species)

# -- Train the model
m <- bmrm(x,y,lossfun=ordinalRegressionLoss,LAMBDA=0.001,EPSILON_TOL=0.0001)
m2 <- bmrm(x,y,impl="quadratic",lossfun=ordinalRegressionLoss,LAMBDA=0.001,EPSILON_TOL=0.0001)

# -- plot predictions
f <- x \%*\% m$w
f2 <- x \%*\% m2$w
layout(1:2)
plot(y,f)
plot(f,f2,main="compare predictions of quadratic and loglin implementations")

# -- Compute accuracy
ij <- expand.grid(i=seq(nrow(x)),j=seq(nrow(x)))
n <- tapply(f[ij$i] - f[ij$j]>0,list(y[ij$i],y[ij$j]),sum)
N <- table(y[ij$i],y[ij$j])
print(n/N)
}
\references{
Teo et al.  Bundle Methods for Regularized Risk
Minimization JMLR 2010
}
\seealso{
bmrm
}

