\name{quantileRegressionLoss}
\alias{quantileRegressionLoss}
\title{The loss function to perform a quantile regression}
\usage{
quantileRegressionLoss(w, x, y, q = 0.5, cache = NULL)
}
\arguments{
  \item{w}{weight vector where the function have to be
  evaluated}

  \item{x}{matrix of training instances (one instance by
  row)}

  \item{y}{numeric vector of values representing the
  training labels for each instance in x}

  \item{q}{a numeric value in the range [0-1] defining
  quantile value to consider}

  \item{cache}{if NULL (which is the case at the first
  call) parameters values are checked}
}
\value{
a 2 element list (value,gradient) where "value" is the
value of the function at point w, and "gradient" is the
gradient of the loss function at w
}
\description{
The loss function to perform a quantile regression
}
\references{
Teo et al.  Bundle Methods for Regularized Risk
Minimization JMLR 2010
}
\seealso{
bmrm
}

