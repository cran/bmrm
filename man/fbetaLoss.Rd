\name{fbetaLoss}
\alias{fbetaLoss}
\title{F beta score loss function}
\usage{
fbetaLoss(w, x, y, beta = 1, cache = NULL)
}
\arguments{
  \item{w}{weight vector where the function have to be
  evaluated}

  \item{x}{matrix of training instances (one instance by
  row)}

  \item{y}{numeric vector of values in (-1,+1) representing
  the training labels for each instance in x}

  \item{beta}{a numeric value setting the beta parameter is
  the f-beta score}

  \item{cache}{if NULL (which is the case at the first
  call) parameters values are checked}
}
\value{
a 2 element list (value,gradient) where "value" is the
value of the function at point w, and "gradient" is the
gradient of the loss function at w
}
\description{
F beta score loss function
}
\references{
Teo et al.  A Scalable Modular Convex Solver for
Regularized Risk Minimization.  KDD 2007
}
\seealso{
bmrm
}

