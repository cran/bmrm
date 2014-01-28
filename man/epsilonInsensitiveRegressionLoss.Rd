\name{epsilonInsensitiveRegressionLoss}
\alias{epsilonInsensitiveRegressionLoss}
\title{The loss function to perform a epsilon-insensitive regression (Vapnik et al. 1997)}
\usage{
epsilonInsensitiveRegressionLoss(w, x, y, epsilon, cache = NULL)
}
\arguments{
  \item{w}{weight vector where the function have to be
  evaluated}

  \item{x}{matrix of training instances (one instance by
  row)}

  \item{y}{numeric vector of values representing the
  training labels for each instance in x}

  \item{epsilon}{a numeric value setting tolerance of the
  epsilon-regression}

  \item{cache}{if NULL (which is the case at the first
  call) parameters values are checked}
}
\value{
a 2 element list (value,gradient) where "value" is the
value of the function at point w, and "gradient" is the
gradient of the loss function at w
}
\description{
The loss function to perform a epsilon-insensitive
regression (Vapnik et al. 1997)
}
\references{
Teo et al.  Bundle Methods for Regularized Risk
Minimization JMLR 2010
}
\seealso{
bmrm
}

