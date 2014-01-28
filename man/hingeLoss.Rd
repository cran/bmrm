\name{hingeLoss}
\alias{hingeLoss}
\title{Hinge Loss function for SVM}
\usage{
hingeLoss(w, x, y, cache = NULL, loss.weights = 1)
}
\arguments{
  \item{w}{weight vector where the function have to be
  evaluated}

  \item{x}{matrix of training instances (one instance by
  row)}

  \item{y}{numeric vector of values in (-1,+1) representing
  the training labels for each instance in x}

  \item{cache}{if NULL (which is the case at the first
  call) parameters values are checked}

  \item{loss.weights}{numeric vector of loss weights to
  incure for each instance of x in case of misprediction.
  Vector length should match length(y), but values are
  cycled if not of identical size.  Default to 1 so we
  define a standard 0/1 loss for SVM classifier.  The
  parameter might be useful to adapt SVM learning in case
  of unbalanced class distribution.}
}
\value{
a 2 element list (value,gradient) where "value" is the
value of the function at point w, and "gradient" is the
gradient of the loss function at w
}
\description{
Hinge Loss function for SVM
}
\references{
Teo et al.  A Scalable Modular Convex Solver for
Regularized Risk Minimization.  KDD 2007
}
\seealso{
bmrm
}

