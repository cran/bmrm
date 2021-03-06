% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loss_convex.R
\name{ontologyLoss}
\alias{ontologyLoss}
\title{Ontology Loss Function}
\usage{
ontologyLoss(x, y, l = 1 - table(seq_along(y), y),
  dag = diag(nlevels(y)))
}
\arguments{
\item{x}{instance matrix, where x(t,) defines the features of instance t}

\item{y}{target vector where y(t) is an integer encoding target of x(t,)}

\item{l}{loss matrix. l(t,p(t)) must be the loss for predicting target p(t) instead of y(t) 
for instance t. By default, the parameter is set to a 0/1 loss matrix.}

\item{dag}{a numeric matrix defining the path in the Direct Acyclic Graph (DAG) to each class label}
}
\value{
a function taking one argument w and computing the loss value and the gradient at point w
}
\description{
Ontology loss function may be used when the class labels are organized has an ontology structure
}
\examples{
  # -- Load the data
  x <- cbind(intercept=100,data.matrix(iris[1:4]))
  dag <- matrix(nrow=nlevels(iris$Species),byrow=TRUE,dimnames=list(levels(iris$Species)),c(
      1,0,0,0,
      0,1,1,0,
      0,1,0,1
  ))
  w <- nrbm(ontologyLoss(x,iris$Species,dag=dag))
  table(predict(w,x),iris$Species)
}
\references{
Teo et al.
  A Scalable Modular Convex Solver for Regularized Risk Minimization.
  KDD 2007
}
