% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_misc.R
\name{balanced.loss.weights}
\alias{balanced.loss.weights}
\title{Compute loss.weights so that total losses of each class is balanced}
\usage{
balanced.loss.weights(y)
}
\arguments{
\item{y}{a object coerced to factor that represent the class labels of each sample of the dataset}
}
\value{
a numeric vector of the same length as y
}
\description{
Compute loss.weights so that total losses of each class is balanced
}
