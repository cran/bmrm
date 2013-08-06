


#' The loss function to perform a least mean square regression
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
lmsRegressionLoss <- function(w,x,y,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y)) stop('y must be a numeric vector')
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    cache <- list()
  }
  
  w <- rep(w,length.out=ncol(x))
  f <- x %*% w
  
  loss <- 0.5*(f-y)^2
  grad <- f-y
  
  # -- convert scalar loss to generic loss
  return(list(
    value=sum(loss),
    gradient=crossprod(x,grad),
    cache=cache
  ))
}



#' The loss function to perform a least absolute deviation regression
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
ladRegressionLoss <- function(w,x,y,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y)) stop('y must be a numeric vector')
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    cache <- list()
  }
  
  w <- rep(w,length.out=ncol(x))
  f <- x %*% w
  
  loss <- abs(f-y)
  grad <- sign(f-y)
  
  # -- convert scalar loss to generic loss
  return(list(
    value=sum(loss),
    gradient=crossprod(x,grad),
    cache=cache
  ))
}



#' The loss function to perform a logistic regression
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
logisticRegressionLoss <- function(w,x,y,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y)) stop('y must be a numeric vector')
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    cache <- list()
  }
  
  w <- rep(w,length.out=ncol(x))
  f <- x %*% w
  
  loss <- log(1+exp(-y*f))
  grad <- -y/(1+exp(-y*f))
  
  # -- convert scalar loss to generic loss
  return(list(
    value=sum(loss),
    gradient=crossprod(x,grad),
    cache=cache
  ))
}



#' The loss function to perform a quantile regression
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param q a numeric value in the range [0-1] defining quantile value to consider
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
quantileRegressionLoss <- function(w,x,y,q=0.5,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y)) stop('y must be a numeric vector')    
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    if (length(q)!=1 || q<0 || q>1) stop('q must be a length one numeric in the range [0-1]')
    cache <- list()
  }
  
  w <- rep(w,length.out=ncol(x))
  f <- x %*% w
  
  loss <- pmax(q*(f-y),(1-q)*(y-f))
  grad <- ifelse(f>y,q,q-1)
  
  # -- convert scalar loss to generic loss
  return(list(
    value=sum(loss),
    gradient=crossprod(x,grad),
    cache = cache
  ))
}



#' The loss function to perform a epsilon-insensitive regression (Vapnik et al. 1997)
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values representing the training labels for each instance in x
#' @param epsilon a numeric value setting tolerance of the epsilon-regression
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
epsilonInsensitiveRegressionLoss <- function(w,x,y,epsilon,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y)) stop('y must be a numeric vector')
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    cache <- list()
  }
  
  w <- rep(w,length.out=ncol(x))
  f <- x %*% w
  
  loss <- pmax(0,abs(f-y)-epsilon)
  grad <- ifelse(abs(f-y)<epsilon,0,sign(f-y))
  
  # -- convert scalar loss to generic loss
  return(list(
    value=sum(loss),
    gradient=crossprod(x,grad),
    cache=cache
  ))
}





