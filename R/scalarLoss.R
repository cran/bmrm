

#' Hinge Loss function for SVM
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values in (-1,+1) representing the training labels for each instance in x
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @seealso bmrm
hingeLoss <- function(w,x,y,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y) || !all(y %in% c(-1,1))) stop('y must be a numeric vector of either -1 or +1')
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    cache <- list()
  }

  w <- rep(w,length.out=ncol(x))
  f <- x %*% w
  
  loss <- pmax(0,1-y*f)
  grad <- (loss>0) * (-y)
  
  # -- convert scalar loss to generic loss
  return(list(
    value=sum(loss),
    gradient=crossprod(x,grad),
    cache = cache
  ))
}



#' The loss function to maximize area under the ROC curve
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values in (-1,+1) representing the training labels for each instance in x
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
rocLoss <- function(w,x,y,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y) || !all(y %in% c(-1,1))) stop('y must be a numeric vector of either -1 or +1')    
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    cache <- list()
  }
  
  w <- rep(w,length.out=ncol(x))
  c <- x %*% w - 0.5*y
  o <- order(c)
  
  sp <- cumsum(y[o]==+1)
  sm <- sum(y==-1) - cumsum(y[o]==-1)
  l <- numeric(length(o))
  l[o] <- ifelse(y[o]==-1,sp,-sm)
  l <- l/(sum(y==-1)*sum(y==+1))
  
  # -- convert scalar loss to generic loss
  return(list(
    value = crossprod(l,c),
    gradient = crossprod(l,x),
    cache = cache
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
    if (!is.numeric(y)) stop('y must be a numeric vector of either -1 or +1')    
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




#' The loss function to perform a mean least square regression
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
mlsRegressionLoss <- function(w,x,y,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y)) stop('y must be a numeric vector of either -1 or +1')    
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



#' The loss function for ordinal regression
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y integer vector of positive values (>=1) representing the training labels for each instance in x
#' @param C the cost matrix to use, C[i,j] being the cost for predicting label i instead of label j.
#' @param impl either the string "loglin" or "quadratic", that define the implementation to use for the computation of the loss.
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   Bundle Methods for Regularized Risk Minimization
#'   JMLR 2010
#' @seealso bmrm
#' @examples
#' # -- Load the data
#' x <- data.matrix(iris[1:4])
#' y <- as.integer(iris$Species)
#' 
#' # -- Train the model
#' m <- bmrm(x,y,lossfun=ordinalRegressionLoss,LAMBDA=0.001,EPSILON_TOL=0.0001)
#' m2 <- bmrm(x,y,impl="quadratic",lossfun=ordinalRegressionLoss,LAMBDA=0.001,EPSILON_TOL=0.0001)
#' 
#' # -- plot predictions
#' f <- x %*% m$w
#' f2 <- x %*% m2$w
#' layout(1:2)
#' plot(y,f)
#' plot(f,f2,main="compare predictions of quadratic and loglin implementations")
#' 
#' # -- Compute accuracy
#' ij <- expand.grid(i=seq(nrow(x)),j=seq(nrow(x)))
#' n <- tapply(f[ij$i] - f[ij$j]>0,list(y[ij$i],y[ij$j]),sum)
#' N <- table(y[ij$i],y[ij$j])
#' print(n/N)
ordinalRegressionLoss <- function(w,x,y,C="0/1",impl=c("loglin","quadratic"),cache=NULL) {
  impl <- match.arg(impl)

  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    cache <- within(list(),{
      m <- length(y)
      mi <- tabulate(y)
      M <- (m*m - sum(mi*mi))/2
      C <- costMatrix(y,C) / M
    })
  }
  
  .loglin <- function(w,x,y,cache) {
    with(cache,{
      w <- rep(w,length.out=ncol(x))  
      f <- x %*% w
      c <- c(f-0.5,f+0.5)
      o <- order(c)
      
      j <- ((o-1)%%m)+1
      
      l <- matrix(0,2*m,length(mi))
      l[cbind(which(o<=m),y[j[o<=m]])] <- 1
      l <- apply(l,2,cumsum)
      
      u <- matrix(0,2*m,length(mi))
      u[cbind(which(o>m),y[j[o>m]])] <- 1
      u <- mi[col(u)] - apply(u,2,cumsum)
      
      Gu <- t(C)[y[j],] * u
      Gu[col(Gu)>=y[j]] <- 0
      Gl <- C[y[j],] * l
      Gl[col(Gl)<=y[j]] <- 0
      
      v <- ifelse(o<=m,-rowSums(Gu), rowSums(Gl))
      r <- sum(v*c[o])
      g <- matrix(NA,m,2)
      g[cbind(j,1 + (o-1)%/%m)] <- v
      g <- rowSums(g)
    
      # -- convert scalar loss to generic loss
      list(
        value = r,
        gradient = crossprod(g,x),
        cache = cache
      )
    })
  }
  
  .quadratic <- function(w,x,y,cache) {
    with(cache,{
      w <- rep(w,length.out=ncol(x)) 
      f <- x %*% w
      
      # alternative computation in quadratic time for debugging purpose only
      z <- expand.grid(i=factor(1:m),j=factor(1:m))
      z <- z[y[z$i] < y[z$j],]
      z <- z[1+f[z$i]-f[z$j]>0,]
      R <- sum(C[cbind(y[z$i],y[z$j])] * (1+f[z$i]-f[z$j]))
      G <- colSums(C[cbind(y[z$i],y[z$j])] * (x[z$i,]-x[z$j,]))
      
      list(value = R,gradient = G,cache = cache)
    })
  }
  
  impl.fun <- switch(impl,loglin=.loglin,quadratic=.quadratic)
  impl.fun(w,x,y,cache)
}




#' F beta score loss function
#' 
#' @param w weight vector where the function have to be evaluated
#' @param x matrix of training instances (one instance by row)
#' @param y numeric vector of values in (-1,+1) representing the training labels for each instance in x
#' @param beta a numeric value setting the beta parameter is the f-beta score
#' @param cache if NULL (which is the case at the first call) parameters values are checked
#' @return a 2 element list (value,gradient) where "value" is the value of the function at point w, and "gradient" is the gradient of the loss function at w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @seealso bmrm
fbetaLoss <- function(w,x,y,beta=1,cache=NULL) {
  
  # check parameters at first call
  if (is.null(cache)) {
    if (!is.matrix(x)) stop('x must be a numeric matrix')
    if (!is.numeric(y) || !all(y %in% c(-1,1))) stop('y must be a numeric vector of either -1 or +1')
    if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
    cache <- list()
  }
  
  .fbeta <- function(TP,TN,P,N,beta) {
    beta2 <- beta*beta
    (1+beta2)*TP / (TP+N-TN+beta2*P)
  }
  
  w <- rep(w,length.out=ncol(x))
  f <- x %*% w
  o <- order(f,decreasing=FALSE)
  op <- o[y[o]==1]
  on <- rev(o[y[o]==-1])
  
  p <- 2*(sum(f[op]) - cumsum(c(0,f[op])))
  n <- 2*(sum(f[on]) - cumsum(c(0,f[on])))
  R <- outer(seq_along(p),seq_along(n),function(i,j) {
    1 - .fbeta(i-1,j-1,length(op),length(on),beta) - p[i] + n[j]
  })

  ij <- arrayInd(which.max(R),dim(R))
  Y <- -y
  Y[op[seq(length=ij[1,1]-1)]] <- 1
  Y[on[seq(length=ij[1,2]-1)]] <- -1
    
  # -- convert scalar loss to generic loss
  return(list(
    value=R[ij],
    gradient=crossprod(x,Y-y),
    cache = cache
  ))
}


