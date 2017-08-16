


#' Non convex loss function for transductive SVM
#' 
#' @param x matrix of training instances (one instance by row)
#' @param y a logical vector representing the training labels for each instance in x. NA are allowed when labels is unknown.
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @seealso nrbm
#' @export
#' @examples
#' x <- cbind(intercept=100,data.matrix(iris[1:2]))
#' y <- iris$Species=="virginica"
#' y[iris$Species=="setosa"] <- NA
#' w <- nrbm(tsvmLoss(x,y),convexRisk=FALSE)
#' table(predict(w,x),iris$Species)
tsvmLoss <- function(x,y,loss.weights=1) {
  if (!is.logical(y)) stop("y must be logical")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  
  function(w) {
    w <- rep(w, length.out = ncol(x))
    f <- x %*% w
    loss <- loss.weights * pmax(ifelse(is.na(y),-abs(f)+1,ifelse(y,-f+1,f+1)),0)
    grad <- loss.weights * (loss>0) * ifelse(is.na(y),-sign(f),ifelse(y,-1,+1))
    lvalue(w) <- sum(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- c("tsvmLoss","binaryClassificationLoss")
    return(w)
  }
}








