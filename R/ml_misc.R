


#' Split a dataset for Cross Validation taking into account class balance
#' 
#' @param y the class labels of each sample of the dataset
#' @param num.cv number of cross validation required
#' @return a factor of num.cv levels that assign to each sample a test fold
#' @export
balanced.cv.fold <- function(y,num.cv=10) {
  fold <- factor(rep_len(sample(seq_len(num.cv)),length.out=length(y)))
  i <- sample(seq_along(y))
  o <- order(y[i])
  fold[i[o]] <- fold
  fold
}




#' Compute statistics for ROC curve plotting
#' 
#' @param f decision value for each instance
#' @param y a logical that specify binary labels
#' @return a data.frame() that compute for each threshold value 'f' roc curve statistics: TP, FP, TN, FN, FPR, TPR, sensitivity, specificity, precision, recall, accuracy
#' @author Julien Prados, adapted from Bob Horton code
#' @export
#' @examples
#'   x <- cbind(data.matrix(iris[1:4]))
#'   w <- nrbmL1(rocLoss(x,iris$Species=="versicolor"),LAMBDA=0.01)
#'   with(roc.stat(x %*% w,iris$Species=="versicolor"),plot(FPR,TPR,type="l"))
#'   with(roc.stat(-x[,2],iris$Species=="versicolor"),lines(FPR,TPR,col="blue"))
roc.stat <- function(f,y) {
  if (!is.logical(y)) stop("y must be a logical vector")
  if (length(f)!=length(y)) stop("f and y must have same length")
  o <- order(f, decreasing=TRUE)
  roc <- data.frame(
    f=c(-Inf,f[o]),
    TP = c(0,cumsum(y[o])),
    FP = c(0,cumsum(!y[o]))
  )
  roc <- roc[rev(!duplicated(rev(roc$f))),]
  roc$TN <- sum(!y) - roc$FP
  roc$FN <- sum(y) - roc$TP
  roc$FPR <- roc$FP/sum(!y)
  roc$sensitivity <- roc$recall <- roc$TPR <- roc$TP/sum(y)
  roc$accuracy <- (roc$TP+roc$TN)/length(y)
  roc$specificity <- roc$TNR <- roc$TN / (roc$TN + roc$FP)
  roc$precision <- roc$TP/(roc$TP+roc$FP)
  
  dx <- diff(roc$FPR)
  dy <- diff(roc$TPR)
  attr(roc,"AUC") <- sum(dx*roc$TPR[seq_along(dx)] + dx*dy/2)
  return(roc)
}




