

#' Compute or check the structure of a cost matrix 
#' 
#' @param y integer vector of positive values (>=1) representing the labels of the instances
#' @param C either a cost matrix to check for consistency with labels in y, or a character string defining the standard matrix to compute. 
#'        If a character string the accepted values are "0/1" for a 0-1 cost matrix or "linear" for linear cost.
#' @return the cost matrix object
#' @export
#' @seealso bmrm, ordinalRegressionLoss
costMatrix <- function(y,C=c("0/1","linear")) {
  if (!is.integer(y) || min(y)!=1) stop('y must be an integer vector of positive values representing the labels')  
  if (is.character(C)) {
    C <- match.arg(C)
    C <- switch(C,
           "0/1" = local({
              C <- matrix(1,max(y),max(y))
              diag(C) <- 0
              return(C)
           }),
           "linear" = abs(outer(1:max(y),1:max(y),'-'))
    )
  } else {
    C <- as.matrix(C)
    if (nrow(C)!=ncol(C)) stop("C must be a square matrix") 
    if (max(y)>nrow(C)) stop("one value in y is out of range for dimensions of C")
  }
  return(C)
}

