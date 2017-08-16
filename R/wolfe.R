


#' Wolfe Line Search
#' 
#' Implements Wolfe Line Search algorithm.
#' The code is inspired from Matlab code of Do and Artiere, but not tested.
#' The function is not used yet, but might be used later to speed up bmrm/nrbm 
#' convergence.
#' 
#' @param f a function to minimize. It must accept as first argument a numeric vector
#'   representing the optimization point and return a numeric value, with 
#'   gradient attribute setted
#' @param x0 initial search point
#' @param s0 direction of the search from x0
#' @param ... additional parameters passed to f()
#' @param a1 first step coefficient guess
#' @param amax max coefficient value
#' @param c1 lower bound
#' @param c2 upper bound
#' @param maxiter maximum number of iteration for this linesearch
#' @param f.adjust an adjustment method to adjust lvalue and gradient of f
#' @return the optimal point
#' @references Do and Artieres
#'   Regularized Bundle Methods for Convex and Non-Convex Risks
#'   JMLR 2012
#' @export
#' @author Julien Prados
#' @seealso \code{\link{nrbm}}
wolfe.linesearch <- function(f, x0, s0, ..., a1=0.5, amax=1.1, c1=1e-4, c2=0.9, maxiter=5L,f.adjust=identity) {
  x0 <- f.adjust(x0)
  g0 <- as.vector(crossprod(gradient(x0),s0))
  
  zoom <- function(alo, ahi, flo, fhi, glo, ghi, maxiter) {
    for(j in seq_len(maxiter)) {
      # find aj in [alo,ahi] using cubic interpolation
      d1 <- glo + ghi - 3*(flo-fhi)/(alo-ahi)
      d2 <- sqrt(max(d1*d1 - glo*ghi,0))
      aj <- ahi - (ahi-alo)*(ghi+d2-d1)/(ghi-glo+2*d2)
      
      if (alo>ahi) 
        if (aj < alo || aj > ahi) aj <- (alo+ahi)/2
      else
        if (aj > alo || aj < ahi) aj = (alo+ahi)/2;
      Xj <- f(x0 + aj*s0, ...)
      xj <- f.adjust(Xj)
      gj <- as.vector(crossprod(gradient(xj),s0))
      if (lvalue(xj) > lvalue(x0) + c1*aj*g0 || lvalue(xj) > flo) {
        ahi <- aj
        fhi <- lvalue(xj)
      } else {
        if (abs(gj) <= -c2*g0) break;
        if (gj*(ahi-alo) >= 0) {
          ahi <- alo
          fhi <- flo
        }
        alo <- aj
        flo <- lvalue(xj)
      }
      if (abs(alo-ahi) <= 0.01*alo) break;
    }
    attr(Xj,"neval") <- j
    Xj
  }
  
  
  ai_1 <- 0;fi_1 <- lvalue(x0);gi_1 <- g0
  ai <- a1
  for(i in seq_len(maxiter)) {
    Xi <- f(x0+ai*s0,...)
    xi <- f.adjust(Xi)
    gi <- as.vector(crossprod(gradient(xi),s0))
    
    # test for end of search
    if ((lvalue(xi) > lvalue(x0)+c1*ai*g0) || (lvalue(xi) >= fi_1 && i > 1)) {Xi <- zoom(ai_1, ai, fi_1, lvalue(xi), gi_1, gi, maxiter=maxiter-i+1);break}
    if (abs(gi) <= -c2*g0) break
    if (gi >= 0) {Xi <- zoom(ai, ai_1, lvalue(xi), fi_1, gi, gi_1, maxiter=maxiter-i+1);break}
    if (abs(ai - amax) <= 0.01*amax) break
    
    # update variables for next iteration
    ai_1 <- ai; fi_1 <- lvalue(xi); gi_1 <- gi
    ai <- (ai + amax)/2
  }
  attr(Xi,"neval") <- i + ifelse(is.null(attr(Xi,"neval")),0,attr(Xi,"neval"))
  Xi
}


