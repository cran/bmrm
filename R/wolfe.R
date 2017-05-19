


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
#' @param f0 initial function value
#' @param s0 direction of the search from x0
#' @param ... additional parameters passed to f()
#' @param a1 first step coefficient guess
#' @param amax max coefficient value
#' @param c1 lower bound
#' @param c2 upper bound
#' @param maxiter maximum number of iteration for this linesearch
#' @return the optimal point as a 3 element list
#' @references Do and Artieres
#'   Regularized Bundle Methods for Convex and Non-Convex Risks
#'   JMLR 2012
#' @author Julien Prados
#' @seealso \code{\link{bmrm}} \code{\link{nrbm}}
wolfe.linesearch <- function(f, x0, f0, s0, ..., a1=0.5, amax=1.1, c1=1e-4, c2=0.9, maxiter=5L) {
  gradient(f0) <- crossprod(gradient(f0),s0)
  
  zoom <- function(alo, ahi, flo, fhi) {
    if (alo>ahi) stop("[alo,ahi] is empty")
    for(i in seq_len(maxiter)) {
      # find aj in [alo,ahi] using cubic interpolation
      d1 <- gradient(flo) + gradient(fhi) - 3*(flo-fhi)/(alo-ahi)
      d2 <- sqrt(max(d1*d1 - gradient(flo)*gradient(fhi),0))
      aj <- ahi - (ahi-alo)*(gradient(fhi)+d2-d1)/(gradient(fhi)-gradient(flo)+2*d2)      
      if (aj < alo || aj > ahi) aj <- (alo+ahi)/2
      
      xj <- x0 + aj*s0
      fj <- f(xj, ...)
      gradient(fj) <- crossprod(gradient(fj),s0)
      if (fj > f0 + c1*aj*gradient(f0) || fj > flo) {
        ahi <- aj
        fhi <- fj
      } else {
        if (abs(gradient(fj)) <= -c2*gradient(f0)) break;
        if (gradient(fj)*(ahi-alo) >= 0) {
          ahi <- alo
          fhi <- flo
        }
        alo <- aj
        flo <- fj
      }
      if (abs(alo-ahi) <= 0.01*alo) break;
    }
    if (i==maxiter) warning("iteration stop after maxiter loops")
    list(astar=aj,xstar=xj,fstar=fj)
  }
  
  ai_1 <- 0
  fi_1 <- f0
  for(i in seq_len(maxiter)) {
    ai <- if (i==1L) a1 else (ai + amax)/2
    xi <- x0+ai*s0
    fi <- f(xi,...)
    gradient(fi) <- crossprod(gradient(fi),s0)
    
    if (fi > (f0+c1*ai*gradient(f0)) || (fi >= fi_1 && i > 1)) return(zoom(ai_1, ai, fi_1, fi))
    if (abs(gradient(fi)) <= -c2*gradient(f0)) break
    if (gradient(fi) >= 0) return(zoom(ai, ai_1, fi, fi_1))
    if (abs(ai - amax) <= 0.01*amax) break
   
    # update variables for next iteration
    ai_1 <- ai
    fi_1 <- fi
  }
  if (i==maxiter) warning("iteration stop after maxiter loops") 
  list(astar=ai, xstar=xi, fstar=fi)
}
