#' Compute expectations via weighted mixtures
#' 
#' Approximates expectations of the form
#' \deqn{E[h(\theta)] = \int h(\theta) f(\theta) d\theta}
#' using a weighted mixture
#' \deqn{E[h(\theta)] \approx \sum_{j=1}^k h(\theta^{(k)}) w_k}
#' 
#' 
#' @export
#' 
#' @param h Function for which the expectation should be taken.  The function 
#'   should be defined so it is can be called via \code{f(params, ...)}.
#'   Additional parameters may be passed to \eqn{h} via \code{...}.
#' @param params Matrix in which each row contains parameters at which
#'   \eqn{h} should be evaluated.  The number of rows in \code{params} should 
#'   match the number of mixture components \eqn{k}.
#' @param wts vector of weights for each mixture component
#' @param errorNodesWts list with elements \code{inds} and \code{weights} that 
#'   point out which \code{params} get used to compute an approximation of the 
#'   quadrature error.
#' @param ... additional arguments to be passed to \code{h}
#' 
#' 
emix = function(h, params, wts, errorNodesWts = NULL, ...){
  
  if(!is.matrix(params)) {
    params = matrix(params, ncol=1)
  }
  
  # initialize posterior mean
  h.theta = 0
  
  # initialize state for quadrature error bound
  if(!is.null(errorNodesWts)) { 
    h.theta.l = 0
    err.ind = 1
    next.err.ind = errorNodesWts$inds[err.ind]
  }
  
  # approximate expectation by summing over mixtures
  # TODO: parallelize step
  for(i in 1:nrow(params)) {
    
    # build posterior mean estimate
    h.eval = h(as.numeric(params[i,]), ...)
    h.theta = h.theta + h.eval * wts[i]
    
    # build quadrature error bound
    if(!is.null(errorNodesWts)) {
      if(i == next.err.ind) {
        
        h.theta.l = h.theta.l + h.eval * errorNodesWts$weights[err.ind]
        
        err.ind = err.ind + 1
        next.err.ind = ifelse(err.ind <= length(errorNodesWts$inds), 
                              errorNodesWts$inds[err.ind],
                              -1)
      }
    }
    
  }
    
  if(!is.null(errorNodesWts)) {
    list( E = h.theta,
          E.coarse = h.theta.l,
          rel.err.bound = (h.theta.l - h.theta)/h.theta * 100 )
  } else {
    h.theta
  }
}
