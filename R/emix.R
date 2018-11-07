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
#' @param ... additional arguments to be passed to \code{h}
#' 
#' @example examples/dmixEx.R
#' 
emix = function(h, params, wts, ...){
  
  if(!is.matrix(params)) {
    params = matrix(params, ncol=1)
  }
  
  # approximate expectation by summing over mixtures
  h.theta = 0
  for(i in 1:nrow(params)) {
    h.theta = h.theta + h(as.numeric(params[i,]), ...) * wts[i]
  }
    
  h.theta
}
