#' Evaluate a mixture density
#'
#' 
#' @export
#' 
#' @param x Points at which the mixture should be evaluated.  If the density 
#'   is multivariate, then each row of \code{x} should contain one set of 
#'   points at which the mixture should be evaluated.
#' @param f Unnormalized density used in the mixture. The function \eqn{f} 
#'   should include a arguments \code{x}, \code{params}, and \code{log}, the
#'   last of which returns \eqn{log(f(x))}.
#' @param pars Matrix in which each row contains parameters for \code{f}.  The 
#'   number of rows in \code{pars} should match the number of mixture
#'   components to evaluate.
#' @param wts vector of weights for each mixture component
#' @param log TRUE to return the log of the mixture density
#' @param ... additional parameters to be passed to \code{f}
#' 
dmix = function(x, f, pars, wts, log = FALSE, ...){
  
  if(!is.matrix(x)) {
    x = matrix(x, ncol=1)
  }
  
  if(is.numeric(pars)) {
    pars = matrix(pars, ncol=1)
  }
  
  res = numeric(nrow(x))
  
  for(i in 1:length(res)) {
    # evaluate mixture components
    lnf = apply(pars, 1, function(params){
      f(x[i,], params = params, log = TRUE, ...)})
    
    # numerically stable evaluation of mixture density on log scale
    lnk = -mean(lnf)
    res[i] = log(sum(exp(lnf + lnk) * wts)) - lnk
  }
  
  if(log) { res } else { exp(res) }
}