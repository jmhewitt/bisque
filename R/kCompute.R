#' Compute integration constants for unnormalized densities
#'
#' 
#' @export
#' 
#' @param f Unnormalized density for which to compute a scaled integration
#'   constant.  the function \eqn{f} should include an argument \code{log}, 
#'   which returns \eqn{log(f(x))}.
#' @param init Initial guess for the density's mode
#' @param maxit maximum number of iterations \code{optim} should use in 
#'   searching for the density's mode
#' @param method method to be used to search for the density's mode
#' @param level accuracy level (typically number of grid points for the 
#'   underlying 1D quadrature rule) [description from mvQuad::createNIGrid]
#' @param log TRUE to return log of integration constant
#' @param link character vector that specifies transformations used during 
#'   optimization and integration of f(theta2 | X).  while theta2 may be 
#'   defined on arbitrary support, \code{wtdMix} performs optimization and 
#'   integration of theta2 on an unconstrained support.  the \code{link} 
#'   vector describes the transformations that must be applied to each 
#'   element of theta2.  Jacobian functions for the transformations will 
#'   automatically be added to the optimization and integration routines.
#'   currently supported link functions are 'log', 'logit', and 'identity'.
#' @param ... additional arguments to pass to \code{f}
#' 
#' @examples
#' kCompute(dgamma, init = 1, shape=2, link='log', level = 5)
#' 
kCompute = function(f, init, method = 'BFGS', maxit=1e4, level = 2, log = FALSE,
                    link = NULL, ...) {
  
  # default is identity links
  if(is.null(link)) {
    link = rep('identity', length(init))
  }
  
  # function to apply link transformations
  tx = function(x) {
    for(i in 1:length(link)) {
      x[i] = switch (link[i],
                     'identity' = x[i],
                     'log' = log(x[i]),
                     'logit' = qlogis(x[i])
      )
    }
    x
  }
  
  # function to invert link transformations
  itx = function(x) {
    for(i in 1:length(link)) {
      x[i] = switch (link[i],
                     'identity' = x[i],
                     'log' = exp(x[i]),
                     'logit' = plogis(x[i])
      )
    }
    x
  }
  
  # function to apply jacobians for link transformations
  logjac = function(x) {
    for(i in 1:length(link)) {
      x[i] = switch (link[i],
                     'identity' = 0,
                     'log' = jac.log(x[i], log = TRUE),
                     'logit' = jac.logit(x[i], log = TRUE)
      )
    }
    x
  }
  
  # find the density's mode
  mode = optim(par = tx(init), fn = function(par, ...) {
    f(itx(par), log = TRUE, ...) + sum(logjac(par))
  }, method = method, control = list(fnscale = -1, maxit=maxit), 
  hessian = TRUE, ...)
  
  # warn if optimization failed
  if(mode$convergence != 0) {
    warning('Mode not found.')
  }
  
  # build integration grid
  grid = createLocScaleGrid(mu = mode$par, prec = -mode$hessian, level = level)
  
  # evaluate the unnormalized log-density at integration points
  lnf = apply(grid$nodes, 1, function(x){
    f(itx(x), log = TRUE, ...) + sum(logjac(x)) })
  lnf = lnf - grid$d

  # initialize return with scaled integration constant
  lnk = -mean(lnf)
  kC = sum(exp(lnf + lnk) * grid$weights)
  lnC = log(kC) - lnk
  
  if(log) { lnC } else { exp(lnC) }
}