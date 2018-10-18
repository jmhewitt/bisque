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
#' @param ... additional arguments to pass to \code{f}
#' 
kCompute = function(f, init, method = 'BFGS', maxit=1e4, level = 2, log = FALSE,
                    ...) {
  
  # find the density's mode
  mode = optim(par = init, fn = f, method = method, 
               control = list(fnscale = -1, maxit=maxit), hessian = TRUE, 
               log = TRUE, ...)
  
  # warn if optimization failed
  if(mode$convergence != 0) {
    warning('Mode not found.')
  }
  
  # build integration grid
  grid = createLocScaleGrid(mu = mode$par, prec = -mode$hessian, level = level)
  
  # evaluate the unnormalized log-density at integration points
  lnf = apply(grid$nodes, 1, function(x){f(x, log = TRUE, ...)})
  lnf = lnf - grid$d

  # initialize return with scaled integration constant
  lnk = -mean(lnf)
  kC = sum(exp(lnf + lnk) * grid$weights)
  lnC = log(kC) - lnk 
  
  if(log) { lnC } else { exp(lnC) }
}