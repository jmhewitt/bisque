#' Assemble a weighted mixture of posteriors
#'
#' 
#' @export
#' 
#' @param f1 evaluates f(theta1 | X, theta2).  f1 must be defined as 
#'   f1 = function(theta1, params, log, ...).  
#'   \describe{
#'     \item{theta1}{a matrix of parameters at which to evaluate 
#'       f(theta1 | X, theta2). each row should be one set of values at which 
#'       the density should be evaluated}
#'     \item{par}{a vector of parameters needed to evaluate 
#'       f(theta1 | X, theta2).  In most cases par will equal theta2, but in 
#'       some cases, f(theta1 | X, theta2) depends on functions of theta2, 
#'       which can be precomputed (to save time) as the weighted mixture 
#'       approximation is constructed.}
#'     \item{log}{TRUE to return log(f(theta1 | X, theta2))}
#'     \item{...}{additional arguments needed for function evaluation}
#'   }
#' @param f2 evaluates f(theta2 | X).  f2 must be defined as 
#'   f2 = function(theta2, log, ...).  The interpretation of these arguments
#'   is the same as for f1
#' @param param function that pre-computes parameters for evaluating 
#'   f(theta1 | X, theta2).  param must be defined as param = function(theta2).
#' @param w.mod TRUE if the last element in the vector returned by the param
#'   function should be added to the log-weights before they are normalized
#'   before they are normalized.  should be a function of theta2
#' @param link character vector that specifies transformations used during 
#'   optimization and integration of f(theta2 | X).  while theta2 may be 
#'   defined on arbitrary support, \code{wtdMix} performs optimization and 
#'   integration of theta2 on an unconstrained support.  the \code{link} 
#'   vector describes the transformations that must be applied to each 
#'   element of theta2.  Jacobian functions for the transformations will 
#'   automatically be added to the optimization and integration routines.
#'   currently supported link functions are 'log', 'logit', and 'identity'.
#' @param level accuracy level (typically number of grid points for the 
#'   underlying 1D quadrature rule) [description from mvQuad::createNIGrid]
#' @param f2.init initial guess for mode of f(theta2 | X)
#' @param control limited list of options to pass to optim.  should be a list 
#'   of lists:
#'   \describe{
#'     \item{f2}{(default) list(method='BFGS', maxit=5e4)
#'     }
#'   }
#' 
wtdMix = function(f1, f2, f2.init, param, w.mod = FALSE, link = NULL, level = 2,
                  control = NULL) {
  
  # default is identity links
  if(is.null(link)) {
    link = rep('identity', length(f2.init))
  }
  
  # default optim options
  if(is.null(control)) { control = list() }
  if(is.null(control$f2)) {
    control$f2 = list( method = 'BFGS', maxit = 5e4)
  }
  
  # function to apply link transformations
  tx = function(x) {
    for(i in 1:length(link)) {
      x[i] = switch (link[i],
                     'identity' = x,
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
                     'identity' = x,
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
  
  # find posterior mode for f(theta2 | X) on transformed scale
  f2.mode = optim(par = tx(f2.init), fn = function(par) {
    f2(itx(par), log = TRUE) + sum(logjac(par))
  }, method = control$f2$method, 
     control = list(fnscale = -1, maxit = control$f2$maxit), hessian = TRUE)
  
  if(f2.mode$convergence != 0) {
    warning('Posterior mode for f(theta2 | X) may not have been found.')
  }
  
  # create integration grid
  grid = createLocScaleGrid(mu = f2.mode$par, prec = -f2.mode$hessian)
  
  # precompute parameters and weights for posterior mixture for theta1
  # TODO: don't double compute p0
  # TODO: allow some of the internal state to be used to evaluate param(.) ?
  p0 = param(itx(grid$nodes[1,]))
  mix = matrix(NA, nrow = nrow(grid$nodes), ncol = length(p0))
  wts = numeric(nrow(mix))
  for(i in 1:nrow(mix)) {
    theta2 = as.numeric(itx(grid$nodes[i,]))
    mix[i,] = param(theta2)
    wts[i] = f2(theta2, log = TRUE) + sum(logjac(grid$nodes[i,])) - grid$d[i]
    if(w.mod == TRUE) { 
      wts[i] = wts[i] + mix[i,ncol(mix)]
    }
  }
  
  # normalize weights
  wts = exp(wts - mean(wts)) * grid$weights
  wts = wts / sum(wts)
  
  # build and return weighted marginal posterior
  list(
    f = function(theta1, log = FALSE) {
      dmix(x = theta1, f = f1, params = mix, wts = wts, log = log)
    },
    mix = mix,
    wts = wts
  )
}