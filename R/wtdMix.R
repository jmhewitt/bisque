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
#'     \item{param}{a vector of parameters needed to evaluate
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
#' @param f1.precompute function that pre-computes parameters for evaluating
#'   f(theta1 | X, theta2).  param must be defined as param = function(theta2).
#' @param w TRUE if the last element in the vector returned by the param
#'   function should be added to the log-weights before they are normalized
#'   before they are normalized.  should be a function of theta2
#' @param f2.link character vector that specifies transformations used during
#'   optimization and integration of f(theta2 | X).  while theta2 may be
#'   defined on arbitrary support, \code{wtdMix} performs optimization and
#'   integration of theta2 on an unconstrained support.  the \code{link}
#'   vector describes the transformations that must be applied to each
#'   element of theta2.  Jacobian functions for the transformations will
#'   automatically be added to the optimization and integration routines.
#'   currently supported link functions are 'log', 'logit', and 'identity'.
#' @param level accuracy level (typically number of grid points for the
#'   underlying 1D quadrature rule) [description from mvQuad::createNIGrid]
#' @param w.init initial guess for mode of f(theta2 | X)
#' @param control limited list of options to pass to optim.  should be a list
#'   of lists:
#'   \describe{
#'     \item{w}{(default) list(method='BFGS', maxit=5e4)
#'     }
#'   }
#' @param ... additional arguments to pass to f1, f2, param
#'
wtdMix = function(f1, f1.precompute = function(x, ...){x}, f2,
                  w = 'direct', w.init, w.link = NULL, level = 2,
                  w.control = NULL, ...) {

  #
  # verify we can compute a gaussian approximation to f(theta2 | X)
  #

  w.fail = FALSE
  if(is.character(w)) {
    if(w == 'direct') {
      # assumption: w is exactly f(theta2 | X) or proportional to f(theta2 | X)
      w.approx = TRUE
    }
  } else if(is.list(w)) {
    # assumption: f(theta2 | X) will be approximated from f(theta1, theta2 | X)
    w.approx = FALSE

    # verify joint posterior is available
    if(!is.function(w$f12)) { w.fail = TRUE }

    # verify indices for marginalizing f^G(theta1, theta2 | X) is provided
    if(is.null(w$theta2.inds)) { w.fail = TRUE }

    # verify or initialize information for computing C1(theta2)
    if(is.null(w$f1.init)) { w.fail = TRUE }
    if(w.fail == FALSE) {
      # default is identity links
      if(is.null(w$f1.link)) {
        w$f1.link = rep('identity', length(w$f1.init))
      }
      # default optim options
      if(is.null(w$f1.control)) {
        w$f1.control = list(method = 'BFGS', maxit = 5e4)
      }
      # default integration level
      if(is.null(w$f1.level)) { w$f1.level = 2 }
    }
  }

  if(w.fail) {
    stop('Gaussian approximation to f(theta2 | X) cannot be
          computed or used with the information passed in the argument w.')
  }


  #
  # initialize gaussian approximation for f(theta2 | X)
  #

  # define g as function to optimize
  if(w.approx) { g = f2 }
  else { g = w$f12 }
  
  # default is identity links
  if(is.null(w.link)) { w.link = rep('identity', length(w.init)) }

  # default optim options
  if(is.null(w.control)) { w.control = list(method = 'BFGS', maxit = 5e4) }


  #
  # develop gaussian approximation to f(theta2 | X) on transformed scale
  #

  # find mode of g
  g.mode = optim(par = tx(w.init,  w.link), fn = function(par) {
    g(itx(par, w.link), log = TRUE, ...) + sum(logjac(par, w.link))
  }, method = w.control$method, hessian = TRUE,
  control = list(fnscale = -1,  maxit = w.control$maxit))

  # check convergence
  if(g.mode$convergence != 0) { warning('Mode for w may not have been found.') }
  
  # extract gaussian approximation for f(theta2 | X) on transformed scale
  if(w.approx) {
    f2.mode = g.mode$par
    f2.prec = - g.mode$hessian
    f2.link = w.link
  } else {
    f2.mode = g.mode$par[w$theta2.inds]
    L1 = matrix(- g.mode$hessian[-w$theta2.inds, -w$theta2.inds],
                nrow = nrow(g.mode$hessian) - length(w$theta2.inds))
    L2 = matrix(- g.mode$hessian[w$theta2.inds, w$theta2.inds],
                nrow = length(w$theta2.inds))
    L12 = matrix(- g.mode$hessian[-w$theta2.inds, w$theta2.inds],
                 ncol = length(w$theta2.inds))
    f2.prec =  L2 - t(L12) %*% solve(L1) %*% L12
    f2.link = w.link[w$theta2.inds]
  }


  #
  # develop weighted mixtures approximation
  #
  
  # create integration grid
  grid = createLocScaleGrid(mu = f2.mode, prec = f2.prec, level = level)

  # preallocate space for C1(theta2), as necessary
  if(w.approx == FALSE) { C1 = numeric(nrow(grid$nodes)) }

  # precompute parameters and weights for posterior mixture for theta1
  # TODO: don't double compute p0
  p0 = f1.precompute(itx(grid$nodes[1,], f2.link), ...)
  mix = matrix(NA, nrow = nrow(grid$nodes), ncol = length(p0))
  wts = numeric(nrow(mix))
  for(i in 1:nrow(mix)) {
    theta2 = as.numeric(itx(grid$nodes[i,], f2.link))
    mix[i,] = f1.precompute(theta2, ...)
    wts[i] = f2(theta2, log = TRUE, ...) +
      sum(logjac(grid$nodes[i,], f2.link)) - grid$d[i]

    # compute C1(theta2) and update weights
    if(w.approx == FALSE) {
      C1[i] = kCompute(f = function(theta1, log = TRUE) {
        res = f1(theta1, theta2, log = log, ...)
        if(log) { res } else { exp(res) }
      }, init = w$f1.init, log = TRUE, level = w$f1.level, link = w$f1.link)
      wts[i] = wts[i] + C1[i]
    }
  }

  # normalize weights
  wts = exp(wts - mean(wts)) * grid$weights
  wts = wts / sum(wts)

  # use C1(theta2)  to build a second-layer function for dmix
  if(w.approx) { h = f1 }
  else { 
  # TODO: Make it safer to pass C1 to the dmix function
    mix = cbind(mix, C1)
    h = function(theta1, params, log, ...) {
      res = f1(theta1, params, log, ...) - params[length(params)]
      if(log) { res } else { exp(res) }
    }
  }

  # build and return weighted marginal posterior
  list(
    f = function(theta1, log = FALSE) {
      dmix(x = theta1, f = h, params = mix, wts = wts, log = log, ...)
    },
    mix = mix,
    wts = wts
  )
}
