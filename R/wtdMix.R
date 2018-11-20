#' Develop a weighted mixture of posteriors
#'
#' For a Bayesian model
#' \deqn{ X ~ f(X | \theta_1, \theta_2)}
#' \deqn{ (\theta_1, \theta_2) ~ f(\theta_1, \theta_2),}
#' the marginal  posterior \eqn{f(\theta_1 | X)} distribution can be
#' approximated via weighted mixtures via
#' \deqn{ f(\theta_1 | X) \approx \sum_{j=1}^K f(\theta_1 | X, \theta_2) w_j }
#' where \eqn{w_j} is based on \eqn{f(\theta_2^{(j)} | X)} and weights
#' \eqn{\tilde w_j}, where \eqn{\theta_2^{(j)}} and \eqn{\tilde w_j} are
#' nodes and weights for a sparse-grid quadrature integration scheme.
#' The quadrature rule is developed by finding the posterior mode of
#' \eqn{f(\theta_2|X)}, after transforming \eqn{\theta_2} to an unconstrained
#' support.  For best results, \eqn{\theta_2} should be a continuous random
#' variable, or be able to be approximated by one.
#'
#' @import foreach
#' @importFrom itertools ichunk
#'
#' @export
#'
#' @param f1 evaluates \eqn{f(\theta_1 | X, \theta_2)}.  \code{f1} must be able
#'   to be called via \code{f1(theta1, params, log, ...)}.
#'   \describe{
#'     \item{\code{theta1}}{a matrix of parameters at which to evaluate
#'       \eqn{f(\theta_1 | X, \theta_2)}. each row should be one set of values
#'       at which the density should be evaluated}
#'     \item{params}{a vector of parameters needed to evaluate
#'       \eqn{f(\theta_1 | X, \theta_2)}.  In most cases \code{params} will
#'       equal \eqn{theta_2}, but in some cases, \eqn{f(\theta_1 | X, \theta_2)}
#'       depends on functions of \eqn{\theta_2}, which can be pre-evaluated
#'       as the weighted mixture approximation is constructed.}
#'     \item{log}{TRUE to return \eqn{ln(f(\theta_1 | X, \theta_2))}}
#'     \item{...}{additional arguments needed for function evaluation}
#'   }
#' @param f1.precompute function that pre-computes parameters for evaluating
#'   \eqn{f(\theta_1 | X, \theta_2)}.  \code{f1.precompute} must be able to
#'   be called via \code{f1.precompute(theta2, ...)} and return the argument
#'   \code{params} for the function \code{f1}.
#' @param f2 evaluates \eqn{f(theta_2 | X)}.  \code{f2} must be able to be
#'   called via \code{f2(theta2, log, ...)}.
#' @param w Argument that specifies how the posterior mode of
#'   \eqn{f(\theta_2| X)} should be evaluated.  The options for \code{w} are
#'   described below.
#'   \describe{
#'     \item{\code{'direct'}}{Instructs \code{wtdMix} to optimize \eqn{f2}
#'       directly, after applying link function transformations to the
#'       parameters \eqn{\theta_2}.}
#'     \item{\code{list()}}{Instructs \code{wtdMix} to determine the posterior
#'       mode of \eqn{f(\theta_2 | X)} using a function that is proportional,
#'       or otherwise a major component of \eqn{f(\theta_2 | X)}.  The list
#'       must specify the following arguments:
#'       \describe{
#'         \item{\code{f12}}{Function that is informative of
#'           \eqn{f(\theta_2 | X)} that will be optimized.  For example,
#'           \code{f12} may include components of \eqn{f(\theta_2 | X)}, or it
#'           may be proportional to the joint posterior distribution
#'           \eqn{f(\theta_1, \theta_2 | X)}.  \code{f12} must be able to be
#'           called via \code{f12(theta12, log, ...)}.}
#'         \item{\code{theta2.inds}}{After optimizing \code{f12}, \code{wtdMix}
#'           will use the hessian at the mode to develop a gaussian
#'           approximation to \code{f12}.  The argument \code{theta2.inds}
#'           specifies the indices of \code{theta12} that correspond to
#'           \eqn{\theta_2}.  \code{wtdMix} will integrate out the
#'           non-\eqn{\theta_2} indices from the Gaussian approximation to
#'           \code{f12}.}
#'         \item{\code{f1.init}}{It is assumed that
#'           \eqn{f(\theta_1|\theta_2, X)} is not known exactly if \code{f12}
#'           is used.  \code{wtdMix} will determine an integration constant for
#'           \eqn{f(\theta_1|\theta_2, X)} in this case by using numerical
#'           integration around the mode of the function.  \code{f1.init}
#'           specifies the initial guess for the mode.}
#'         \item{\code{f1.link}}{The integration of \code{f1} will be done on
#'           a transformed scale.  \code{f1.link} is a character vector that
#'           specifies the transformations to use.}
#'          \item{\code{f1.control}}{Passes arguments to \code{optim} to
#'           determine the optimization scheme used.}
#'          \item{\code{f1.level}}{Determines the quality of the approximate
#'           integration.  Higher numbers are more accurate approximations.}
#'       }
#'      }
#'   }
#' @param w.link character vector that specifies transformations used during
#'   optimization and integration of \eqn{f(\theta_2 | X)}.  While
#'   \eqn{\theta_2} may be defined on arbitrary support, \code{wtdMix} performs
#'   optimization and integration of \eqn{\theta_2} on an unconstrained support.
#'   The \code{link} vector describes the transformations that must be applied
#'   to each element of \eqn{\theta_2}.  Jacobian functions for the
#'   transformations will automatically be added to the optimization and
#'   integration routines. Currently supported link functions are \code{'log'},
#'   \code{'logit'}, and \code{'identity'}.
#' @param w.init initial guess for mode of \eqn{f(\theta_2 | X)}.  Default is
#'   \code{'identity'} for all parameters.
#' @param level accuracy level (typically number of grid points for the
#'   underlying 1D quadrature rule) [description from mvQuad::createNIGrid]
#' @param w.control Limited list of options to pass to optim.
#' @param ... Additional arguments to pass to \code{f1}, \code{f1.precompute},
#'   \code{f12}, and \code{f2}.
#' @param ncores number of cores used to parallelize computation of parameters
#'   for \eqn{f(\theta_1 | \theta_2, X)}.
#'
#' @return A list with class \code{wtdMix}, which contains the following items.
#'   \describe{
#'     \item{\code{f}}{Function for evaluating the posterior density
#'      \eqn{f(\theta_1|X)}.  \code{f} is callable  via
#'      \code{f(theta1, log, ...)}.}
#'     \item{\code{mix}}{A matrix containing the pre-computed parameters for
#'       evaluating the mixture components \eqn{f(\theta_1 | \theta_2, X)}.
#'       Each row of the matrix contains parameters for one of the \eqn{K}
#'       mixture components.}
#'     \item{\code{wts}}{Integration weights for each of the mixture components.
#'       Some of the weights may be negative.}
#'     \item{\code{expectation}}{List containing additional tools for computing
#'       posterior expectations of \eqn{f(\theta_2|X)}.  However, posterior 
#'       expectations of \eqn{f(\theta_1|X)} can also be computed when 
#'       expectations of \eqn{f(\theta_1|\theta_2, X)} are known.  The elements
#'       of \code{expectation} are
#'       \describe{
#'         \item{\code{Eh}}{Function to compute \eqn{E[h(\theta_2)|X]}.  
#'           \code{Eh} is callable via \code{Eh(h, ...)}, where \code{h} is a
#'           function callable via \code{h(theta2, ...)} and \code{...} are 
#'           additional arguments to the function.  The function \code{h} is  
#'           evaluated at the quadrature nodes \eqn{\theta_2^{(j)}}.}
#'         \item{\code{Eh.precompute}}{Exactly the same idea as \code{Eh}, but
#'           the function \code{h} is evalauted at the quadrature nodes after 
#'           being passed through the function \code{f1.precompute}.}
#'         \item{\code{grid}}{The sparse-quadrature integration grid used.  
#'           Helpful for seeing the quadrature nodes \eqn{\theta_2^{(j)}}.}
#'         \item{\code{wts}}{The integration weights for approximating the 
#'           expectation \eqn{E[h]}.  Note that these integration weights may
#'           differ from the main integration weights for evaluating the 
#'           posterior density \eqn{f(\theta_1|X)}.}
#'       }}
#'   }
#'
wtdMix = function(f1, f1.precompute = function(x, ...){x}, f2,
                  w = 'direct', w.init, w.link = NULL, level = 2,
                  w.control = NULL, ncores = 1, ...) {

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
    if(!is.na(w$theta2.inds)) {
      f2.mode = g.mode$par[w$theta2.inds]
      L1 = matrix(- g.mode$hessian[-w$theta2.inds, -w$theta2.inds],
                  nrow = nrow(g.mode$hessian) - length(w$theta2.inds))
      L2 = matrix(- g.mode$hessian[w$theta2.inds, w$theta2.inds],
                  nrow = length(w$theta2.inds))
      L12 = matrix(- g.mode$hessian[-w$theta2.inds, w$theta2.inds],
                   ncol = length(w$theta2.inds))
      Delta = forwardsolve(t(chol(L1)), L12)
      f2.prec = L2 - t(Delta) %*% Delta
      # f2.prec =  L2 - t(L12) %*% solve(L1) %*% L12
      f2.link = w.link[w$theta2.inds]
    } else {
      f2.mode = g.mode$par
      f2.prec = - g.mode$hessian
      f2.link = w.link
    }
  }

  #
  # develop weighted mixtures approximation
  #

  # create integration grid
  grid = createLocScaleGrid(mu = f2.mode, prec = f2.prec, level = level)

  # preallocate space for C1(theta2), as necessary
  if(w.approx == FALSE) { C1 = numeric(nrow(grid$nodes)) }

  # precompute parameters and weights for posterior mixture for theta1
  p0 = f1.precompute(itx(grid$nodes[1,], f2.link), ...)
  nodes  = nrow(grid$nodes)
  chunkSize = ceiling(nodes/ncores)
  pc = foreach(inds = ichunk(1:nodes, chunkSize = chunkSize, mode = 'numeric'),
               .combine = mergePars, 
               .export = c('itx', 'logjac', 'kCompute')) %dopar% {
                             
    # initialize return objects
    mix = matrix(NA, nrow = length(inds), ncol = length(p0))
    wts = numeric(nrow(mix))
    wts.e = numeric(nrow(mix))
    C1 = numeric(nrow(mix))
    nodes.backtransformed = grid$nodes[inds,]

    for(i in 1:nrow(mix)) {
      # back-transform parameters
      theta2 = as.numeric(itx(grid$nodes[inds[i],], f2.link))

      # compute mixture parameters
      mix[i,] = f1.precompute(theta2, ...)

      # compute base weights
      wts[i] = f2(theta2, log = TRUE, ...) +
        sum(logjac(grid$nodes[inds[i],], f2.link)) - grid$d[inds[i]]

      # compute weights for evaluating expectations in secondary analyses
      wts.e[i] = wts[i]
      
      # as necessary, compute C1(theta2) and update weights
      if(w.approx == FALSE) {
        C1[i] = kCompute(f = function(theta1, log = TRUE) {
          res = f1(theta1, theta2, log = log, ...)
          if(log) { res } else { exp(res) }
        }, init = w$f1.init, log = TRUE, level = w$f1.level, link = w$f1.link,
        method = w$f1.control$method)
        wts[i] = wts[i] + C1[i]
      }

      # store back-transformed integration nodes in grid
      nodes.backtransformed[i,] = theta2
    }

    # package results
    list(mix = mix, wts = wts, wts.e = wts.e, C1 = C1,
         nodes.backtransformed = nodes.backtransformed
    )
  }
  
  # unwrap results
  mix = pc$mix
  wts = pc$wts
  wts.e = pc$wts.e
  C1 = pc$C1
  grid$nodes = pc$nodes.backtransformed
  
  # normalize weights
  wts = exp(wts - mean(wts)) * grid$weights
  wts = wts / sum(wts)

  # normalize base weights for computing expectations
  wts.e = exp(wts.e - mean(wts.e)) * grid$weights
  wts.e = wts.e / sum(wts.e)

  # use C1(theta2) to build a second-layer function for dmix
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
  res = list(
    f = function(theta1, log = FALSE, ...) {
      dmix(x = theta1, f = h, params = mix, wts = wts, log = log, ...)
    },
    mix = mix,
    wts = wts,
    expectation = list(
      Eh = function(h, ...) { emix(h, grid$nodes, wts.e, ...) },
      Eh.precompute = function(h, ...) { emix(h, mix, wts.e, ...) },
      grid = grid,
      wts = wts.e
    )
  )
  class(res) = 'wtdMix'
  res
}
