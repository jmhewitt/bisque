#' Construct a weighted mixture object
#'
#'
#' @export
#'
#'
#'
wMix = function(f1, f2, w, f1.precompute = function(x, ...){x}, spec = 'ff', 
                level = 2, c.int = NULL, c.level = 2, c.init = NULL,
                c.link = rep('identity', length(c.init)),
                c.link.params = rep(list(NA), length(c.init)),
                c.optim.control = list(maxit = 5e3, method = 'BFGS'),
                ncores = 1, quadError = TRUE, ...) {
  
  # determine if intermediate integration constants need to be computed
  spec.split = strsplit(spec, character(0))[[1]]
  f1.cst = ifelse(spec.split[1] == 'g', TRUE, FALSE)
  f2.cst = ifelse(spec.split[2] == 'g', TRUE, FALSE)
  f.cst = any(f1.cst, f2.cst)
  if(f.cst) {
    if(any(is.null(c.int), is.null(c.init))) {
      stop('Must provide c.int AND c.init.')
    }
  }

  # create integration grid
  grid = w$gridfn(level = level, quadError = quadError)

  op = ifelse(ncores > 1, `%dopar%`, `%do%`)

  # precompute parameters and weights for posterior mixture for theta1
  p0 = f1.precompute(itx(grid$nodes[1,], w$link, w$link.params), ...)
  nodes  = nrow(grid$nodes)
  chunkSize = ceiling(nodes/ncores)
  pc = op(foreach(
    inds = ichunk(1:nodes, chunkSize = chunkSize, mode = 'numeric'),
    .combine = mergePars, .export = c('itx', 'logjac', 'kCompute')), {

    # initialize return objects
    mix = matrix(NA, nrow = length(inds), ncol = length(p0))
    wts = numeric(nrow(mix))
    wts.e = numeric(nrow(mix))
    C1 = numeric(nrow(mix))
    nodes.backtransformed = grid$nodes[inds, , drop = FALSE]

    for(i in 1:nrow(mix)) {

      # back-transform parameters
      theta2 = as.numeric(itx(grid$nodes[inds[i],], w$link, w$link.params))

      # compute mixture parameters
      mix[i,] = f1.precompute(theta2, ...)

      # compute base weights (i.e., the weight function ratios)
      wts[i] = f2(theta2, log = TRUE, ...) +
        sum(logjac(grid$nodes[inds[i],], w$link, w$link.params)) -
        grid$d[inds[i]]

      # compute weights for evaluating expectations in secondary analyses
      wts.e[i] = wts[i]

      # as necessary, compute C1(theta2) and update weights
      if(f.cst) {
        C1[i] = kCompute(f = function(theta1, log = TRUE) {
          res = c.int(theta1, theta2, log = log, ...)
          if(log) { res } else { exp(res) }
        }, init = c.init, log = TRUE, level = c.level, link = c.link,
        linkparams = c.link.params, method = c.optim.control$method,
        maxit = c.optim.control$maxit)
        
        if(f2.cst) { wts[i] = wts[i] + C1[i] }
      }

      # store back-transformed integration nodes in grid
      nodes.backtransformed[i,] = theta2
    }

    # package results
    list(mix = mix, wts = wts, wts.e = wts.e, C1 = C1,
         nodes.backtransformed = nodes.backtransformed
    )
  })

  # unwrap results
  mix = pc$mix
  wts = pc$wts
  wts.e = pc$wts.e
  C1 = pc$C1
  grid$nodes = pc$nodes.backtransformed

  # standardize quadError weights before the main quadrature weights
  if(quadError) {
    grid$errorNodes$weights = grid$errorNodes$weights *
      exp(wts[grid$errorNodes$inds]-mean(wts[grid$errorNodes$inds]))
    grid$errorNodes$weights =
      grid$errorNodes$weights / sum(grid$errorNodes$weights)
  }

  # normalize weights
  wts = exp(wts - mean(wts)) * grid$weights
  wts = wts / sum(wts)

  # use C1(theta2) to build a second-layer function for dmix
  if(!f1.cst) { h = f1 }
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
    f = function(theta1, log = FALSE, quadError = FALSE, ...) {
      if(quadError) {
        dmix(x = theta1, f = h, params = mix, wts = wts, log = log,
             errorNodesWts = grid$errorNodes, ...)
      } else {
        dmix(x = theta1, f = h, params = mix, wts = wts, log = log, ...)
      }
    },
    mix = mix,
    wts = wts,
    expectation = list(
      Eh = function(h, ncores = 1, quadError = FALSE, ...) {
        if(quadError) {
          emix(h = h, params = grid$nodes, wts = wts, ncores = ncores,
               errorNodesWts = grid$errorNodes, ...)
        } else {
          emix(h = h, params = grid$nodes, wts = wts, ncores = ncores, ...) }
      },
      Eh.precompute = function(h, ncores = 1, quadError = FALSE, ...) {
        if(quadError) {
          emix(h = h, params = mix, wts = wts, ncores = ncores,
               errorNodesWts = grid$errorNodes, ...)
        } else {
          emix(h = h, params = mix, wts = wts, ncores = ncores, ...) }
      },
      grid = grid
    ),
    ratios = wts.e
  )
  class(res) = 'wMix'
  res
}
