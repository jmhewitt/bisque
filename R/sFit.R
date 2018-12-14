#' Fit a spatially varying coefficient model
#'
#' @import Rcpp
#' 
#' @export
#'
#' @useDynLib smolBayes, .registration = TRUE
#' 
#' @example examples/spatial.R
#' 

sFit = function(x, coords, init, nSamples, thin=1, rw.initsd=.1, inits = list(),
                C=1, alpha=.44, priors = list(sigmasq = list(a=2, b=1), 
                rho = list(L=0, U=1), nu = list(L=0, U=1))) {
  
  d = as.matrix(dist(coords))
    
  if(is.null(inits)) {
    inits = list()
  }
  
  res <- .Call(`t_sfit`, as.double(x), as.matrix(d), as.double(priors$sigmasq$a),
            as.double(priors$sigmasq$b), as.double(priors$rho$L),
            as.double(priors$rho$U), as.double(priors$nu$L),
            as.double(priors$nu$U), as.integer(nSamples), as.integer(thin),
            as.double(rw.initsd), as.list(inits), as.double(C), 
            as.double(alpha))
  
  
  reslist = list(
    parameters = list(samples = res),
    priors = priors,
    coords = coords
  )
  
  class(reslist) = 'sFit'
  
  reslist
}