#' Fit a spatially varying coefficient model
#'
#' @import Rcpp
#' @import foreach
#' @import doRNG
#' @importFrom itertools ichunk
#'
#' @export
#'
#' @useDynLib smolBayes, .registration = TRUE
#' 
#' @example examples/spatial.R
#' 

sKrig = function(x, sFit, coords.krig, coords = sFit$coords, burn = 0,
                 ncores = 1) {
  
  if(burn > 0) {
    sFit$parameters$samples$sigmasq = sFit$parameters$samples$sigmasq[-(1:burn)]
    sFit$parameters$samples$rho = sFit$parameters$samples$rho[-(1:burn)]
    sFit$parameters$samples$nu = sFit$parameters$samples$nu[-(1:burn)]
  }
  
  d = as.matrix(dist(rbind(as.matrix(coords.krig), as.matrix(coords))))
  
  n0 = nrow(coords.krig)
  
  d00 = d[1:n0, 1:n0]
  d01 = d[1:n0, -(1:n0)]
  d11 = d[-(1:n0), -(1:n0)]
  
  # make looping more efficient
  mcoptions = list(preschedule=FALSE)
  
  # estimate chunksize that will minimize number of function calls
  nSamples = length(sFit$parameters$samples$sigmasq)
  chunkSize = ceiling((nSamples+1)/ncores)
  
  op = ifelse(ncores>1, `%dorng%`, `%do%`)
  
  res = op(foreach(inds = ichunk(1:nSamples, chunkSize = chunkSize, 
                              mode = 'numeric'), 
                .combine = 'rbind', .options.multicore = mcoptions, 
                .packages = 'Rcpp'), {
    
    .Call(`t_spredict`, as.numeric(x), as.matrix(d00), as.matrix(d01),
          as.matrix(d11), 
          as.numeric(sFit$parameters$samples$sigmasq[inds]),
          as.numeric(sFit$parameters$samples$rho[inds]),
          as.numeric(sFit$parameters$samples$nu[inds]))$x0
  })
  
  reslist = list(
    samples = res,
    coords.krig = coords.krig
  )
  
  class(reslist) = 'sKrig'
  
  reslist
}
