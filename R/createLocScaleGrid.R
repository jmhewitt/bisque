#' Create a centered and scaled sparse integration grid
#' 
#' Enhances mvQuad::createNIGrid by shifting and scaling a sparse integration
#' grid, and evaluating the weight function at each of the grid nodes.
#' 
#' @export
#' 
#' @importFrom mvQuad createNIGrid
#' 
#' @param mu location at which grid should be centered
#' @param prec "precision matrix" associated with the integration grid.  When 
#'   building a sparse integration grid for a density, \code{prec} is often 
#'   the negative of the hessian at the mode.
#' @param level accuracy level (typically number of grid points for the 
#'   underlying 1D quadrature rule) [description from mvQuad::createNIGrid]
#' 
createLocScaleGrid = function(mu = 0, prec = 1, level = 2) {
  # determine standardized quadrature points
  grid = createNIGrid(dim = length(mu), level = level, type = 'nHN', 
                      ndConstruction = 'sparse')
  # evaluate the weight function at each quadrature point
  grid$d = apply(dnorm(grid$nodes, log = TRUE), 1, sum)
  # center and scale integration grid around posterior mode
  prec.chol = chol(prec)
  grid$nodes = sweep(t(solve(prec.chol, t(grid$nodes))), 2, - mu)
  # add jacobian
  grid$d = grid$d + sum(log(diag(prec.chol)))

  grid
}