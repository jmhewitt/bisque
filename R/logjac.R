#' Wrapper to abstractly evaluate log-Jacobian functions for transforms
#' 
#' @param x values at which to evaluate \eqn{J(x)}
#' @param link Character vector specifying link function for which the 
#'   inverse link function should be evaluated.  Supports \code{'identity'},
#'   \code{'log'}, and \code{'logit'}.
#' 
#' @examples 
#' smolBayes:::logjac(1, 'logit')
#' 
#' @seealso \code{\link{jac.log}}, \code{\link{jac.logit}}
logjac = function(x, link) {
  for(i in 1:length(link)) {
    x[i] = switch (link[i],
                   'identity' = 0,
                   'log' = jac.log(x[i], log = TRUE),
                   'logit' = jac.logit(x[i], log = TRUE)
    )
  }
  x
}