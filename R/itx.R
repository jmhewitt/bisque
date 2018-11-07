#' Named inverse transformation functions
#' 
#' Evaluates the inverse of the named link function at the locations
#' \code{x}.
#' 
#' @param x Values at which to evaluate the inverse link function
#' @param link Character vector specifying link function for which the 
#'   inverse link function should be evaluated.  Supports \code{'identity'},
#'   \code{'log'}, and \code{'logit'}.
#'   
#' @examples 
#' itx(0, 'logit')
#' 
itx = function(x, link) {
  # function to invert link transformations
  for(i in 1:length(link)) {
    x[i] = switch (link[i],
                   'identity' = x[i],
                   'log' = exp(x[i]),
                   'logit' = plogis(x[i])
    )
  }
  x
}