#' Named transformation functions
#' 
#' Evaluates the named link function at the locations \code{x}.
#' 
#' @param x Values at which to evaluate the link function
#' @param link Character vector specifying link function to evaluate.  Supports 
#'   \code{'identity'}, \code{'log'}, and \code{'logit'}.
#'   
#' @examples 
#' tx(0.5, 'logit')
#' 
tx = function(x, link) {
  # function to apply link transformations
  for(i in 1:length(link)) {
    x[i] = switch (link[i],
                   'identity' = x[i],
                   'log' = log(x[i]),
                   'logit' = qlogis(x[i])
    )
  }
  x
}