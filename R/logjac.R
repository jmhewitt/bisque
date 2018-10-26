logjac = function(x, link) {
  # function to apply jacobians for link transformations
  for(i in 1:length(link)) {
    x[i] = switch (link[i],
                   'identity' = 0,
                   'log' = jac.log(x[i], log = TRUE),
                   'logit' = jac.logit(x[i], log = TRUE)
    )
  }
  x
}