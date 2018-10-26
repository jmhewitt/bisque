itx = function(x, link) {
  # function to invert link transformations
  for(i in 1:length(link)) {
    x[i] = switch (link[i],
                   'identity' = x,
                   'log' = exp(x[i]),
                   'logit' = plogis(x[i])
    )
  }
  x
}