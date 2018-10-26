tx = function(x, link) {
  # function to apply link transformations
  for(i in 1:length(link)) {
    x[i] = switch (link[i],
                   'identity' = x,
                   'log' = log(x[i]),
                   'logit' = qlogis(x[i])
    )
  }
  x
}