# evaluate mixture density at these locations
x = seq(0, 1, length.out = 100)

# density will be a mixture of beta distributions
f = function(x, theta, log = FALSE) {
  dbeta(x, shape1 = theta[1], shape2 = theta[2], log = log)
}

# beta parameters are randomly assigned
params = matrix(exp(2*runif(10)), ncol=2)

# mixture components are equally weighted
wts = rep(1/nrow(params), nrow(params))

# evaluate mixture density
fmix = dmix(x = x, f = f, params = params, wts = wts)

# plot mixture density
plot(x, fmix, type='l', ylab = expression(f(x)), 
     ylim = c(0, 4))

# plot component densities
for(i in 1:length(wts)){
  curve(f(x, params[i,]), col = 2, add = TRUE)
}
