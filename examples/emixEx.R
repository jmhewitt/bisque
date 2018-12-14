# density will be a mixture of betas
params = matrix(exp(2*runif(10)), ncol=2)

# mixture components are equally weighted
wts = rep(1/nrow(params), nrow(params))

# compute mean of distribution by cycling over each mixture component
h = function(p) { p[1] / sum(p) }

# compute mixture mean
mean.mix = emix(h, params, wts)

# (comparison) Monte Carlo estimate of mixture mean
nsamples = 1e4
component = sample(x = 1:length(wts), size = nsamples, prob = wts, 
                   replace = TRUE)
x = sapply(component, function(cmp) {
  rbeta(n = 1, shape1 = params[cmp, 1], shape2 = params[cmp, 2])
})
mean.mix.mc = mean(x)

# compare estimates
c(emix = mean.mix, MC = mean.mix.mc)
