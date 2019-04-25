library(fields)

simulate.field = function(n = 100, range = .3, smoothness = .5, phi = 1){
  # Simulates a mean-zero spatial field on the unit square
  #
  # Parameters:
  #  n - number of spatial locations
  #  range, smoothness, phi - parameters for Matern covariance function
  
  coords = matrix(runif(2*n), ncol=2)
  
  Sigma = Matern(d = as.matrix(dist(coords)), 
                 range = range, smoothness = smoothness, phi = phi)
  
  list(coords = coords,
       params = list(n=n, range=range, smoothness=smoothness, phi=phi),
       x = t(chol(Sigma)) %*%  rnorm(n))
}

# simulate data
x = simulate.field()

# configure gibbs sampler  
it = 100

# run sampler using default posteriors
post.samples = sFit(x = x$x, coords = x$coords, nSamples = it)

# build kriging grid
cseq = seq(0, 1, length.out = 10)
coords.krig = expand.grid(x = cseq, y = cseq)

# sample from posterior predictive distribution
burn = 75
samples.krig = sKrig(x$x, post.samples, coords.krig = coords.krig, burn = burn)
