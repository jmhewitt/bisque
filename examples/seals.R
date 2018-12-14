# load data
data('furseals')
dat = furseals
rm(furseals)

# compute constants
r = sum(dat$m)
nC = nrow(dat)

# set basic initialization for parameters
init = list(alpha=rep(.5, nC), theta=rep(.5,2), N = r+1)


#
# weighted mixtures posterior for N
#

post.alpha_theta = function(theta2, log = TRUE, ...) {
  # Function proportional to f(alpha, theta1, theta2 | c, r) 
  
  alpha = theta2[1:nC]
  theta = theta2[-(1:nC)]
  
  res = - sum(theta)/1000 - r * log(1-prod(1-alpha))
  for(i in 1:nC) {
    res = res + dbeta(alpha[i], theta[1] + dat$c[i], theta[2] + r - dat$c[i], 
                      log = TRUE) - lbeta(theta[1], theta[2]) + 
      lbeta(theta[1] + dat$c[i], theta[2] + r - dat$c[i])
  }
  
  if(log) { res } else { exp(res) }
}


post.N.mixtures = function(N, params, log = TRUE) {
  # The mixture component of the weighted mixtures for f(N | c, r)
  dnbinom(x = N-r, size = r, prob = params, log = log)
}

mixparams.N = function(theta2) {
  # compute parameters for post.N.mixtures
  1 - prod(1 - theta2[1:nC])
}

# build a weighted mixture of posteriors
post.N = wtdMix(
  f1 = post.N.mixtures,
  f1.precompute = mixparams.N,
  f2 = post.alpha_theta,
  w.init = c(init$alpha, init$theta),
  w.link = c(rep('logit', nC), rep('log', 2))
)

# compute posterior mean
post.N$expectation$Eh.precompute(h = function(p) {((1-p)*r/p + r)}, 
                                 quadError = TRUE)
