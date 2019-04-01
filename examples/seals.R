# Use BISQuE to approximate the marginal posterior distribution for unknown
# population f(N|c, r) for the fur seals capture-recapture data example in 
# Givens and Hoeting (2013), example 7.10.

data('furseals')

# define theta transformation and jacobian
tx.theta = function(theta) { 
  c(log(theta[1]/theta[2]), log(sum(theta[1:2]))) 
}
itx.theta = function(u) { 
  c(exp(sum(u[1:2])), exp(u[2])) / (1 + exp(u[1])) 
}
lJ.tx.theta = function(u) {
  log(exp(u[1] + 2*u[2]) + exp(2*sum(u[1:2]))) - 3 * log(1 + exp(u[1]))
}

# compute constants
r = sum(furseals$m)
nC = nrow(furseals)

# set basic initialization for parameters
init = list(U = c(-.7, 5.5))
init = c(init, list(
  alpha = rep(.5, nC),
  theta = itx.theta(init$U),
  N = r + 1
))


post.alpha_theta = function(theta2, log = TRUE, ...) {
  # Function proportional to f(alpha, U1, U2 | c, r) 
  
  alpha = theta2[1:nC]
  u = theta2[-(1:nC)]
  theta = itx.theta(u)
  p = 1 - prod(1-alpha)
  
  res = - sum(theta)/1e3 - r * log(p) + lJ.tx.theta(u) - 
    nC * lbeta(theta[1], theta[2])
  for(i in 1:nC) {
    res = res + (theta[1] + furseals$c[i] - 1)*log(alpha[i]) + 
      (theta[2] + r - furseals$c[i] - 1)*log(1-alpha[i])
  }
  
  if(log) { res } else { exp(res) }
}

post.N.mixtures = function(N, params, log = TRUE, ...) {
  # The mixture component of the weighted mixtures for f(N | c, r)
  dnbinom(x = N-r, size = r, prob = params, log = log)
}

mixparams.N = function(theta2, ...) {
  # compute parameters for post.N.mixtures
  1 - prod(1 - theta2[1:nC])
}


w.N = wBuild(f = post.alpha_theta, init = c(init$alpha, init$U), 
             approx = 'gauss', link = c(rep('logit', nC), rep('identity', 2)))

m.N = wMix(f1 = post.N.mixtures, f1.precompute = mixparams.N, 
           f2 = post.alpha_theta, w = w.N)



# compute posterior mean
m.N$expectation$Eh.precompute(h = function(p) ((1-p)*r/p + r), 
                                   quadError = TRUE)

# compute posterior density
post.N.dens = data.frame(N = r:105)
post.N.dens$d = m.N$f(post.N.dens$N)

# plot posterior density
plot(d~N, post.N.dens, ylab = expression(f(N~'|'~bold(c),r)))

