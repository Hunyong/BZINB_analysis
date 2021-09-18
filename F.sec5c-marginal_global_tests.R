
# marginal mean estimators
mu1 = function(theta) {
  unlist.theta(theta) # Now a0, a1, a2, ... are individually available.
  (p1 + p2) * (a0 + a1) * b1
}


mu2 = function(theta) {
  unlist.theta(theta) # Now a0, a1, a2, ... are individually available.
  (p1 + p3) * (a0 + a2) * b2
}


# marginal mean estimators

del.mu.1 = function(theta) {
  unlist.theta(theta) # Now a0, a1, a2, ... are individually available.
  p12 = p1 + p2
  a01 = a0 + a1
  c(a0 = p12 * b1, a1 = p12 * b1, a2 = 0, b1 = p12 * a01, b2 = 0,
    p1 = a01 * b1, p2 = a01 * b1, p3 = 0)
}

del.mu.2 = function(theta) {
  unlist.theta(theta) # Now a0, a1, a2, ... are individually available.
  p13 = p1 + p3
  a02 = a0 + a2
  c(a0 = p13 * b2, a1 = 0, a2 = p13 * b2, b1 = 0, b2 = p13 * a02,
    p1 = a02 * b2, p2 = 0, p3 = a02 * b2)
}

est.mu = function(theta, covmat) {
  mu1 = mu1(theta)
  mu2 = mu2(theta)
  d1 = del.mu.1(theta)
  d2 = del.mu.2(theta)
  d  = cbind(d1, d2)
  
  list(mu = c(mu1 = mu1, mu2 = mu2), 
       cov = t(d) %*% covmat %*% d)
}

zinb.mu = function(vec) {
  a = pscl::zeroinfl(vec ~ 1, dist = "negbin") # ....
  mu.nz  = a$coefficients$count %>% exp %>% as.numeric
  # sig = a$theta %>% as.numeric
  pi  = a$coefficients$zero %>% plogis %>% as.numeric
  mu = (1-pi) * mu.nz
  del.g = mu * c(1, pi)
  covar = t(del.g) %*% a$vcov %*% del.g
  
  c(mu = mu, se = sqrt(covar))
}
