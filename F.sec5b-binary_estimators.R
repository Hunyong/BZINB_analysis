### Binary probability and their SE estimators
# Pr(Y1 = 0), Pr(Y1 = Y2 = 0), ...

# unlist.theta() function defined in F002-overall_rho.R
# unlist.theta = function(theta) {
#   theta = as.list(theta)
#   names(theta) = c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3")
#   theta[["p4"]] = 1 - (theta[["p1"]] + theta[["p2"]] + theta[["p3"]])
#   list2env(theta, envir = parent.frame())
#   NULL
# }


J1 = function(theta) {
  unlist.theta(theta) # Now a0, a1, a2, ... are individually available.
  p1 / (b1 + 1)^a1 / (b2 + 1)^a2 / (b1 + b2 + 1)^a0
}
J2 = function(theta) {
  unlist.theta(theta);   p2 / (b1 + 1)^(a0 + a1)
}
J3 = function(theta) {
  unlist.theta(theta);   p3 / (b2 + 1)^(a0 + a2)
}
J4 = function(theta) {
  unlist.theta(theta);   p4
}

K1 = function(theta) {
  unlist.theta(theta);   (p1 + p2) / (b1 + 1)^(a0 + a1)
}
K2 = function(theta) {
  unlist.theta(theta);   (p3 + p4)
}
L1 = function(theta) {
  unlist.theta(theta);   (p1 + p3) / (b2 + 1)^(a0 + a2)
}
L2 = function(theta) {
  unlist.theta(theta);   (p2 + p4)
}
f00 = function(theta) J1(theta) + J2(theta) + J3(theta) + J4(theta)  # Pr(Y1 = Y2 = 0)
f10 = function(theta) K1(theta) + K2(theta)      # Pr(Y1 = 0)
f20 = function(theta) L1(theta) + L2(theta)      # Pr(Y2 = 0) 
f0n = function(theta) - f00(theta) + f10(theta)  # Pr(Y1 = 0, Y2 > 0)
fn0 = function(theta) - f00(theta) + f20(theta)  # Pr(Y2 = 0, Y1 > 0)
fnn = function(theta) 1 + f00(theta) - f10(theta) - f20(theta)  # Pr(Y2 > 0, Y1 > 0)

delA = function(theta) {
  unlist.theta(theta)
  j1 = J1(theta)
  j2 = J2(theta)
  j3 = J3(theta)
  
  -c(a0 = (j1 + j2 + j3) * log(b1 + b2 + 1),
     a1 = (j1 + j2) * log(b1 + 1),
     a2 = (j1 + j3) * log(b2 + 1),
     b1 = j1 * (a0 / (b1 + b2 + 1) + a1 / (b1 + 1)) + j2 * (a0 + a1)/(b1 + 1),
     b2 = j1 * (a0 / (b1 + b2 + 1) + a1 / (b1 + 1)) + j3 * (a0 + a2)/(b2 + 1),
     p1 = 1 - j1/p1,
     p2 = 1 - j2/p2,
     p3 = 1 - j3/p3)
}
delB = function(theta) {
  unlist.theta(theta)
  k1 = K1(theta)
  k2 = K2(theta)
  
  -c(a0 = k1 * log(b1 + 1),
     a1 = k1 * log(b1 + 1),
     a2 = 0,
     b1 = k1 * (a0 + a1)/(b1 + 1),
     b2 = 0,
     p1 = 1 - k1/p1,
     p2 = 1 - k1/p1,
     p3 = 0)
}
delC = function(theta) {
  unlist.theta(theta)
  l1 = L1(theta)
  l2 = L2(theta)
  
  -c(a0 = l1 * log(b2 + 1),
     a1 = 0,
     a2 = l1 * log(b2 + 1),
     b1 = 0,
     b2 = l1 * (a0 + a2)/(b2 + 1),
     p1 = 1 - l1/p1,
     p2 = 0,
     p3 = 1 - l1/p1)
}
delD = function(theta) - delA(theta) + delB(theta) 
delE = function(theta) - delA(theta) + delC(theta) 
delF = function(theta) - delA(theta) + delB(theta) + delC(theta)

s00 = function(theta, covmat) sqrt(t(delA(theta)) %*% covmat %*% delA(theta))  # Pr(Y1 = Y2 = 0)
s10 = function(theta, covmat) sqrt(t(delB(theta)) %*% covmat %*% delB(theta))  # Pr(Y1 = 0)
s20 = function(theta, covmat) sqrt(t(delC(theta)) %*% covmat %*% delC(theta))  # Pr(Y2 = 0) 
s0n = function(theta, covmat) sqrt(t(delD(theta)) %*% covmat %*% delD(theta))  # Pr(Y1 = 0, Y2 > 0)
sn0 = function(theta, covmat) sqrt(t(delE(theta)) %*% covmat %*% delE(theta))  # Pr(Y2 = 0, Y1 > 0)
snn = function(theta, covmat) sqrt(t(delF(theta)) %*% covmat %*% delF(theta))  # Pr(Y2 > 0, Y1 > 0)

summary.bin = function(theta, covmat) {
  
  se.tab = cont.tab = matrix(NA, 3, 3)
  cont.tab[1,1] = f00(theta)
  cont.tab[1,2] = f0n(theta) # cont.tab[1,3] - cont.tab[1,1]
  cont.tab[2,1] = fn0(theta) # cont.tab[3,1] - cont.tab[1,1]
  cont.tab[2,2] = fnn(theta) # 1 - sum(cont.tab[1:2, 1:2], na.rm = TRUE)
  cont.tab[1:2,3] = apply(cont.tab[1:2,1:2], 1, sum)
  cont.tab[3,] = apply(cont.tab[1:2,], 2, sum)
  dimnames(cont.tab) = list(c("x=0", "x>0", "total"), c("y=0", "y>0", "total"))
  
  se.tab[1,1] = s00(theta, covmat)
  se.tab[1,2] = s0n(theta, covmat) 
  se.tab[1,3] = s10(theta, covmat) 
  se.tab[2,3] = se.tab[1,3]
  se.tab[2,1] = sn0(theta, covmat) 
  se.tab[2,2] = snn(theta, covmat)
  se.tab[3,1] = s20(theta, covmat) 
  se.tab[3,2] = se.tab[3,1]
  se.tab[3,3] = 0
  dimnames(se.tab) = list(c("x=0", "x>0", "total"), c("y=0", "y>0", "total"))
  
  list(est = cont.tab, se = se.tab)
}


### This code is wrong. Please see summary.bin() in F0201-binary_estimators.R
# # 2 x 2 contingency table
# cont.table = function(coef) {
#   cont.tab = matrix(NA, 3, 3)
#   cont.tab[1,1] = bzinb:::lik.bzinb(0,0, param = coef) %>% exp  #0,0: 0.122
#   cont.tab[1,2] = dnbinom(0, size = coef[c(1,2)] %>% sum, prob = 1/(coef["b1"] + 1)) ## Something here is wrong.
#   cont.tab[2,1] = dnbinom(0, size = coef[c(1,3)] %>% sum, prob = 1/(coef["b2"] + 1)) ## Something here is wrong.
#   cont.tab[1,2] = cont.tab[1,2] - cont.tab[1,1]
#   cont.tab[2,1] = cont.tab[2,1] - cont.tab[1,1]
#   cont.tab[2,2] = 1 - sum(cont.tab[1:2, 1:2], na.rm = TRUE)
#   cont.tab[1:2,3] = apply(cont.tab[1:2,1:2], 1, sum)
#   cont.tab[3,] = apply(cont.tab[1:2,], 2, sum)
#   dimnames(cont.tab) = list(c("x=0", "x>0", "total"), c("y=0", "y>0", "total"))
#   cont.tab
# }
