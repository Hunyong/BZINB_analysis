# density plots
densitymap <- function(xvec, yvec) {
  xy.reduced <- as.data.frame(table(xvec,yvec))
  names(xy.reduced) <- c("x", "y","freq")
  xy.reduced <- xy.reduced[xy.reduced$freq != 0,]
  xy.reduced$x <- as.numeric(as.character(xy.reduced$x))
  xy.reduced$y <- as.numeric(as.character(xy.reduced$y))
  xy.reduced$density <- as.numeric(as.character(xy.reduced$freq))/length(xvec)
  xy.reduced
}

densitymap.parm <- function(densityFUN, ...) {
  obj <- expand.grid(0:50, 0:50)
  names(obj) = c("x", "y")
  obj$density = densityFUN(obj$x, obj$y, ...)
  obj
}

# measures
entropy.generic <- function(..., densityfn, summation = 100) {
  if (sum(is.na(c(...)))>0) {return(data.frame(H1 = NA, H2 = NA, H12 = NA, MI = NA))} else {
    joint <- sapply(0:summation, function(x) {sapply(0:summation, function(y) {
      densityfn (x, y, ...)})})
    marginal1 <- apply(joint,1,sum)
    marginal2 <- apply(joint,2,sum)
    
    H1 <- - marginal1 * log(marginal1); H1[!is.finite(H1)] <- 0; H1 <- sum(H1)
    H2 <- - marginal2 * log(marginal2); H2[!is.finite(H2)] <- 0; H2 <- sum(H2)
    H12 <- - joint * log(joint);        H12[!is.finite(H12)] <- 0; H12 <- sum(H12)
    M12 <- max(- H12 + H1 + H2,0)
    return(data.frame(H1 = H1, H2 = H2, H12 = H12, MI = M12)) 
  }
} 

dbnb = function(...) exp(bzinb::lik.bnb(...))

measures <- function(xvec, yvec) {
  require(bzinb)
  BZINB = bzinb(xvec, yvec, maxiter = 50000)
  abp = BZINB$coefficients[, 1]
  list(measure =
         data.frame (
           PC = cor(xvec, yvec),
           SC = cor(xvec, yvec, method = "spearman"),
           tau = cor(xvec, yvec, method = "kendall"),
           MI = entropy::mi.empirical(table(xvec,yvec)),
           BZINB.cor = BZINB$rho[1, 1],
           BZINB.MI = entropy.generic(abp[1], abp[2], abp[3], abp[4], abp[5], 
                                      densityfn=dbnb)["MI"]
         ),
       parameters = BZINB$coefficients,
       rho = BZINB$rho, lik = BZINB$lik, iter = BZINB$iter,
       rho.ci = ci(rho[1, 1], rho[1, 2], logit = TRUE, alpha = 0.05)) %>% print
}

RHO = function(a0, a1, a2, b1, b2) {(a0^2 / (a0 + a1) / (a0 + a2) * b1 / (b1 + 1) * b2 / (b2 + 1))^.5}


ci <- function(rho, se.rho, logit = FALSE, alpha = 0.05, ...) {
  CI <- rho + c(lb = -1, ub = 1) * qnorm(1 - alpha/2) * se.rho
  if (logit) CI <- plogis(CI)
  
  return(CI)
}

longform <- function(dat, name) {
  require(tidyverse)
  mutate(dat, subgroup = factor(rep(1:2, 20), labels = c("-a", "-b"))) %>% # subgroup for shape
    mutate(zi = factor(zi, levels = c("low", "mid-B", "mid-U", "high-B", "high-U"))) %>% 
    gather(key = n, value = !!name, -id, -rho, -zi, -subgroup) %>% 
    mutate(n = factor(gsub("n\\.", "", n), levels = sampsize), rho = factor(round(rho,2), levels = c(.6, .3, .1, .01)))
}