zi.test = 
  function(vec) {
    vec = as.integer(vec)
    zp  = mean(vec == 0, na.rm = TRUE)
    if (zp > 1 - 1e-7) return(c(chibsq = NA, p = NA, LL_ZINB = NA, LL_NB = NA, zero.prop = zp))
    llnb   = MASS::glm.nb(vec ~ 1)$twologlik / 2
    llzinb = pscl::zeroinfl(vec ~ 1, dist = "negbin")$loglik
    statistic = 2 * (llzinb - llnb)
    p = pchisq(q = statistic, df = 1, lower.tail = FALSE) / 2 # 50:50 chisq 0 and 1 distributions
    if (statistic <= 0) p = 1.0
    return(c(chibsq = statistic, p = p, LL_ZINB = llzinb, LL_NB = llnb, zero.prop = zp))
  }
zi.test.group = 
  function(vec, Tx) {
    vec = as.integer(vec)
    zp  = mean(vec == 0, na.rm = TRUE)
    if (zp > 1 - 1e-7) return(c(chibsq = NA, p = NA, LL_ZINB = NA, LL_NB = NA, zero.prop = zp))
    llnb   = MASS::glm.nb(vec ~ Tx)$twologlik / 2
    llzinb = pscl::zeroinfl(vec ~ Tx, dist = "negbin")$loglik
    statistic = 2 * (llzinb - llnb)
    p = pchisq(q = statistic, df = 1, lower.tail = FALSE) / 2 # 50:50 chisq 0 and 1 distributions
    if (statistic <= 0) p = 1.0
    return(c(chibsq = statistic, p = p, LL_ZINB = llzinb, LL_NB = llnb, zero.prop = zp))
  }

# for (i in 1:10)
#   zi.test(data[i,-1]) %>% print

if (0) {
  zi.result = sapply(1:1000, function(i) zi.test(data[i,-1])) %>% t
  hist(zi.result[,"p"]); mean(zi.result[, "p"] < 0.05, na.rm = TRUE)
}



# ZI plots
param.zinb2 <- function(zinb.obj) {
  #if (class(zinb.obj)[1] != "zeroinfl") stop("Designed only for zeroinfl object.")
  if (class(zinb.obj)[1] == "zeroinfl") {
    mu  = zinb.obj$coefficients$count %>% exp %>% as.numeric
    sig = zinb.obj$theta %>% as.numeric
    pi  = zinb.obj$coefficients$zero %>% plogis %>% as.numeric
  } else {
    mu  = zinb.obj$coefficients %>% exp %>% as.numeric
    sig = zinb.obj$theta %>% as.numeric
    pi = NA
  }
  return(c(mu = mu, sig = sig, pi = pi))
}
dnb <- function(par, rng = 0:10) {
  mu    = par["mu"]
  theta = par["sig"]  # =alpha
  beta  = mu/theta
  phi   = 1/(beta + 1)
  
  dnbinom(rng, size = theta, prob = phi)
}
plot.nb <- function (y, par = NULL, rng = 0: max(y)) {
  if (is.null(par)) par = param.zinb2(MASS::glm.nb(y ~ 1))
  n = length(y)
  data = data.frame(y = rng, count = dnb(par, rng = rng) * n)
  ggplot(data = data.frame(y = y), mapping = aes(y)) + 
    stat_count() +
    # geom_line(data = data, aes(y, count), col = "red") +
    geom_point(data = data, aes(y, count), col = "red", size = 2) +
    theme_bw() +
    ggtitle("NB empirical (hisogram) vs NB estimates (points)")
}
plot.zinb <- function (y = NULL, par = NULL, rng = 0:max(y)) {
  if (is.null(par)) par = param.zinb2(pscl::zeroinfl(y ~ 1, dist = "negbin"))
  n = length(y)
  data = data.frame(y = rng, count = dnb(par, rng = rng) * n * (1 - par["pi"]))
  data[data$y == 0, "count"] = data[data$y == 0, "count"] + n * par["pi"]  # zero inflation
  ggplot(data = data.frame(y = y), mapping = aes(y)) + 
    stat_count() +
    # geom_line(data = data, aes(y, count), col = "red") +
    geom_point(data = data, aes(y, count), col = "red", size = 2)  +
    theme_bw() +
    ggtitle("NB empirical (hisogram) vs NB estimates (points)")
}
plot.2nb <- function (y, xlab = "") { # combined plot (ZINB + NB)
  rng = 0:max(y)
  n = length(y)
  par.nb = param.zinb2(MASS::glm.nb(y ~ 1))
  par.zinb = param.zinb2(pscl::zeroinfl(y ~ 1, dist = "negbin"))
  
  # ZINB model
  data = data.frame(y = rng, 
                    count = dnb(par.zinb, rng = rng) * n * (1 - par.zinb["pi"]),
                    model = "ZINB")
  data[data$y == 0 & data$model == "ZINB", "count"] = 
    data[data$y == 0 & data$model == "ZINB", "count"] + 
    n * par.zinb["pi"]  # zero inflation
  
  # NB model
  data.nb = data.frame(y = rng, count = dnb(par.nb, rng = rng) * n,
                       model = "NB")
  data =  rbind(data, data.nb) # %>% print
  
  ggplot(data = data.frame(y = rep(y, times = 2), model = rep(c("ZINB", "NB"), each = n)), 
         mapping = aes(y)) + 
    facet_grid(. ~ model) +
    stat_count(fill = "lemonchiffon2") +
    # geom_line(data = data, aes(y, count), col = "red") +
    geom_point(data = data, aes(y, count, col = model), size = 2)  +
    guides(col = FALSE) +
    xlab(xlab) +
    theme_bw() # +
    #ggtitle("NB empirical (hisogram) vs NB estimates (points)")
}