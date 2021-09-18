library(sas7bdat); library(dplyr); library(devtools); library(bzinb); library(ggplot2)
dat.dental <- read.sas7bdat("../Data-XACT-AGE-HENKE/SASDATA/le0220402.sas7bdat")
source("F.sec5a-zeroInflationTest.R")
source("F.sec5b-binary_estimators.R")
source("F.sec5c-marginal_global_tests.R")

### 0. data
  # smooth surface: front and back faces
  # occlusal surface: top wiggly surfaces
  # proximal surface: side, neighboring surfaces
  
  # xylitol may increase the association between smooth and occlusal
  # However, since proximal surface may not be affected much, the associaion might decrease.

### 1. Zero-inflation test for each of the three variables
  zi.result = data.frame()
  zi.result = zi.test(dat.dental$CUM_D2_PROXIMAL_THREE_YEARS) %>% t %>% as.data.frame
  zi.result[2, ] = zi.test(dat.dental$CUM_D2_OCCLUSAL_THREE_YEARS)
  zi.result[3, ] = zi.test(dat.dental$CUM_D2_SMOOTH_THREE_YEARS)
  rownames(zi.result) = c("CUM_D2_PROXIMAL_THREE_YEARS", "CUM_D2_OCCLUSAL_THREE_YEARS", "CUM_D2_SMOOTH_THREE_YEARS")
  print(zi.result)
  
  ### Tx-conditional test: The same conclusion at 5% significance level
  zi.result[1, ] = zi.test.group(dat.dental$CUM_D2_PROXIMAL_THREE_YEARS, dat.dental$treatment) %>% t %>% as.data.frame
  zi.result[2, ] = zi.test.group(dat.dental$CUM_D2_OCCLUSAL_THREE_YEARS, dat.dental$treatment)
  zi.result[3, ] = zi.test.group(dat.dental$CUM_D2_SMOOTH_THREE_YEARS, dat.dental$treatment)
  rownames(zi.result) = c("CUM_D2_PROXIMAL_THREE_YEARS", "CUM_D2_OCCLUSAL_THREE_YEARS", "CUM_D2_SMOOTH_THREE_YEARS")
  print(zi.result)
  
  # 1.1 proximal (side by side): ZINB
  plot.2nb(y = dat.dental$CUM_D2_PROXIMAL_THREE_YEARS, xlab = "number of proximal-surface caries") + 
    geom_text(data = data.frame(model = "NB"), parse = TRUE, 
              aes(x =10, y = 150, label = "bar(chi)^2 == 1.22~(p == '0.13')")) + # ZINB has a much better fit
    ylab("frequency")
  ggsave("../figure/xact/Figure0dentalZI1.png")
  
  # 1.2 occlusal (top): NB is enough
  plot.2nb(y = dat.dental$CUM_D2_OCCLUSAL_THREE_YEARS, xlab = "number of occlusal-surface caries")  + 
    geom_text(data = data.frame(model = "NB"), parse = TRUE, 
              aes(x = 10, y = 300, label = "bar(chi)^2 == 0.00~(p == '1.00')")) + # NB is good enough
    ylab("frequency")
  ggsave("../figure/xact/Figure0dentalZI2.png")
  
  # 1.3 smooth (front/back): ZINB
  plot.2nb(y = dat.dental$CUM_D2_SMOOTH_THREE_YEARS, xlab = "number of smooth-surface caries") + 
    geom_text(data = data.frame(model = "NB"), parse = TRUE, 
              aes(x = 10, y = 140, label = "bar(chi)^2 == 7.88~(p == '0.0025')")) + # ZINB has a much better fit
    ylab("frequency")
  ggsave("../figure/xact/Figure0dentalZI3.png")
  

### 2. Model fitting
  trt = dat.dental$treatment == 1
  
  # 2.2 PS: proximal & smooth   overall 0.597, control 0.599, case = 0.609
  x = dat.dental$CUM_D2_PROXIMAL_THREE_YEARS
  y = dat.dental$CUM_D2_SMOOTH_THREE_YEARS
  bnb.ps = bnb(x, y)
  bzinb.ps = bzinb(x, y, maxiter = 50000, showFlag = TRUE, vcov = TRUE)
  
  # Non-xylitol group rho = 0.599 (se = 0.054)
  bzinb.ps.xy0 = bzinb(x[!trt], y[!trt], maxiter = 100000, showFlag = TRUE, vcov = TRUE) # 67308 iterations
  # Xylitol group rho = 0.609 (se = 0.058)
  bzinb.ps.xy1 = bzinb(x[trt], y[trt], maxiter = 100000, showFlag = TRUE, vcov = TRUE)
  
  # binary contingency table - model-based
  cont.table(bzinb.ps.xy0$coefficients[,1])  # 12.19 / 5.97 / 3.83 / 78.01
  cont.table(bzinb.ps.xy1$coefficients[,1])  # 11.91 / 9.82 / 4.78 / 73.50  (9.81)
  
  res0 = summary.bin(bzinb.ps.xy0$coefficients[,1], bzinb.ps.xy0$vcov)
  res1 = summary.bin(bzinb.ps.xy1$coefficients[,1], bzinb.ps.xy1$vcov)
  
  # Unequal-variance z-test
  t.diff = (res0$est - res1$est) / sqrt(res0$se^2 + res1$se^2) %>% print
  p.diff = 2*(1 - pnorm(abs(t.diff))) %>% print
  
  # binary contingency table - empirical
  {table(ifelse(x[!trt], "x>0", "x=0") , ifelse(y[!trt], "y>0", "y=0")) /length(x[!trt])} %T>% 
    print %>% {sqrt(. * (1-.) / dim(dat.dental)[1])}
  table(ifelse(x[trt], "x>0", "x=0") , ifelse(y[trt], "y>0", "y=0")) /length(x[trt]) 
  
  
### marginal mean test  
  
  # 1. BZINB-based marginal test
  mu.control = est.mu(bzinb.ps.xy0$coefficients[,1], bzinb.ps.xy0$vcov)
  mu.interv  = est.mu(bzinb.ps.xy1$coefficients[,1], bzinb.ps.xy1$vcov)
  
  # diff in surface 1: 0.212, z = 0.831, p = 0.203 (two-sided), 0.101 (one-sided)
  (mu.control$mu["mu1"] - mu.interv$mu["mu1"])/sqrt(mu.interv$cov[1, 1] + mu.control$cov[1, 1])
  # diff in surface 2: 0.380, z = 1.266, p = 0.103 (two-sided), 0.051 (one-sided)
  (mu.control$mu["mu2"] - mu.interv$mu["mu2"])/sqrt(mu.interv$cov[2, 2] + mu.control$cov[2, 2])
  
  # 1B. BZINB-based global test
  t(mu.interv$mu - mu.control$mu) %*% solve(mu.interv$cov + mu.control$cov) %*% (mu.interv$mu - mu.control$mu)
  # chi2 = 1.97, 2 df, p = 0.373
  
  # 2. ZINB-based univariate test 
  # 2.1 diff in surface 1: 2.683 - 2.480 = 0.203, z = 0.628, p = 0.265, 0.132 (one-sided)
  mu.control.zinb = zinb.mu(x[trt]) %>% print
  mu.interv.zinb = zinb.mu(x[!trt]) %>% print
  (mu.interv.zinb["mu"] - mu.control.zinb["mu"])/sqrt(mu.interv.zinb["se"]^2 + mu.control.zinb["se"]^2)
  
  # 2.1 diff in surface 2: 3.311 - 2.947 = 0.364, z = 0.801, p = 0.212, 0.106 (one-sided)
  mu.control.zinb = zinb.mu(y[trt]) %>% print
  mu.interv.zinb = zinb.mu(y[!trt]) %>% print
  (mu.interv.zinb["mu"] - mu.control.zinb["mu"])/sqrt(mu.interv.zinb["se"]^2 + mu.control.zinb["se"]^2)