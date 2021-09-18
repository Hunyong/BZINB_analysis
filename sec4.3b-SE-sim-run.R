#### 0. library ####
  # install.packages("bzinb")
  library(bzinb)
  library(tidyverse); library(magrittr)
  source("sec4.3a-SE-sim-parm.R")
  
  args = commandArgs(trailingOnly=TRUE)  # passed from script
  print(args)
  i = as.numeric(args[1]) # i = 1..40: param. ex param = abp[i,]
  j = as.numeric(args[2]) # j = 1..5: n       ex n = sampsize[j]
  maxiter <- 50000L
  n.sim <- 1000L
  
  # output directory
  output_folder <- "../../output/C0301SE-MI50K/"
  if (!dir.exists(output_folder)) dir.create(output_folder)
  
  # picking up simulation scenario
  param <- abp[i, c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4")] %>% unlist
  rho <- abp[i, "rho"]
  zi <- abp[i, "zi"]
  n <- sampsize[j]
  abp.names <- c("a0", "a1", "a2", "b1", "b2", "p1", "p2", "p3", "p4") # global variable
  est.names <- list(1:n.sim, c(abp.names, "rho", "logit.rho", paste0("se.", abp.names), 
                               "se.rho", "se.logit.rho", 
                               "lb.95","ub.95","cov.95","lb.n.95","ub.n.95","cov.n.95",
                               "iter", "time"))
  fn <- paste0(output_folder, "result-", i, "-", j , ".sim")
  fn.tmp <- paste0(output_folder, "tmp.result-", i, "-", j , ".sim")
  # if final output is there, stop
  if (file.exists(fn)) {
    cat("done already.\n")
    stop("done already.")
  }
  
  # if the intermediate output is there, use it and start from the last row.
  if (file.exists(fn.tmp)) {
    cat("Starting from the tmp file.\n")
    result <- readRDS(fn.tmp)
    start <- max(which(!is.na(result$est[,1]))) + 1
  } else {
    # otherwise, build from skeleton.
    # skeleton for sim
    cat("Starting from an empty object.\n")
    
    result <- 
      list( # i: param, j: n
        setting.summary = c(i = i, j = j, `n (j)` = n, n.sim = n.sim),
        `param (i)` = c(param, rho = rho, zi = zi),
        est = data.frame(matrix(NA, n.sim, 2 * (9 + 2) + 8, #9: abp, 2: rho & logit.rho, 6: ci
                                dimnames = est.names))
      )
    start = 1L
  }
  
  for (k in start:n.sim) {
    # set.seed
    sn <- i + j * 1000 + k * 100000
    set.seed(sn)
    tmp.bgn <- Sys.time()
    xy <- rbzinb(n = n, param = param)
    tmp <- bzinb(xy[,1], xy[,2], maxiter = maxiter, initial = param, showFlag=FALSE)
    print(tmp$rho[1, ])
    result$est[k, 1:(2 * (9 + 2))] <- c(tmp$coefficients[,1], tmp$rho[,1], tmp$coefficients[,2], tmp$rho[,2])
    CI.naive <- ci(tmp$rho[1,1], tmp$rho[1,2], logit = FALSE)
    CI.logit <- ci(tmp$rho[2,1], tmp$rho[2,2], logit = TRUE)
    result$est[k, c("lb.95", "ub.95", "cov.95")] <- c(CI.logit, ifelse(prod(CI.logit - rho) < 0, 1, 0))
    result$est[k, c("lb.n.95", "ub.n.95", "cov.n.95")]  <- c(CI.naive, ifelse(prod(CI.naive - rho) < 0, 1, 0))
    result$est[k, "iter"] <- tmp$iter
    result$est[k, "time"] <- difftime(Sys.time(), tmp.bgn, units="mins")
    saveRDS(result, fn.tmp)
    # if (k %% 30 == 1) 
    cat("k = ", k, " out of ", n.sim, "\n")
  }
  
  
  result$statistic <-
    c(mean.rho = mean(result$est$rho, na.rm = TRUE),
      sd.rho = sd(result$est$rho, na.rm = TRUE),
      se.rho.mean = mean(result$est$se.rho, na.rm = TRUE),
      se.rho.median = median(result$est$se.rho, na.rm = TRUE),
      mean.logit.rho = mean(qlogis(result$est$rho), na.rm = TRUE),
      sd.logit.rho = sd(qlogis(result$est$rho), na.rm = TRUE),
      se.logit.rho.mean = mean(result$est$se.logit.rho, na.rm = TRUE),
      se.logit.rho.median = median(result$est$se.logit.rho, na.rm = TRUE),
      na.prop.rho = mean(is.na(result$est$rho)),
      na.prop.se = mean(is.na(result$est$se.rho)),  # added later
      mse.rho = mean((result$est$rho - rho)^2, na.rm = TRUE),
      CP = mean(result$est$cov.95, na.rm = TRUE),
      CIL.mean = mean(result$est$ub.95 - result$est$lb.95, na.rm = TRUE),
      CP.naive = mean(result$est$cov.n.95, na.rm = TRUE),
      CIL.naive.mean = mean(result$est$ub.n.95 - result$est$lb.n.95, na.rm = TRUE),
      mean.iter = mean(result$est$iter, na.rm = TRUE),
      mean.time = mean(result$est$time, na.rm = TRUE))
  
  saveRDS(result, fn) 
  
  result$statistic
  result$est %>% apply(2, mean, na.rm=TRUE); result$est[,1:11] %>% apply(2, sd, na.rm=TRUE)
  