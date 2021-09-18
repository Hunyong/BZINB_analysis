#### 0. library ####
  library(bzinb)
  library(tidyverse); library(magrittr)
  source("F.sec4-functions.R")
  sampsize = c(250, 500, 800, 1500, 2500)
  
  
#### 1. parameters ####
  ## 1.1 alpha-beta
  ab <- c(0.01, 0.1, 1.0, 0.5, 0.5,
          0.05, 2.0, 3.0, 3.0, 0.5,
          0.2, 0.3, 3.0, 2.0, 1.5,
          0.5, 2.0, 2.0, 0.5, 3.0,
          1.0, 1.0, 1.0, 1.5, 1.5,
          3.0, 2.0, 1.0, 1.5, 0.5,
          0.2, 0.05, 0.05, 3.0, 3.0,
          2.0, 0.7, 0.1, 2.5, 2.5)
  ab <- matrix(ab, 5) %>% t %>% data.frame  
  names(ab) <- c("a0", "a1", "a2", "b1", "b2")
  ab$rho <- apply(ab, 1, function(s) do.call(RHO, as.list(s)))
  
  ## 1.2 pi
  pp <- c(0.7, 0.1, 0.1, 0.1,
          0.5, 0.15, 0.15, 0.2,
          0.5, 0.1, 0.3, 0.1,
          0.2, 0.2, 0.2, 0.4,
          0.2, 0.1, 0.4, 0.3)
  pp <- matrix(pp, 4) %>% t %>% data.frame  
  names(pp) <- c("p1", "p2", "p3", "p4")
  pp$zi <- c("low", "mid-B", "mid-U", "high-B", "high-U")
  
  ## 1.3 combination of ab & pi
  abp <- expand.grid(1:nrow(ab), 1:nrow(pp))
  abp <- cbind(ab[abp$Var1,], pp[abp$Var2,])
  abp <- abp %>% 
    mutate(id = 1:n()) %>% 
    select(id, rho, zi, everything())