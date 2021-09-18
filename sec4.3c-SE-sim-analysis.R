### 0. library
  source("sec4.3a-SE-sim-parm.R")
  library(ggplot2); library(gridExtra)

  ## file names
  path <- "../../output/SE-MI50K/result-"; 
  tmp.path <- gsub("result", "tmp.result", path)
  rds.name <- ".sim"
  
  cp.type <- "logit"
  # cp.type <- "naive"
  fig.name1 <- paste0("../figure/SE-50K-",cp.type, "-CP.png")
  fig.name2 <- paste0("../figure/SE-50K-",cp.type, "-SE.png")
  
  ## skeletons
  tmp2 <-
    cbind(abp[, c("id", "rho", "zi")], 
                matrix(NA, nrow = dim(abp)[1], ncol = length(sampsize), 
                       dimnames = list(NULL, paste0("n.", sampsize))))
  stat <- list(mean.rho = tmp2,
               sd.rho = tmp2,
               se.rho.mean = tmp2,
               mean.logit.rho = tmp2,
               sd.logit.rho = tmp2,
               se.logit.rho.mean = tmp2,
               na.prop.rho = tmp2,
               na.prop.se = tmp2,
               mse.rho = tmp2,
               CP = tmp2, 
               CP.naive = tmp2, 
               CIL = tmp2,
               mean.iter = tmp2,
               mean.time = tmp2)

### 1. collecting statistics
  
  for (i in 1:40) {
    for (j in 1:5) {
      if (file.exists(paste0(path, i, "-", j , rds.name))) {
        tmp <- readRDS(paste0(path, i, "-", j , rds.name))$statistic
        stat$mean.rho[i, 3 + j] <- tmp["mean.rho"]
        stat$sd.rho[i, 3 + j] <- tmp["sd.rho"]
        stat$se.rho.mean[i, 3 + j] <- tmp["se.rho.mean"]
        stat$mean.logit.rho[i, 3 + j] <- tmp["mean.logit.rho"]
        stat$sd.logit.rho[i, 3 + j] <- tmp["sd.logit.rho"]
        stat$se.logit.rho.mean[i, 3 + j] <- tmp["se.logit.rho.mean"]
        stat$na.prop.rho[i, 3 + j] <- tmp["na.prop.rho"]
        stat$na.prop.se[i, 3 + j] <- tmp["na.prop.se"]
        stat$mse.rho[i, 3 + j] <- tmp["mse.rho"]
        stat$CP[i, 3 + j] <- tmp["CP"]
        stat$CP.naive[i, 3 + j] <- tmp["CP.naive"]
        stat$CIL[i, 3 + j] <- tmp["CIL"]
        stat$mean.iter[i, 3 + j] <- tmp["mean.iter"]
        stat$mean.time[i, 3 + j] <- tmp["mean.time"]
      }
    }
  }
  stat$CP
  stat$CP.naive

  ## longform transformation
  stat2 <- Reduce(function(...) merge(..., by=c("id", "rho", "zi", "subgroup", "n"), all.x=TRUE), 
                  map2(stat, names(stat), longform))

  ## create separate name vectors
  # run `demo(plotmath)` for more examples of mathematical annotation in R
  zi_names <- c(
    `low` = "{'i. '}~pi==(0.7~0.1~0.1~0.1)",
    `mid-B` = "{'ii. '}~pi==(0.5~0.15~0.15~0.2)",
    `mid-U` = "{'iii. '}~pi==(0.5~0.1~0.3~0.1)",
    `high-B` = "{'iv. '}~pi==(0.2~0.2~0.2~0.4)",
    `high-U` = "{'v. '}~pi==(0.2~0.1~0.4~0.3)"
  )
  rho_names <- c(
    `0.01` = "{'4. '}~rho==0.01",
    `0.1` = "{'3. '}~rho==0.1",
    `0.3` = "{'2. '}~rho==0.3",
    `0.6` = "{'1. '}~rho==0.6"
  )
  
### 2. plots      
  stat2 %>% 
    ggplot(aes(x = n, y = mean.rho, col = rho:subgroup, shape = subgroup, group = id)) +
    facet_grid(rho ~ zi, scale = "free_y",
               labeller = labeller(rho = as_labeller(rho_names, label_parsed),
                                   zi = as_labeller(zi_names, label_parsed))
               ) + guides(col = FALSE) +
    geom_line() + geom_point() + 
    geom_line(aes(x = n, y = as.numeric(as.character(rho)), col = rho:subgroup), linetype = 2) +
    scale_x_discrete(breaks = sampsize, labels = parse(text = paste0("n == ", sampsize))) + #sampsize = c(250, 500, 800, 1500, 2500)
    theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.title = element_blank(), legend.position = "bottom") + 
    xlab ("sample size") + ylab(expression({"mean"}~hat(rho)~{"(line),   CP (number)"})) +
    # 1. labels of mean(CP)
    geom_text(aes(x = as.numeric(n) + ifelse(subgroup=="-a", -0.4, +0.4), 
                  y = mean.rho + ifelse(subgroup=="-a", -0.04, +0.04) * ifelse(rho %in% c(0.1, 0.3), -1, 1), 
                  col = rho:subgroup, label =  formatC(if(cp.type == "logit") stat2$CP else stat2$CP.naive, digits = 2, format = 'f', flag= '0')),
              position = ggstance::position_dodgev(height = 0.03), size =3, angle = 90) #+
    
  ggsave(fig.name1, width = 20, height = 20, units="cm")


# technical summary
  stat2 %>% 
    ggplot(aes(x = n, col = rho:subgroup, shape = subgroup, group = id)) +
    # SE
    geom_point(aes(x = n, y = se.rho.mean, col = rho:subgroup)) +
    geom_line(aes(x = n, y = se.rho.mean, col = rho:subgroup), linetype = 1) +
    # SD
    geom_point(aes(x = n, y = sd.rho, col = rho:subgroup)) +
    geom_line(aes(x = n, y = sd.rho, col = rho:subgroup), linetype = 2) +
    facet_grid(rho ~ zi,
               labeller = labeller(rho = as_labeller(rho_names, label_parsed),
                                   zi = as_labeller(zi_names, label_parsed))
    ) + guides(col = FALSE) + 
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_discrete(breaks = sampsize, labels = parse(text = paste0("n == ", sampsize))) + #sampsize = c(250, 500, 800, 1500, 2500)
    theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.title = element_blank(), legend.position = "bottom") + 
    xlab ("sample size") + ylab("SE (solid),   SD (dashed)")
  ggsave(fig.name2, width = 20, height = 20, units="cm")

