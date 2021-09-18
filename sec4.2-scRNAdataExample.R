##############################################################################
# Codes for Section 4. Model and measure comparisons based on real data
##############################################################################

#### 0. library #### 
# functions and libraries
  # install.packages("bzinb")
  library(bzinb)
  library(tidyverse); library(magrittr)
  library(ggplot2); library(ggrepel)
  library(grid); library(gridExtra)
  
  source("F.sec4-functions.R")
  
  # data
  data <-read.delim("../../Data/sox9ko_paneth_800_cell_counts.txt" , sep="\t")

#### 1. distribution #### 

# distribution of zero-proportions
  zero.prop <-
    data[, -1] %>% 
    apply(1, function(x) {mean(x==0)})

  zero.prop %>% hist
  zero.prop %>% (function(x) {mean(x >= .8)}) # % genes with more than 80% zeros = 98.3%
  zero.prop %>% (function(x) {mean(x >= .9)}) # % genes with more than 90% zeros = 93.8%
  zero.prop %>% mean                          # mean zero-proportions            = 97.3%

  
#### 2. stratefied sampling according to zero-proportion #### 
  ## 2.1 grouping data with zero-proportion
  zero.prop <- data.frame(id = 1:length(zero.prop),
                          prop = zero.prop,
                          class = cut(zero.prop, breaks = c(0, 0.6, 0.8, 0.9, 0.98, 1), 
                                      labels = c("very low", "low", "moderate", "high", "very high"), 
                                      include.lowest = TRUE))
  zero.prop$max.expr <- apply(data[, -1], 1, max)
  
  zero.prop$class %>% table %>% print %>% "/"(sum(.)) %>% round(4) 
    # very low      83    0.4%
    # low          327    1.4%
    # moderate   1,060    4.5%
    # high       6,375   27.2%
    # very high 15,580   66.5%
  
  {zero.prop %>% filter(max.expr <= 10000, class != "very high")}$class %>% 
    table %>% print %>% "/"(sum(.)) %>% round(4) 
    # very low      73    0.3%
    # low          326    1.4%
    # moderate   1,060    4.5%
    # high       6,375   27.2%
    # very high 15,580   66.5%

  n.candidate = 100 # number of genes to be considered for each class
  n.sample = 5  # final sample size for each class
  pairs <- data.frame(geneA = rep(NA, n.sample * 10), 
                      geneB = NA,
                      class = rep(c("HH", "MM", "LL", "VV", 
                                    "HM", "HL", "HV", "ML", "MV", "LV"), 
                                  each = n.sample),
                      propA = NA, propB = NA) 
  pairs$nameA = paste0(1:50, "A (", substr(pairs$class, 1, 1), ")")
  pairs$nameB = paste0(1:50, "B (", substr(pairs$class, 2, 2), ")")


  ## 2.2 sampling candidates
  set.seed(1)
  high <- zero.prop %>% 
    dplyr::filter(class == "high" & max.expr <= 10000) %>% 
    sample_n(n.candidate)
  
  set.seed(2)
  moderate <- zero.prop %>% 
    dplyr::filter(class == "moderate" & max.expr <= 10000) %>% 
    sample_n(n.candidate)
  
  set.seed(3)
  low <- zero.prop %>% 
    dplyr::filter(class == "low" & max.expr <= 10000) %>% 
    sample_n(n.candidate)
  
  set.seed(4)
  Vlow <- zero.prop %>% 
    dplyr::filter(class == "very low" & max.expr <= 10000) %>% 
    sample_n(40)

  comb = tibble(comb = 1:10, 
                name = c("HH", "MM", "LL", "VV", "HM", "HL", "HV", "ML", "MV", "LV"),
                setA = substr(name, 1, 1),
                setB = substr(name, 2, 2))

  
  ## 2.3 sampling the final data
  for (l in 1:10) {
    cat("finding pairs from ", unlist(comb[l, "name"]), "\n")
    datA <- switch(comb[l, "setA"] %>% unlist, H = high, M = moderate, L = low, V = Vlow)
    datB <- switch(comb[l, "setB"] %>% unlist, H = high, M = moderate, L = low, V = Vlow)
    nA = dim(datA)[1]
    nB = dim(datB)[1]
    sameSet = comb[l, "setA"] == comb[l, "setB"]
    
    k = 1
    
    for (i in 1:nA) {
      for (j in 1:nB) {
        if (sameSet) {if (j <= i | i == nA) next}
        idA = datA$id[i]
        idB = datB$id[j]
        if (idA %in% c(pairs$geneA, pairs$geneB)) break
        if (idB %in% c(pairs$geneA, pairs$geneB)) next
        xvec <- data[idA, -1] %>% as.numeric
        yvec <- data[idB, -1] %>% as.numeric
        {contingency <- bzinb:::bin.profile(xvec, yvec)} # %>% print # both nonzero / x and y only nonzero / both zero
        {both.nonzero.prop <- contingency[1] / sum(contingency[1:3])} # %>% print
        # if (both.nonzero.prop >= 0.05 & contingency[1] >= 4) {  #second = [/ge 0.005]
          cat(i, "-", j, " pair selected. both nonzero = ", contingency[1], " out of ", sum(contingency[1:3]), "bothnonzero % =", both.nonzero.prop, "\n")
          pairs[k + n.sample * (l-1), "geneA"] <- idA
          pairs[k + n.sample * (l-1), "geneB"] <- idB
          pairs[k + n.sample * (l-1), "propA"] <- datA$prop[i]
          pairs[k + n.sample * (l-1), "propB"] <- datB$prop[j]
          k <- k + 1
          break # avoid multiple j's for the same i.
        #}
      }
      if (k > n.sample) break
    }
  }
  pairs
  pairs$pair.name = paste0(pairs$class, 1:5)
  # filename.pairs = "../output/C0201.pairs.rds"
  # saveRDS(pairs, filename.pairs)
  # pairs <- readRDS(measure.est.filename)


#### 3. Density estimation (PC, SC, tau, MI, rho*, MI*) #### 

  # skeleton  
  measure.i = setNames(vector("list", dim(pairs)[1]), 1:dim(pairs)[1])
  measure.est.filename = "../../output/C0201.measure.est.rds"
  
  for (i in 1:dim(pairs)[1]) {  
    cat("i = ", i, "\n")
    xvec <- data[pairs[i, "geneA"], -1] %>% as.numeric
    yvec <- data[pairs[i, "geneB"], -1] %>% as.numeric
    measure.i[[i]] <- measures(xvec, yvec)
    parm.BZNB <- measure.i[[i]]$parameters[, 1]
    
    # adding other model parameters
    # BP
    measure.i[[i]]$param.BP <- bp(xvec, yvec)
    
    # BZIP
    parm.BZIP <- bzip.b(xvec, yvec, maxiter = 10000)
    
    # parm.BZIP[7] <- 1 - sum(parm.BZIP[4:6])    # overcome small numerical errors.
    measure.i[[i]]$param.BZIP <- parm.BZIP
    
    # BNB
    measure.i[[i]]$param.BNB <- bnb(xvec, yvec)
    
    # saveRDS(measure.i, measure.est.filename)
    gc()
  }
  saveRDS(measure.i, measure.est.filename)
  # measure.i <- readRDS(measure.est.filename)


#### 4. Plots for model comparisons #### 
  
  ## 4.1 graphical params
  xlim = c(0,80)
  ylim = c(0,100)
  ratio = ylim[2]/xlim[2]
  midpoint = 0.4
  size = 2
  fig_folder = "../../figure/realData_20190626/"
  fig_folder2 = "../../figure/realData_20190626/raw/"
  if (!dir.exists(fig_folder)) dir.create(fig_folder)
  if (!dir.exists(fig_folder2)) dir.create(fig_folder2)
  
  ## 4.2 generating raw plots
  for (seed.no in 1:3) {
    for (i in 1:dim(pairs)[1]) {
      cat("i: ", i, "\n")
      xvec <- data[pairs[i, "geneA"], -1] %>% as.numeric
      yvec <- data[pairs[i, "geneB"], -1] %>% as.numeric
      
      # 3.2.1 estimation
      parm.BP <- measure.i[[i]]$param.BP[1:3]
      parm.BZIP <- measure.i[[i]]$param.BZIP[1:7]
      parm.BNB <- measure.i[[i]]$param.BNB$coefficients[, 1]
      parm.BZNB <- measure.i[[i]]$parameters[,1]
      
      # 3.2.2 Empirical distn
      densitymap(xvec, yvec) -> a.emp
      a.emp %>% arrange(desc(x+y)) %>% 
        ggplot(aes(x,y,col=density)) + geom_point(size=size, shape=15) + xlim(xlim) + ylim(ylim) + ggtitle(NULL) +
        scale_colour_gradient2(limits = c(0,1), midpoint = midpoint, low = "white", mid = "grey", high = "black") +
        xlab(bquote(.(paste0("Gene ", pairs[i, "nameA"]," (")) ~ Y[1] ~.(")"))) + 
        ylab(bquote(.(paste0("Gene ", pairs[i, "nameB"]," (")) ~ Y[2] ~.(")"))) + coord_equal() + 
        theme_bw() + theme(panel.border = element_rect(colour = "gray"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #panel.border = element_blank(), 
                           plot.margin = rep(unit(0, "null"),4), panel.margin = unit(0.1,"cm"),
                           legend.position="bottom") +
        ggtitle(pairs$pair.name[i]) -> pr0
      ggsave(paste0(fig_folder2, "plot4-1-",i,"a.png"), pr0, width = 2.7, height = 4, dpi = 150)
      
      
      # 3.2.3 Empirical distns of models' random realizations
      # random data generation for each model and combine.
      a.tmp.comb = data.frame(x = NULL, y = NULL, freq = NULL, density = NULL, sim = NULL)
      set.seed(seed.no)
      rbp(800, param = parm.BP[1:3]) -> a.tmp
      a.tmp.comb = rbind(a.tmp.comb, cbind(densitymap(a.tmp[,1], a.tmp[,2])%>% arrange(desc(x+y)), sim = "1. BP"))
      set.seed(seed.no)
      rbzinb(800, param = c(parm.BNB %>% as.numeric,1,0,0,0)) -> a.tmp
      a.tmp.comb = rbind(a.tmp.comb, cbind(densitymap(a.tmp[,1], a.tmp[,2])%>% arrange(desc(x+y)), sim = "3. BNB"))
      set.seed(seed.no)
      rbzip.b(800, param = parm.BZIP) -> a.tmp
      a.tmp.comb = rbind(a.tmp.comb, cbind(densitymap(a.tmp[,1], a.tmp[,2])%>% arrange(desc(x+y)), sim = "2. BZIP"))
      set.seed(seed.no)
      rbzinb(800, param = parm.BZNB) -> a.tmp
      a.tmp.comb = rbind(a.tmp.comb, cbind(densitymap(a.tmp[,1], a.tmp[,2])%>% arrange(desc(x+y)), sim = "4. BZINB"))
      
      a.tmp.comb %>% 
        ggplot(aes(x,y,col=density)) + geom_point(size=size, shape=15) + xlim(xlim) + ylim(ylim) + ggtitle(NULL) +
        facet_wrap(~sim, 2, 2) +
        scale_colour_gradient2(limits = c(0,1), midpoint = midpoint, low = "white", mid = "grey", high = "black") +
        xlab(NULL) + #xlab(bquote(.("Gene 1 (") ~ Y[1] ~.(")"))) + 
        ylab(NULL) + #ylab(bquote(.("Gene 2 (") ~ Y[2] ~.(")"))) +
        coord_equal() + guides(col=FALSE) + theme_bw() + 
        theme(panel.border = element_rect(colour = "gray"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              plot.margin = rep(unit(0.05,"null"),4), panel.margin = unit(0.1,"cm")) -> p
    #  ggsave(paste0(fig_folder2, "plot4-1-",i,"b_seed",seed.no, ".png"), p, width = 4, height = 5, dpi = 500)
      
      
      # 3.2.4 Putting them all together
      grid.arrange(
        grobs = list(pr0, p),
        widths = c(2, 3)
      ) -> p
      ggsave(paste0(fig_folder2, "Figure1-",i,"_seed",seed.no, ".png"), p, width = 6, height = 5, dpi = 600)
      saveRDS(p, paste0(fig_folder2, "rds-Figure1-",i,"_seed",seed.no, ".rds"))
      gc()
    }
  }
  
  ## 4.3 combining plots
  # Figure1. 10 figures from first pair of each combination (only for seed 1)
  fig.tab <- lapply(1:10, function(i) {
    readRDS(paste0(fig_folder2, "rds-Figure1-",5*i - 4,"_seed1.rds")) %>% 
      gtable::gtable_add_grob(rectGrob(gp=gpar(lwd=5, fill=NA)), 1, 1, 1, 2)})
  grid.arrange(
    grobs = fig.tab,
    ncol = 4, nrow = 3
  ) -> p
  ggsave(paste0(fig_folder, "Figure1.png"), p, width = 20, height = 15, dpi = 300)
  
  # Appendix-Figure1A ~ 1C. 10 x 5 table figure (seed 1, 2, 3)
  for (seed.no in 1:3) {
    fig.tab <- lapply(1:50, function(i) {
      readRDS(paste0(fig_folder2, "rds-Figure1-",i,"_seed", seed.no, ".rds")) %>% 
        gtable::gtable_add_grob(rectGrob(gp=gpar(lwd=5, fill=NA)), 1, 1, 1, 2)})
    grid.arrange(
      grobs = fig.tab,
      ncol = 5, nrow = 10
    ) -> p
    ggsave(paste0(fig_folder, "Figure1app_seed", seed.no, ".png"), p, width = 25, height = 35, dpi = 300)
  }
  
  
  
#### 5. Measure comparisons  ####
  ## measure comparison figure
  xlim = c(0,100) ## 100 now!!
  ylim = c(0,100)
  midpoint = 0.4
  ratio = 1
  measure.plots <- vector("list", dim(pairs)[1])
  for (i in 1:dim(pairs)[1]) {
    xvec = data[pairs[i, "geneA"],-1] %>% as.numeric
    yvec = data[pairs[i, "geneB"],-1] %>% as.numeric
    densitymap(xvec, yvec) -> a.emp
    # measures(xvec, yvec) %>% round(3) -> measure.tmp
    measure.tmp <- measure.i[[i]]$measure
    measure.tmp <- cbind(measure.tmp, lb = measure.i[[i]]$rho.ci[1], ub  = measure.i[[i]]$rho.ci[2])
    measure.tmp %<>% format(digits = 1, nsmall=3)
    measure.tmp <- c(measure.tmp[1:5], 
                     rho.ci = paste0("(", measure.tmp$lb, ", ", measure.tmp$ub, ")"),
                     measure.tmp[6])
    a.emp %>% arrange(desc(x+y)) %>% 
      ggplot(aes(x,y,col=density), col = "turquoise4") + 
      geom_point(size=size, shape=15) + xlim(xlim) + ylim(ylim) + 
      # ggtitle(paste0("Pair ", i)) +
      scale_colour_gradient2(limits = c(0,1), midpoint = midpoint, low = "white", mid = "grey", high = "black") +
      guides(col = FALSE) + 
      xlab(bquote(.(paste0("Gene ", pairs[i, "nameA"]," (")) ~ Y[1] ~.(")"))) + 
      ylab(bquote(.(paste0("Gene ", pairs[i, "nameB"]," (")) ~ Y[2] ~.(")"))) + coord_equal() + 
      annotate("text", x = 40, y = 100 - (1:7)*8, hjust = 0, size = 2.7,
               label = paste0(c("PC", "SC", "tau", "EMI", "rho~{'*'}", "", "MI~{'*'}")), parse = TRUE) +
      annotate("text", x = 55, y = 100 - (1:7)*8, hjust = 0, size = 2.7,
               label = paste0(c(rep("=", 5), "", "=")), parse = FALSE) +
      annotate("text", x = 90, y = 100 - (1:7)*8, hjust = 1, size = 2.7,
               label = measure.tmp) +
      annotate("text", x = 15, y = 90, size = 5, label = pairs$pair.name[i]) +
      theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.border = element_rect(colour = "gray"),
                         plot.margin = rep(unit(0,"null"),4),
                         panel.margin = unit(0,"null")) -> measure.plots[[i]]
    
    ggsave(paste0(fig_folder2, "Figure2-", i,".png"), width = 2, height = ratio*2, dpi = 600)
  }
  saveRDS(measure.plots, paste0(fig_folder2, "rds-Figure2.rds"))
  grid.arrange(
    grobs = measure.plots,
    ncol = 5, nrow = 10
  ) -> p
  ggsave(paste0(fig_folder, "Figure2App.png"), p, width = 15, height = 25, dpi = 300, limitsize = FALSE)
  
  
  measure.summary <- as.data.frame(t(sapply(measure.i, function(s) s$measure)))
  measure.summary <- as.data.frame(sapply(measure.summary, unlist))
  names(measure.summary) <- c("PC", "SC", "tau", "EMI", "rho*", "MI*")
  measure.summary$class <- rep(c("HH", "MM", "LL", "VV", "HM", "HL", "HV", "ML", "MV", "LV"), each = n.sample)
  measure.summary$pair <- paste0(measure.summary$class, 1:5)
  
  #PC vs rho*
  measure.summary %>% ggplot(aes(PC, `rho*`, col = class)) + ylab(bquote(rho~{'*'})) + 
    geom_smooth(method = "lm", se = FALSE, col = "grey", fullrange = TRUE) +
    geom_point(aes(PC, `rho*`), col = "red", shape = 1, size = 10, stroke = 1,
               data = measure.summary %>% filter(pair %in% c("HL4", "HL5"))) +
    geom_point() + guides(col = FALSE) +
    geom_text_repel(mapping = aes(label = pair), size = 3) -> p1
  ggsave(paste0(fig_folder2, "Figure2B-1.png"))
  
  #MI vs MI*
  measure.summary %>% ggplot(aes(EMI, `MI*`, col = class)) + ylab("MI*") + 
    geom_smooth(method = "lm", se = FALSE, col = "grey", fullrange = TRUE) +
    geom_point(aes(EMI, `MI*`), col = "red", shape = 1, size = 10, stroke = 1,
               data = measure.summary %>% filter(pair %in% c("MV1"))) +
    geom_point() + guides(col = FALSE) +
    geom_text_repel(mapping = aes(label = pair), size = 3) -> p2
  ggsave(paste0(fig_folder2, "Figure2B-2.png"))
  
  grid.arrange(
    grobs = list(p1, p2),
    width = c(1, 1),
    layout_matrix = matrix(c(1, 2) , 1, byrow = TRUE)
  ) -> p
  ggsave(paste0(fig_folder, "Figure2.png"), p, width = 8, height = 4, dpi = 600)
