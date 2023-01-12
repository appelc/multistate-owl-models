## multistate occupancy models with 2018 owl data

## create marginal plots from post-hoc model without STUDY_AREA factor covariate:
## p(NOISE + EFFORT) δ(.) ψ(NR500 + BO_TOTAL) R(NR500 + BO_TOTAL)

## CONTENTS:

## 1) View estimated and derived parameters
## 2) Explore covariate effects on occupancy (psi)
## 3) Explore covariate effects on conditional pair occupancy (R)
## 4) Calculate predicted/fitted values (occupancy)
## 5) Create marginal plots (occupancy) 

library(dplyr)
library(magrittr)
library(knitr)
library(jagsUI)
library(MCMCvis)
library(kableExtra)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

## -----------------------------------------------------------------------------
## 1) View estimated and derived parameters to report ####

## Import model
noArea_model <- readRDS('results/06_model_posthoc_norm1priors_marginals/model_output.rds')

  noArea_model$summary %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')


## -----------------------------------------------------------------------------  
## 2) Explore covariate effects on landscape use (psi) ####

## Plot posterior means and CI (50% = thick lines) (95% = thin lines):
  MCMCplot(noArea_model, params = 'beta', main = '"psi" parameters', ref_ovl = TRUE)
  #[1]: intercept, [2]: NR500, [3]: study area, [4]: BO_total

## Calculate means, SDs, and odds ratios:

  #intercepts
  psi_intMean_na  = mean(noArea_model$sims.list$beta[,1])
  psi_intSD_na    = sd(noArea_model$sims.list$beta[,1])
  psi_intQuant_na = quantile(noArea_model$sims.list$beta[,1], probs = c(0.025, 0.975))
  
  
  #forest suitability mean 500m
  psi_forestSlope_na = mean(noArea_model$sims.list$beta[,2])
  psi_forestSD_na    = sd(noArea_model$sims.list$beta[,2])
  psi_forestQuant_na = quantile(noArea_model$sims.list$beta[,2], probs = c(0.025, 0.975))
  psi_forestOR_na    = median(exp(noArea_model$sims.list$beta[,2]))       #odds ratio - use median
  psi_forestOR_LCI_na = quantile(exp(noArea_model$sims.list$beta[,2]), probs = c(0.025)) #odds ratio LCI
  psi_forestOR_UCI_na = quantile(exp(noArea_model$sims.list$beta[,2]), probs = c(0.975)) #odds ratio UCI
  
  #total barred owl
  psi_BOSlope_na = mean(noArea_model$sims.list$beta[,3])
  psi_BOSD_na    = sd(noArea_model$sims.list$beta[,3])
  psi_BOQuant_na = quantile(noArea_model$sims.list$beta[,3], probs = c(0.025, 0.975))
  psi_BOOR_na    = median(exp(noArea_model$sims.list$beta[,3]))           #odds ratio - use median
  psi_BOOR_LCI_na = quantile(exp(noArea_model$sims.list$beta[,3]), probs = c(0.025)) #odds ratio LCI
  psi_BOOR_UCI_na = quantile(exp(noArea_model$sims.list$beta[,3]), probs = c(0.975)) #odds ratio UCI
  
  psi_params_noArea <- data.frame('parameter' = rep('psi', 3),
                                  'covariate' = c('Intercept','Forest suitability mean 200m',
                                                  'Total barred owl calling'),
                                  'mean' = c(psi_intMean_na, psi_forestSlope_na, psi_BOSlope_na),
                                  'SD' = c(psi_intSD_na, psi_forestSD_na, psi_BOSD_na),
                                  'LCI_95' = c(psi_intQuant_na[[1]], psi_forestQuant_na[[1]],
                                               psi_BOQuant_na[[1]]),
                                  'UCI_95' = c(psi_intQuant_na[[2]], psi_forestQuant_na[[2]],
                                               psi_BOQuant_na[[2]]),
                                  'F_score' = noArea_model$f$beta,
                                  'Odds_ratio' = c(NA, psi_forestOR_na, psi_BOOR_na),
                                  'OR_LCI' = c(NA, psi_forestOR_LCI_na, psi_BOOR_LCI_na),
                                  'OR_UCI' = c(NA, psi_forestOR_UCI_na, psi_BOOR_UCI_na))
  ## FOR SI TABLES:  
  psi_params_noArea %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')


## -----------------------------------------------------------------------------  
## 3) Explore covariate effects on probability of pair given use (R) ####

  ## Plot posterior means and CI (50% = thick lines) (95% = thin lines):
  MCMCplot(noArea_model, params = 'beta2', main = '"R" parameters', ref_ovl = TRUE)
  #[1]: intercept, [2]: NR500, [3]: BO_total

  ## Calculate means, SDs, and odds ratios:
  
  #intercepts
  r_intMean_na  = mean(noArea_model$sims.list$beta2[,1])
  r_intSD_na    = sd(noArea_model$sims.list$beta2[,1])
  r_intQuant_na = quantile(noArea_model$sims.list$beta2[,1], probs = c(0.025, 0.975))
  
  #forest suitability weighted mean 200m
  r_forestSlope_na = mean(noArea_model$sims.list$beta2[,2])
  r_forestSD_na    = sd(noArea_model$sims.list$beta2[,2])
  r_forestQuant_na = quantile(noArea_model$sims.list$beta2[,2], probs = c(0.025, 0.975))
  r_forestOR_na    = median(exp(noArea_model$sims.list$beta2[,2]))        #odds ratio - use median
  r_forestOR_LCI_na = quantile(exp(noArea_model$sims.list$beta2[,2]), probs = c(0.025)) #odds ratio LCI
  r_forestOR_UCI_na = quantile(exp(noArea_model$sims.list$beta2[,2]), probs = c(0.975)) #odds ratio UCI
  
  #total barred owl
  r_BOSlope_na = mean(noArea_model$sims.list$beta2[,3])
  r_BOSD_na    = sd(noArea_model$sims.list$beta2[,3])
  r_BOQuant_na = quantile(noArea_model$sims.list$beta2[,3], probs = c(0.025, 0.975))
  r_BOOR_na    = median(exp(noArea_model$sims.list$beta2[,3]))            #odds ratio - use median
  r_BOOR_LCI_na = quantile(exp(noArea_model$sims.list$beta2[,3]), probs = c(0.025)) #odds ratio LCI
  r_BOOR_UCI_na = quantile(exp(noArea_model$sims.list$beta2[,3]), probs = c(0.975))     #odds ratio UCI
  
  r_params_na <- data.frame('parameter' = rep('R', 3),
                            'covariate' = c('Intercept','Forest suitability mean 200m',
                                            'Total barred owl calling'),
                            'mean' = c(r_intMean_na, r_forestSlope_na, r_BOSlope_na),
                            'SD' = c(r_intSD_na, r_forestSD_na, r_BOSD_na),
                            'LCI_95' = c(r_intQuant_na[[1]], r_forestQuant_na[[1]], r_BOQuant_na[[1]]),
                            'UCI_95' = c(r_intQuant_na[[2]], r_forestQuant_na[[2]], r_BOQuant_na[[2]]),
                            'F_score' = noArea_model$f$beta2,
                            'Odds_ratio' = c(NA, r_forestOR_na, r_BOOR_na),
                            'OR_LCI' = c(NA, exp(r_forestQuant_na[[1]]), exp(r_BOQuant_na[[1]])),
                            'OR_UCI' = c(NA, exp(r_forestQuant_na[[2]]), exp(r_BOQuant_na[[2]])))
  ## FOR SI TABLES:  
  r_params_na %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')


## -----------------------------------------------------------------------------
## 4) Calculate predicted/fitted values (use/occupancy) ####

## FOREST SUITABILITY MEAN ####
  
  #import raw covariate values and calculate mean and SD
  nr500Raw <- fread('data2018/covariates/raw/nr500_raw_2018.csv', header = TRUE)
  (zmean_nr500 <- mean(nr500Raw$MEAN))
  (zsd_nr500 <- sd(nr500Raw$MEAN))
  
  #import the scaled covariate and convert to vector
  nr500std <- fread('data2018/covariates/nr500_2018.csv')
  nr500std <- nr500std$MEAN_STD[order(nr500std$MEAN_STD)]
  
  #create a vector of values within the range of the real (scaled) values  
  nr500scaled <- seq(min(nr500std), max(nr500std), length.out = 100)
  
  #extract relevant columns from the posterior draws (intercept/s and slope for this covariate)
  simsNRpsi <- noArea_model$sims.list$beta[,c(1,2)]
  simsNRr   <- noArea_model$sims.list$beta2[,c(1,2)]
  
  #set up matrices
  l1 <- matrix(nrow = nrow(simsNRpsi), ncol = length(nr500scaled))
  l2 <- matrix(nrow = nrow(simsNRr), ncol = length(nr500scaled))
  
  #multiply intercept by 1; slope by value of the covariate
  for(i in 1:nrow(simsNRpsi)){
    for(j in 1:length(nr500scaled)){
      l1[i, j] <- sum(simsNRpsi[i, ] * c(1, nr500scaled[j])) #l1: psi
    }
  }
  for(i in 1:nrow(simsNRr)){
    for(j in 1:length(nr500scaled)){
      l2[i, j] <- sum(simsNRr[i, ] * c(1, nr500scaled[j]))   #l2: R
    }
  }
  
  #back-transform from logit
  psiBTnr <- plogis(l1)        #prob use regardless of state (psi)
  rBTsuit <- plogis(l2)        #prob of being paired (R)
  
  #calculate derived parameters
  pairBTnr  <- psiBTnr * rBTsuit        #prob occ by pairs (psi*R)
  nopairBTnr <- psiBTnr * (1 - rBTsuit) #prob occ by non-pair (psi * (1-R)) 
  unoccBTnr <- 1 - psiBTnr              #prob unoccupied (1 - psi)
  
  #convert the scaled covariate values that I used to their real scale
  nrBT <- (nr500scaled * zsd_nr500) + zmean_nr500
  
  #set up dataframe
  NRplotNoArea <- data.frame(x = rep(nrBT,4),
                             y = c(apply(psiBTnr, 2, mean),          #psi
                                   apply(pairBTnr, 2, mean),         #psi*R
                                   apply(nopairBTnr, 2, mean),       #psi*(1-R)
                                   apply(rBTsuit, 2, mean)),         #R
                             lo = c(apply(psiBTnr, 2, quantile, probs = 0.025), 
                                    apply(pairBTnr, 2, quantile, probs = 0.025),
                                    apply(nopairBTnr, 2, quantile, probs = 0.025),
                                    apply(rBTsuit, 2, quantile, probs = 0.025)),
                             hi = c(apply(psiBTnr, 2, quantile, probs = 0.975), 
                                    apply(pairBTnr, 2, quantile, probs = 0.975),
                                    apply(nopairBTnr, 2, quantile, probs = 0.975),
                                    apply(rBTsuit, 2, quantile, probs = 0.975)),
                             grp = c(rep('Landscape\nuse', length(nrBT)), rep('Pair\noccupancy', length(nrBT)),
                                     rep('Non-pair use', length(nrBT)), rep('R', length(nrBT))))


## BARRED OWL TOTAL ####

  #import raw covariate values and calculate mean and SD
  boTotalRaw <- fread('data2018/covariates/raw/bo_total_raw_2018.csv', header = TRUE)
  (zmean_boTotal <- mean(boTotalRaw$total))
  (zsd_boTotal <- sd(boTotalRaw$total))
  
  #import the scaled covariate and convert to vector
  boTotalStd <- fread('data2018/covariates/bo_total_2018.csv')
  
  #create a vector of values within the range of the real (scaled) values  
  boTotalScaled <- seq(min(boTotalStd$total_std), max(boTotalStd$total_std), length.out = 100)
  
  #extract relevant columns from the posterior draws (intercept/s and slope for this covariate)
  simsBOtotalPsi <- noArea_model$sims.list$beta[,c(1,3)]  
  simsBOtotalR   <- noArea_model$sims.list$beta2[,c(1,3)] 
  
  #set up matrices
  n1 <- matrix(nrow = nrow(simsBOtotalPsi), ncol = length(boTotalScaled))
  n2 <- matrix(nrow = nrow(simsBOtotalR), ncol = length(boTotalScaled))
  
  #multiply intercept by 1; slope by value of the covariate
  for(i in 1:nrow(simsBOtotalPsi)){
    for(j in 1:length(boTotalScaled)){
      n1[i, j] <- sum(simsBOtotalPsi[i, ] * c(1, boTotalScaled[j])) #n1: psi
    }
  }
  for(i in 1:nrow(simsBOtotalR)){
    for(j in 1:length(boTotalScaled)){
      n2[i, j] <- sum(simsBOtotalR[i, ] * c(1, boTotalScaled[j]))   #n2: R
    }
  }
  
  #back-transform from logit
  psiBTbo <- plogis(n1)         #psi
  rBTbo   <- plogis(n2)         #R

  #calculate derived parameters
  pairBTbo  <- psiBTbo * rBTbo         #prob occ by pairs (psi*R)
  unoccBTbo <- 1 - psiBTbo             #prob unoccupied (1 - psi)
  nopairBTbo <- psiBTbo * (1 - rBTbo)  #prob occ by non-pair (psi * (1-R))
  
  #convert the scaled covariate values that I used to their real scale
  boTotalBT <- (boTotalScaled * zsd_boTotal) + zmean_boTotal
  
  #set up dataframe
  boTotalPlotNoArea <- data.frame(x = rep(boTotalBT, 4),
                                  y = c(apply(psiBTbo, 2, mean),     #psi
                                        apply(pairBTbo, 2, mean),    #psi*R
                                        apply(nopairBTbo, 2, mean),  #psi*(1-R)
                                        apply(rBTbo, 2, mean)),      #R
                                  lo = c(apply(psiBTbo, 2, quantile, probs = 0.025), 
                                         apply(pairBTbo, 2, quantile, probs = 0.025),
                                         apply(nopairBTbo, 2, quantile, probs = 0.025),
                                         apply(rBTbo, 2, quantile, probs = 0.025)),
                                  hi = c(apply(psiBTbo, 2, quantile, probs = 0.975), 
                                         apply(pairBTbo, 2, quantile, probs = 0.975),
                                         apply(nopairBTbo, 2, quantile, probs = 0.975),
                                         apply(rBTbo, 2, quantile, probs = 0.975)),
                                  grp = c(rep('Landscape\nuse', length(boTotalBT)), rep('Pair\noccupancy', length(boTotalBT)),
                                          rep('Non-pair use', length(boTotalBT)), rep('R', length(boTotalBT))))
  
  

## -----------------------------------------------------------------------------
## 5) Create marginal plots (use/occupancy) ####

## NR 500
  rr <- ggplot(NRplotNoArea[NRplotNoArea$grp %in% c('Landscape\nuse','Pair\noccupancy'),], 
               aes(x/10000, y)) + 
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = grp), alpha = 0.3) +
    geom_line(size = 1.5, aes(color = grp, linetype = grp)) +
    ylab('Probability \u00B1 95 CI') + xlab('Mean NR suitability index (500m)') +
    ylim(c(0,1)) +
    # ylim(c(0,0.7)) +
    scale_fill_manual(values = c('Landscape\nuse' = 'black', 'Pair\noccupancy' = 'darkblue')) +
    scale_colour_manual(values = c('Landscape\nuse' = 'black', 'Pair\noccupancy' = 'darkblue')) +
    scale_linetype_manual(values = c('Landscape\nuse' = 'solid', 'Pair\noccupancy' = 'twodash')) +
    geom_rug(data = nr500Raw, mapping = aes(x = MEAN/10000), inherit.aes = FALSE) + 
    theme(panel.background = element_rect(fill = 'transparent'),
          axis.line = element_line(),
          axis.title = element_text(size = 18),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 16),
          legend.title = element_blank(),
          #legend.position = 'none',
          legend.background = element_rect(fill='transparent')
    )   
  rr
 

## BARRED OWL TOTAL
bb <- ggplot(boTotalPlotNoArea[boTotalPlotNoArea$grp %in% c('Landscape\nuse','Pair\noccupancy'),], aes(x, y)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = grp), alpha = 0.3) +
  geom_line(size = 1.5, aes(color = grp, linetype = grp)) +
  ylab('Probability \u00B1 95 CI') + xlab('Barred owl detections (total)') +
  # scale_x_continuous(breaks = seq(0, 8000, 2000), limits = c(0, 5000)) +
  ylim(c(0,1)) +
  # ylim(0,0.7) +
  scale_fill_manual(values = c('Landscape\nuse' = 'black', 'Pair\noccupancy' = 'darkblue')) +
  scale_colour_manual(values = c('Landscape\nuse' = 'black', 'Pair\noccupancy' = 'darkblue')) +
  scale_linetype_manual(values = c('Landscape\nuse' = 'solid', 'Pair\noccupancy' = 'twodash')) +
  geom_rug(data = boTotalRaw, mapping = aes(x = total), inherit.aes = FALSE) +
  theme(panel.background = element_rect(fill = 'transparent'),
        axis.line = element_line(),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        #legend.position = 'none',
        legend.background = element_rect(fill='transparent')
  )   
bb

## Combine and export
  Fig4 <- ggarrange(rr + rremove("ylab"), bb + rremove("ylab"), 
                    labels = c("a", "b"), font.label = list(size = 24), 
                    vjust = c(1.5,1.5), hjust = c(-4,-3.5), 
                    ncol = 1, nrow = 2, 
                    common.legend = TRUE, legend = 'right') +
    theme(plot.margin = margin(0.1,0.1,0.1,1, "cm"))

  Fig4_shared_axis <- annotate_figure(Fig4,
                                      left = textGrob('Probability \u00B1 95 CI',
                                                      rot = 90, vjust = 1,
                                                      gp = gpar(cex = 1.6)))

  ## FOR FIGURE 3
  tiff(filename = 'figures/fig4.tif', height = 18, width = 16, units = 'cm',
       res = 600, compression = 'lzw')
  print(Fig4_shared_axis)
  dev.off()  


## -----------------------------------------------------------------------------

