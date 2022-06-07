## multistate occupancy models with 2018 owl data

## evaluate full model for inference on covariate effects:
## p(NOISE + EFFORT) δ(.) ψ(NR500 + STUDY_AREA + BO_TOTAL) R(NR500 + STUDY_AREA + BO_TOTAL)

## CONTENTS:
## 1) Calculate derived parameters and compare with naive estimates
## 2) Explore covariate effects on detection (calculate odds-ratios, etc.)
## 3) Calculate predicted/fitted values (detection)
## 4) Create marginal plots (detection) 
## 5) Explore covariate effects on occupancy (psi)
## 6) Explore covariate effects on conditional pair occupancy (R) 
## 7) Calculate predicted/fitted values (occupancy)
## 8) Create marginal plots (occupancy)

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
library(stringr)
library(ggpubr)

## -----------------------------------------------------------------------------
## 1) View estimated and derived parameters to report ####

## Import model
inference_model <- readRDS('results/05a_model_inference_norm1priors/model_output.rds')
# inference_model <- readRDS('C:/users/caral/Documents/_RESEARCH/Models/results2018_Jan2022/06_model_reduced_norm1priors/model_output.rds')

  ## FOR SI TABLES:
  inference_model$summary %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')

  
## -----------------------------------------------------------------------------
## 2) Explore covariate effects on detection  ####
  
## Plot posterior means and CI (50% = thick lines) (95% = thin lines):
  MCMCplot(inference_model, params = 'gamma', main = '"p" parameters', ref_ovl = TRUE)
  #[1]: intercept, [2]: intercept adjustment for pairs, [3]: noise, [4]: effort
  
## Calculate means, SDs, and odds ratios:
  
  #intercepts
  p_intMean_inf  = mean(inference_model$sims.list$gamma[,1])
  p_intSD_inf    = sd(inference_model$sims.list$gamma[,1])
  p_intQuant_inf = quantile(inference_model$sims.list$gamma[,1], probs = c(0.025, 0.975))
  
  p_intPairMean_inf  = mean(inference_model$sims.list$gamma[,2])
  p_intPairSD_inf    = sd(inference_model$sims.list$gamma[,2])
  p_intPairQuant_inf = quantile(inference_model$sims.list$gamma[,2], probs = c(0.025, 0.975))
  
  #noise
  p_noiseSlope_inf = mean(inference_model$sims.list$gamma[,3])
  p_noiseSD_inf    = sd(inference_model$sims.list$gamma[,3])
  p_noiseQuant_inf = quantile(inference_model$sims.list$gamma[,3], probs = c(0.025, 0.975))
  p_noiseOR_inf    = median(exp(inference_model$sims.list$gamma[,3]))      #odds ratio - use median
  p_noiseOR_LCI_inf = quantile(exp(inference_model$sims.list$gamma[,3]), probs = c(0.025)) #odds ratio LCI
  p_noiseOR_UCI_inf = quantile(exp(inference_model$sims.list$gamma[,3]), probs = c(0.975)) #odds ratio UCI
  
  #effort
  p_effortSlope_inf = mean(inference_model$sims.list$gamma[,4])
  p_effortSD_inf    = sd(inference_model$sims.list$gamma[,4])
  p_effortQuant_inf = quantile(inference_model$sims.list$gamma[,4], probs = c(0.025, 0.975))
  p_effortOR_inf    = median(exp(inference_model$sims.list$gamma[,4]))    #odds ratio - use median
  p_effortOR_LCI_inf = quantile(exp(inference_model$sims.list$gamma[,4]), probs = c(0.025))   #odds ratio LCI 
  p_effortOR_UCI_inf = quantile(exp(inference_model$sims.list$gamma[,4]), probs = c(0.975))   #odds ratio UCI
  
  p_params_inference <- data.frame('parameter' = rep('p', 4),
                                'covariate' = c('Intercept','Intercept for pairs',
                                                'Background noise', 'Recording effort'),
                                'mean' = c(p_intMean_inf, p_intPairMean_inf, p_noiseSlope_inf, p_effortSlope_inf),
                                'SD' = c(p_intSD_inf, p_intPairSD_inf, p_noiseSD_inf, p_effortSD_inf),
                                'LCI_95' = c(p_intQuant_inf[[1]], p_intPairQuant_inf[[1]], 
                                             p_noiseQuant_inf[[1]], p_effortQuant_inf[[1]]),
                                'UCI_95' = c(p_intQuant_inf[[2]], p_intPairQuant_inf[[2]], 
                                             p_noiseQuant_inf[[2]], p_effortQuant_inf[[2]]),
                                'F_score' = inference_model$f$gamma,
                                'Odds_ratio' = c(NA, NA, p_noiseOR_inf, p_effortOR_inf), 
                                'OR_LCI' = c(NA, NA, p_noiseOR_LCI_inf, p_effortOR_LCI_inf),
                                'OR_UCI' = c(NA, NA, p_noiseOR_UCI_inf, p_effortOR_UCI_inf))

  ## FOR TABLE 4:
  p_params_inference %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')
  

## -----------------------------------------------------------------------------
## 3) Calculate predicted/fitted values (detection) ####
  
## NOISE ####

  #import raw covariate values and calculate mean and SD
  noiseRaw <- fread('data2018/covariates/raw/noise_raw_2018.csv', header = TRUE)
  noiseRaw <- data.frame('noise' = unlist(noiseRaw[,-1]))
    (zmean_noise <- mean(noiseRaw$noise, na.rm = TRUE))
    (zsd_noise   <- sd(noiseRaw$noise, na.rm = TRUE))
  
  #import the scaled covariate and convert to vector
  noiseStd <- fread('data2018/covariates/noise_2018.csv')
  noiseStd <- noiseStd[,-c(1,2)]
    noiseStd <- unlist(noiseStd, use.names = F)
    noiseStd <- noiseStd[order(noiseStd)]
  
  #create a vector of values within the range of the real (scaled) values  
  noiseScaled <- seq(min(noiseStd), max(noiseStd), length.out = 100)  
  
  #extract relevant columns from the simulations list (intercept/s and slope for this covariate)
  simsNoise <- inference_model$sims.list$gamma[,c(1:3)]
  
  #set up matrices
  f1 <- matrix(nrow = nrow(simsNoise), ncol = length(noiseScaled))
  f2 <- matrix(nrow = nrow(simsNoise), ncol = length(noiseScaled))
  
  # multiply intercept by 1; slope by value of the covariate
  for(i in 1:nrow(simsNoise)){
    for(j in 1:length(noiseScaled)){
      f1[i, j] <- sum(simsNoise[i, ] * c(1, 0, noiseScaled[j])) #non-pair
      f2[i, j] <- sum(simsNoise[i, ] * c(1, 1, noiseScaled[j])) #pair
    }
  }
  
  #back-transform from logit
  f1bt <- plogis(f1)
  f2bt <- plogis(f2)
  
  #convert the scaled covariate values to their real scale
  noiseBT <- (noiseScaled * zsd_noise) + zmean_noise
  
  #set up dataframe
  noisePlotinference <- data.frame(x = rep(noiseBT, 2),
                                 y = c(apply(f1bt, 2, mean), apply(f2bt, 2, mean)),
                                 lo = c(apply(f1bt, 2, quantile, probs = 0.025), 
                                        apply(f2bt, 2, quantile, probs = 0.025)),
                                 hi = c(apply(f1bt, 2, quantile, probs = 0.975), 
                                        apply(f2bt, 2, quantile, probs = 0.975)),
                                 grp = c(rep('Non-pairs', length(noiseBT)), rep('Pairs', length(noiseBT))))
  
  
#EFFORT ####  

  #import raw covariate values and calculate mean and SD
  effortRaw <- fread('data2018/covariates/raw/effort_raw_2018.csv', header = TRUE)
  effortRaw <- data.frame('effort' = unlist(effortRaw[,-1]))
    (zmean_effort <- mean(effortRaw$effort, na.rm = TRUE))
    (zsd_effort   <- sd(effortRaw$effort, na.rm = TRUE))
  
  #import the scaled covariate and convert to vector
  effortStd <- fread('data2018/covariates/effort_2018.csv')
  effortStd <- effortStd[,-c(1:2)]
    effortStd <- unlist(effortStd, use.names = FALSE)
    effortStd <- effortStd[order(effortStd)]
  
  #create a vector of values within the range of the real (scaled) values  
  effortScaled <- seq(min(effortStd), max(effortStd), length.out = 100)
  
  #extract relevant columns from the posterior draws (intercept/s and slope for this covariate)
  simsEffort <- inference_model$sims.list$gamma[,c(1:2,4)]
  
  #set up matrices
  g1 <- matrix(nrow = nrow(simsEffort), ncol = length(effortScaled))
  g2 <- matrix(nrow = nrow(simsEffort), ncol = length(effortScaled))
  
  #multiply intercept by 1; slope by value of the covariate
  for(i in 1:nrow(simsEffort)){
    for(j in 1:length(effortScaled)){
      g1[i, j] <- sum(simsEffort[i, ] * c(1, 0, effortScaled[j])) #non-pair
      g2[i, j] <- sum(simsEffort[i, ] * c(1, 1, effortScaled[j])) #pair
    }
  }
  
  #back-transform from logit
  g1bt <- plogis(g1)
  g2bt <- plogis(g2)
  
  #convert the scaled covariate values to their real scale
  effortBT <- (effortScaled * zsd_effort) + zmean_effort
  
  #set up dataframe
  effortPlotinference <- data.frame(x = rep(effortBT, 2),
                                  y = c(apply(g1bt, 2, mean), apply(g2bt, 2, mean)),
                                  lo = c(apply(g1bt, 2, quantile, probs = 0.025), 
                                         apply(g2bt, 2, quantile, probs = 0.025)),
                                  hi = c(apply(g1bt, 2, quantile, probs = 0.975), 
                                         apply(g2bt, 2, quantile, probs = 0.975)),
                                  grp = c(rep('Non-pairs', length(effortBT)), rep('Pairs', length(effortBT))))
  
  
## -----------------------------------------------------------------------------
## 4) Create marginal plots (detection) ####
  
## NOISE
  nnInf <- ggplot(noisePlotinference, aes(x, y)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = grp), alpha = 0.3) +
    geom_line(size = 1.5, aes(color = grp, linetype = grp)) +
    ylab('Weekly detection probability (p) \u00B1 95% CI') + xlab('Noise (dBFS)') +
    #scale_x_continuous(breaks = seq(0, 8000, 2000), limits = c(0, 5000)) +
    ylim(c(0,1)) +
    geom_rug(data = noiseRaw, mapping = aes(x = noise), inherit.aes = FALSE) +
    scale_fill_manual(values = c('Non-pairs' = 'black', 'Pairs' = 'darkblue')) +
    scale_colour_manual(values = c('Non-pairs' = 'black', 'Pairs' = 'darkblue')) +
    # scale_linetype_manual(values = c('Non-pairs' = 'solid', 'Pairs' = 'twodash')) +
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
  nnInf
  
  
  ## Predict a couple values for reporting:
  
    #p at min noise
    noisePlotinference$y[noisePlotinference$x == min(noisePlotinference$x) & 
                         noisePlotinference$grp %in% 'Non-pairs']         #0.13 for non-pairs
    noisePlotinference$y[noisePlotinference$x == min(noisePlotinference$x) & 
                         noisePlotinference$grp %in% 'Pairs']             #0.45 for non-pairs
    
    #p at mean noise
    noiseMeanIndex <- noisePlotinference$x[which.min(abs(noisePlotinference$x- zmean_noise))] #find x value closest to mean
    noisePlotinference$y[noisePlotinference$x == noiseMeanIndex & 
                         noisePlotinference$grp %in% 'Non-pairs']         #0.03 for non-pairs
    noisePlotinference$y[noisePlotinference$x == noiseMeanIndex & 
                         noisePlotinference$grp %in% 'Pairs']             #0.16 for non-pairs
    
  
## EFFORT
  eeInf <- ggplot(effortPlotinference, aes(x, y)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = grp), alpha = 0.3) +
    geom_line(size = 1.5, aes(color = grp, linetype = grp)) +
    ylab('Weekly detection probability (p) \u00B1 95% CI') + 
    xlab('Effort (recording minutes)') +
    #scale_x_continuous(breaks = seq(0, 8000, 2000), limits = c(0, 5000)) +
    ylim(c(0,1)) +
    geom_rug(data = effortRaw, mapping = aes(x = effort), inherit.aes = FALSE) +
    scale_fill_manual(values = c('Non-pairs' = 'black', 'Pairs' = 'darkblue')) + 
    scale_colour_manual(values = c('Non-pairs' = 'black', 'Pairs' = 'darkblue')) +
    #scale_linetype_manual(values = c('Non-pairs' = 'solid', 'Pairs' = 'twodash')) +
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
  eeInf
  
  ## Predict a couple values for reporting (intended 3360 min):
  effortPlotinference$y[effortPlotinference$x == 3360 &
                        effortPlotinference$grp %in% 'Non-pairs']  #0.04 for non-pairs
  effortPlotinference$y[effortPlotinference$x == 3360 &
                        effortPlotinference$grp %in% 'Pairs']      #0.19 for pairs
  
  
## Combine and export
  Fig3 <- ggarrange(nnInf + rremove("ylab"), eeInf + rremove("ylab"), 
                    labels = c("a", "b"), font.label = list(size = 24), 
                    vjust = c(1.5,1.5), hjust = c(-4,-3.5), 
                    ncol = 1, nrow = 2, 
                    common.legend = TRUE, legend = 'right') +
    theme(plot.margin = margin(0.1,1,0.1,1, "cm"))
  
  Fig3_shared_axis <- annotate_figure(Fig3, 
                                      left = textGrob('Weekly detection probability (p) \u00B1 95% CI', 
                                                      rot = 90, vjust = 2, hjust = 0.5,
                                                      gp = gpar(fontsize = 18)))
  
  ## FOR FIGURE 3
  tiff(filename = 'figures/fig3new.tif', height = 5600, width = 5200, units = 'px',
       res = '800', compression = 'lzw')
  print(Fig3_shared_axis)
  dev.off()  
  

  
## -----------------------------------------------------------------------------  
## 5) Explore covariate effects on occupancy (psi) ####
  
## Plot posterior means and CI (50% = thick lines) (95% = thin lines):
  MCMCplot(inference_model, params = 'beta', main = '"psi" parameters', ref_ovl = TRUE)
  #[1]: intercept, [2]: NR500, [3]: study area, [4]: BO_total
  
## Calculate means, SDs, and odds ratios:
  
  #intercepts
  psi_intMean_inf  = mean(inference_model$sims.list$beta[,1])
  psi_intSD_inf    = sd(inference_model$sims.list$beta[,1])
  psi_intQuant_inf = quantile(inference_model$sims.list$beta[,1], probs = c(0.025, 0.975))
  
  #forest suitability mean 500m
  psi_forestSlope_inf = mean(inference_model$sims.list$beta[,2])
  psi_forestSD_inf    = sd(inference_model$sims.list$beta[,2])
  psi_forestQuant_inf = quantile(inference_model$sims.list$beta[,2], probs = c(0.025, 0.975))
  psi_forestOR_inf    = median(exp(inference_model$sims.list$beta[,2]))       #odds ratio - use median
  psi_forestOR_LCI_inf = quantile(exp(inference_model$sims.list$beta[,2]), probs = c(0.025)) #odds ratio LCI
  psi_forestOR_UCI_inf = quantile(exp(inference_model$sims.list$beta[,2]), probs = c(0.975)) #odds ratio UCI
  
  #study area
  psi_areaSlope_inf = mean(inference_model$sims.list$beta[,3])
  psi_areaSD_inf    = sd(inference_model$sims.list$beta[,3])
  psi_areaQuant_inf = quantile(inference_model$sims.list$beta[,3], probs = c(0.025, 0.975))
  psi_areaOR_inf    = median(exp(inference_model$sims.list$beta[,3]))         #odds ratio - use median
  psi_areaOR_LCI_inf = quantile(exp(inference_model$sims.list$beta[,3]), probs = c(0.025))  #odds ratio LCI
  psi_areaOR_UCI_inf = quantile(exp(inference_model$sims.list$beta[,3]), probs = c(0.975))  #odds ratio UCI
  
  #total barred owl
  psi_BOSlope_inf = mean(inference_model$sims.list$beta[,4])
  psi_BOSD_inf    = sd(inference_model$sims.list$beta[,4])
  psi_BOQuant_inf = quantile(inference_model$sims.list$beta[,4], probs = c(0.025, 0.975))
  psi_BOOR_inf    = median(exp(inference_model$sims.list$beta[,4]))           #odds ratio - use median
  psi_BOOR_LCI_inf = quantile(exp(inference_model$sims.list$beta[,4]), probs = c(0.025)) #odds ratio LCI
  psi_BOOR_UCI_inf = quantile(exp(inference_model$sims.list$beta[,4]), probs = c(0.975)) #odds ratio UCI
  
  psi_params_inference <- data.frame('parameter' = rep('psi', 4),
                                  'covariate' = c('Intercept','Forest suitability mean 200m',
                                                  'Study area', 'Total barred owl calling'),
                                  'mean' = c(psi_intMean_inf, psi_forestSlope_inf, psi_areaSlope_inf,
                                             psi_BOSlope_inf),
                                  'SD' = c(psi_intSD_inf, psi_forestSD_inf, psi_areaSD_inf, psi_BOSD_inf),
                                  'LCI_95' = c(psi_intQuant_inf[[1]], psi_forestQuant_inf[[1]],
                                               psi_areaQuant_inf[[1]], psi_BOQuant_inf[[1]]),
                                  'UCI_95' = c(psi_intQuant_inf[[2]], psi_forestQuant_inf[[2]],
                                               psi_areaQuant_inf[[2]], psi_BOQuant_inf[[2]]),
                                  'F_score' = inference_model$f$beta,
                                  'Odds_ratio' = c(NA, psi_forestOR_inf, psi_areaOR_inf, psi_BOOR_inf),
                                  'OR_LCI' = c(NA, psi_forestOR_LCI_inf, psi_areaOR_LCI_inf, psi_BOOR_LCI_inf),
                                  'OR_UCI' = c(NA, psi_forestOR_UCI_inf, psi_areaOR_UCI_inf, psi_BOOR_UCI_inf))
  ## FOR SI TABLES:  
  psi_params_inference %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')
  
  
## -----------------------------------------------------------------------------  
## 6) Explore covariate effects on conditional pair occupancy (R) ####
  
## Plot posterior means and CI (50% = thick lines) (95% = thin lines):
  MCMCplot(inference_model, params = 'beta2', main = '"R" parameters', ref_ovl = TRUE)
  #[1]: intercept, [2]: NR500, [3]: study area, [4]: BO_total
  
## Calculate means, SDs, and odds ratios:
  
  #intercepts
  r_intMean_inf  = mean(inference_model$sims.list$beta2[,1])
  r_intSD_inf    = sd(inference_model$sims.list$beta2[,1])
  r_intQuant_inf = quantile(inference_model$sims.list$beta2[,1], probs = c(0.025, 0.975))
  
  #forest suitability weighted mean 200m
  r_forestSlope_inf = mean(inference_model$sims.list$beta2[,2])
  r_forestSD_inf    = sd(inference_model$sims.list$beta2[,2])
  r_forestQuant_inf = quantile(inference_model$sims.list$beta2[,2], probs = c(0.025, 0.975))
  r_forestOR_inf    = median(exp(inference_model$sims.list$beta2[,2]))        #odds ratio - use median
  r_forestOR_LCI_inf = quantile(exp(inference_model$sims.list$beta2[,2]), probs = c(0.025)) #odds ratio LCI
  r_forestOR_UCI_inf = quantile(exp(inference_model$sims.list$beta2[,2]), probs = c(0.975)) #odds ratio UCI
  
  #study area
  r_areaSlope_inf = mean(inference_model$sims.list$beta2[,3])
  r_areaSD_inf    = sd(inference_model$sims.list$beta2[,3])
  r_areaQuant_inf = quantile(inference_model$sims.list$beta2[,3], probs = c(0.025, 0.975))
  r_areaOR_inf    = median(exp(inference_model$sims.list$beta2[,3]))          #odds ratio - use median
  r_areaOR_LCI_inf = quantile(exp(inference_model$sims.list$beta2[,3]), probs = c(0.025)) #odds ratio LCI
  r_areaOR_UCI_inf = quantile(exp(inference_model$sims.list$beta2[,3]), probs = c(0.975)) #odds ratio UCI
  
  #total barred owl
  r_BOSlope_inf = mean(inference_model$sims.list$beta2[,4])
  r_BOSD_inf    = sd(inference_model$sims.list$beta2[,4])
  r_BOQuant_inf = quantile(inference_model$sims.list$beta2[,4], probs = c(0.025, 0.975))
  r_BOOR_inf    = median(exp(inference_model$sims.list$beta2[,4]))            #odds ratio - use median
  r_BOOR_LCI_inf = quantile(exp(inference_model$sims.list$beta2[,4]), probs = c(0.025)) #odds ratio LCI
  r_BOOR_UCI_inf = quantile(exp(inference_model$sims.list$beta2[,4]), probs = c(0.975))     #odds ratio UCI
  
  r_params_inference <- data.frame('parameter' = rep('R', 4),
                                'covariate' = c('Intercept','Forest suitability mean 200m','Study area',
                                                'Total barred owl calling'),
                                'mean' = c(r_intMean_inf, r_forestSlope_inf, r_areaSlope_inf,
                                           r_BOSlope_inf),
                                'SD' = c(r_intSD_inf, r_forestSD_inf, r_areaSD_inf, r_BOSD_inf),
                                'LCI_95' = c(r_intQuant_inf[[1]], r_forestQuant_inf[[1]],
                                             r_areaQuant_inf[[1]], r_BOQuant_inf[[1]]),
                                'UCI_95' = c(r_intQuant_inf[[2]], r_forestQuant_inf[[2]],
                                             r_areaQuant_inf[[2]], r_BOQuant_inf[[2]]),
                                'F_score' = inference_model$f$beta2,
                                'Odds_ratio' = c(NA, r_forestOR_inf, r_areaOR_inf, r_BOOR_inf),
                                'OR_LCI' = c(NA, exp(r_forestQuant_inf[[1]]),
                                             exp(r_areaQuant_inf[[1]]), exp(r_BOQuant_inf[[1]])),
                                'OR_UCI' = c(NA, exp(r_forestQuant_inf[[2]]),
                                             exp(r_areaQuant_inf[[2]]), exp(r_BOQuant_inf[[2]])))
  ## FOR SI TABLES:  
  r_params_inference %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')
  

## -----------------------------------------------------------------------------
## 7) Calculate predicted/fitted values (occupancy) ####
  
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
  simsNRpsi <- inference_model$sims.list$beta[,c(1,2,3)]  #3 is the col for study area
  simsNRr   <- inference_model$sims.list$beta2[,c(1,2,3)] #3 is the col for study area
  
  #set up matrices
  l1 <- matrix(nrow = nrow(simsNRpsi), ncol = length(nr500scaled))
  l2 <- matrix(nrow = nrow(simsNRr), ncol = length(nr500scaled))
  l3 <- matrix(nrow = nrow(simsNRpsi), ncol = length(nr500scaled))
  l4 <- matrix(nrow = nrow(simsNRr), ncol = length(nr500scaled))
  
  #multiply intercept by 1; slope by value of the covariate
  for(i in 1:nrow(simsNRpsi)){
    for(j in 1:length(nr500scaled)){
      l1[i, j] <- sum(simsNRpsi[i, ] * c(1, nr500scaled[j], 0)) #psi for OLY (0)
      l3[i, j] <- sum(simsNRpsi[i, ] * c(1, nr500scaled[j], 1)) #psi for COA (1)   
    }
  }
  for(i in 1:nrow(simsNRr)){
    for(j in 1:length(nr500scaled)){
      l2[i, j] <- sum(simsNRr[i, ] * c(1, nr500scaled[j], 0))  #R for OLY (0)
      l4[i, j] <- sum(simsNRr[i, ] * c(1, nr500scaled[j], 1))  #R for COA (1)
    }
  }
  
  #back-transform from logit
  occBTnr_OLY <- plogis(l1)
  rBTsuit_OLY   <- plogis(l2)
  
  occBTnr_COA <-plogis(l3)
  rBTsuit_COA <- plogis(l4)
  
  #calculate pair occupancy, etc.
  pairBTnr_OLY  <- occBTnr_OLY * rBTsuit_OLY      #prob occ by pairs (psi*R)
  unoccBTnr_OLY <- 1 - occBTnr_OLY                #prob unoccupied (1 - psi)
  #occNoPairBTsuit <- occBTsuit * (1 - rBTsuit)        #prob occ by non-pair (psi * (1-R))
  
  pairBTnr_COA  <- occBTnr_COA * rBTsuit_COA
  unoccBTnr_COA <- 1 - occBTnr_COA
  
  #convert the scaled covariate values that I used to their real scale
  nrBT <- (nr500scaled * zsd_nr500) + zmean_nr500
  
  #set up dataframe
  NRplotinference <- data.frame(x = rep(nrBT,4),
                              y = c(apply(occBTnr_OLY, 2, mean), apply(pairBTnr_OLY, 2, mean),
                                    apply(occBTnr_COA, 2, mean), apply(pairBTnr_COA, 2, mean)),
                              lo = c(apply(occBTnr_OLY, 2, quantile, probs = 0.025), 
                                    apply(pairBTnr_OLY, 2, quantile, probs = 0.025),
                                    apply(occBTnr_COA, 2, quantile, probs = 0.025), 
                                    apply(pairBTnr_COA, 2, quantile, probs = 0.025)),
                              hi = c(apply(occBTnr_OLY, 2, quantile, probs = 0.975), 
                                    apply(pairBTnr_OLY, 2, quantile, probs = 0.975),
                                    apply(occBTnr_COA, 2, quantile, probs = 0.975),
                                    apply(pairBTnr_COA, 2, quantile, probs = 0.975)),
                              grp = c(rep('Any owl', length(nrBT)), rep('Pair', length(nrBT)),
                                     rep('Any owl', length(nrBT)), rep('Pair', length(nrBT))),
                              study_area = c(rep('OLY', length(nrBT)*2), rep('COA', length(nrBT)*2)))

  
## STUDY AREA ####
  
  #factor so just need both values (0 and 1)
  areaScaled <- c(0,1)
  
  #extract relevant columns from the posterior draws (intercept/s and slope for this covariate)
  simsAreaPsi <- inference_model$sims.list$beta[,c(1,3)]  
  simsAreaR <- inference_model$sims.list$beta2[,c(1,3)]  
  
  #set up matrices
  m1 <- matrix(nrow = nrow(simsAreaPsi), ncol = length(areaScaled))
  m2 <- matrix(nrow = nrow(simsAreaR), ncol = length(areaScaled))
  
  #multiply intercept by 1; slope by value of the covariate
  for(i in 1:nrow(simsAreaPsi)){
    for(j in 1:length(areaScaled)){
      m1[i, j] <- sum(simsAreaPsi[i, ] * c(1, areaScaled[j]))
    }
  }
  for(i in 1:nrow(simsAreaR)){
    for(j in 1:length(areaScaled)){
      m2[i, j] <- sum(simsAreaR[i, ] * c(1, areaScaled[j]))
    }
  }
  
  #back-transform from logit
  occBTarea <- plogis(m1)
  rBTarea   <- plogis(m2)
  
  #calculate pair occupancy etc.
  pairBTarea  <- occBTarea * rBTarea             #prob occ by pairs (psi*R)
  #unoccBTarea <- 1 - occBTarea                  #prob unoccupied (1 - psi)
  occNoPairBTarea <- occBTarea * (1 - rBTarea)  #prob occ by non-pair (psi * (1-R))
  
  #set up dataframe
  areaPlotInf <- data.frame(x = rep(areaScaled,3),
                           y = c(apply(occBTarea, 2, mean), apply(pairBTarea, 2, mean), 
                                 apply(occNoPairBTarea, 2, mean)),
                           lo = c(apply(occBTarea, 2, quantile, probs = 0.025), 
                                  apply(pairBTarea, 2, quantile, probs = 0.025),
                                  apply(occNoPairBTarea, 2, quantile, probs = 0.025)),
                           hi = c(apply(occBTarea, 2, quantile, probs = 0.975), 
                                  apply(pairBTarea, 2, quantile, probs = 0.975), 
                                  apply(occNoPairBTarea, 2, quantile, probs = 0.975)),
                           grp = c(rep('Any owl', length(areaScaled)), rep('Pair', length(areaScaled)),
                                   rep('Non-pair', length(areaScaled))))
  

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
  simsBOtotalPsi <- inference_model$sims.list$beta[,c(1,3,4)]  #3 is the column for study area
  simsBOtotalR   <- inference_model$sims.list$beta2[,c(1,3,4)] #3 is the column for study area
  
  #set up matrices
  n1 <- matrix(nrow = nrow(simsBOtotalPsi), ncol = length(boTotalScaled))
  n2 <- matrix(nrow = nrow(simsBOtotalR), ncol = length(boTotalScaled))
  n3 <- matrix(nrow = nrow(simsBOtotalPsi), ncol = length(boTotalScaled))
  n4 <- matrix(nrow = nrow(simsBOtotalR), ncol = length(boTotalScaled))
  
  #multiply intercept by 1; slope by value of the covariate
  for(i in 1:nrow(simsBOtotalPsi)){
    for(j in 1:length(boTotalScaled)){
      n1[i, j] <- sum(simsBOtotalPsi[i, ] * c(1, 0, boTotalScaled[j])) #psi for OLY 
      n3[i, j] <- sum(simsBOtotalPsi[i, ] * c(1, 1, boTotalScaled[j])) #psi for COA
    }
  }
  for(i in 1:nrow(simsBOtotalR)){
    for(j in 1:length(boTotalScaled)){
      n2[i, j] <- sum(simsBOtotalR[i, ] * c(1, 0, boTotalScaled[j])) #R for OLY
      n4[i, j] <- sum(simsBOtotalR[i, ] * c(1, 1, boTotalScaled[j])) #R for COA
    }
  }
  
  #back-transform from logit
  occBTbo_OLY <- plogis(n1)
  rBTbo_OLY   <- plogis(n2)
  occBTbo_COA <- plogis(n3)
  rBTbo_COA   <- plogis(n4)
  
  #calculate pair occupancy etc.
  pairBTbo_OLY  <- occBTbo_OLY * rBTbo_OLY      #prob occ by pairs (psi*R)
  # unoccBTbo_OLY <- 1 - occBTbo_OLY                #prob unoccupied (1 - psi)
  # occNoPairBTbo <- occBTbo * (1 - rBTbo)   #prob occ by non-pair (psi * (1-R))
  pairBTbo_COA  <- occBTbo_COA * rBTbo_COA
  # unoccBTbo_COA <- 1 - occBTbo_COA  
  
  #convert the scaled covariate values that I used to their real scale
  boTotalBT <- (boTotalScaled * zsd_boTotal) + zmean_boTotal
  
  #set up dataframe
  boTotalPlotinference <- data.frame(x = rep(boTotalBT,2),
                                    y = c(apply(occBTbo_OLY, 2, mean), apply(pairBTbo_OLY, 2, mean),
                                          apply(occBTbo_COA, 2, mean), apply(pairBTbo_COA, 2, mean)),
                                    lo = c(apply(occBTbo_OLY, 2, quantile, probs = 0.025), 
                                           apply(pairBTbo_OLY, 2, quantile, probs = 0.025),
                                           apply(occBTbo_COA, 2, quantile, probs = 0.025),
                                           apply(pairBTbo_COA, 2, quantile, probs = 0.025)),
                                    hi = c(apply(occBTbo_OLY, 2, quantile, probs = 0.975), 
                                           apply(pairBTbo_OLY, 2, quantile, probs = 0.975),
                                           apply(occBTbo_COA, 2, quantile, probs = 0.975),
                                           apply(pairBTbo_COA, 2, quantile, probs = 0.975)),
                                    grp = c(rep('Any NSO', length(boTotalBT)), rep('Pair', length(boTotalBT)),
                                            rep('Any NSO', length(boTotalBT)), rep('Pair', length(boTotalBT))),
                                    study_area = c(rep('OLY', length(boTotalScaled)*2), rep('COA', length(boTotalScaled)*2)))
  
  
## -----------------------------------------------------------------------------
## 8. Create marginal plots (occupancy) ####
  
## NR 500
  rr <- ggplot(NRplotinference, aes(x/10000, y)) +     #if index mean: x/10000
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = grp), alpha = 0.3) +
    geom_line(size = 1.5, aes(color = grp, linetype = grp)) +
    facet_grid(~study_area) +
    ylab('Probability of occupancy \u00B1 95 CI') + xlab('Mean NR forest suitability index (500m)') +
    # ylab('Probability of occupancy \u00B1 95 CI') + xlab('NR forest suitability (weighted mean)') +
    #scale_x_continuous(breaks = seq(0, 8000, 2000), limits = c(0, 5000)) +
    ylim(c(0,1)) +
    scale_fill_manual(values = c('Any owl' = 'black', 'Pair' = 'darkblue')) + 
    scale_colour_manual(values = c('Any owl' = 'black', 'Pair' = 'darkblue')) +
    scale_linetype_manual(values = c('Any owl' = 'solid', 'Pair' = 'twodash')) +
    geom_rug(data = nr500Raw, mapping = aes(x = MEAN/10000), inherit.aes = FALSE) +   #if index mean
    # geom_rug(data = suitRaw, mapping = aes(x=SUM_AREA_WT_INDEX), inherit.aes=FALSE) + #if wt  mean
    # geom_rug(data = suitRaw, mapping = aes(x = HABa), inherit.aes = FALSE) +  #if HABa
    theme(panel.background = element_rect(fill = 'transparent'),
          axis.line = element_line(),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 16),
          legend.title = element_blank(),
          #legend.position = 'none',
          legend.background = element_rect(fill='transparent')
    )   
  rr
  #don't include this figure; use one from "no_area" model instead
  
  
## STUDY AREA 
  
  #add study area names to dataframe
  areaPlotInf$x <- ifelse(areaPlotInf$x == 0, 'OLY', 'COA')
  areaPlotInf$grp2 <- ifelse(areaPlotInf$grp %in% 'Pair', 'Pair\noccupancy', 
                            ifelse(areaPlotInf$grp %in% 'Any owl', 'Use', 'Unoccupied'))
  areaPlotInf <- areaPlotInf[areaPlotInf$grp2 %in% c('Use','Pair\noccupancy'),]

  aa <- ggplot(areaPlotInf, aes(x, y, linetype = grp2)) +
    geom_pointrange(aes(ymin = lo, ymax = hi, color = grp2, shape = grp2), 
                    size = 1, position = position_dodge(width = 0.3)) +
    ylab('Probability \u00B1 95 CI') + xlab('Study Area') +
    ylim(c(0,0.5)) +
    scale_shape_manual(values = c('Use' = 16, 'Pair\noccupancy' = 17)) +
    scale_colour_manual(values = c('Use' = 'black', 'Pair\noccupancy' = 'darkblue')) +
    scale_linetype_manual(values = c('Use' = 'solid', 'Pair\noccupancy' = 'dashed')) +
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
  aa

## Format margins to match other figures
  Fig5 <- ggarrange(aa + rremove("ylab"), vjust = c(1.5,1.5), hjust = c(-4,-3.5), 
                    ncol = 1, nrow = 1, 
                    common.legend = TRUE, legend = 'right') +
    theme(plot.margin = margin(0.1,1,0.1,1, "cm"))
  
  Fig5_shared_axis <- annotate_figure(Fig5,
                                      left = textGrob('Probability \u00B1 95 CI',
                                                      rot = 90, vjust = 1,
                                                      gp = gpar(cex = 1.6)))
  
  
  ## FOR FIGURE 5
  tiff(filename = 'figures/fig5.tif', height = 4000, width = 5200, units = 'px',
       res = '800', compression = 'lzw')
  print(Fig5_shared_axis)
  dev.off()  
  
  
## BARRED OWL TOTAL
  bb <- ggplot(boTotalPlotinference, aes(x, y, linetype = grp)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = grp), alpha = 0.3) +
    geom_line(size = 1.5, aes(color = grp)) +
    facet_grid(~study_area) +
    ylab('Probability of occupancy \u00B1 95 CI') + xlab('Barred owl detections (total)') +
    #scale_x_continuous(breaks = seq(0, 8000, 2000), limits = c(0, 5000)) +
    ylim(c(0,1)) +
    scale_fill_manual(values = c('Any NSO' = 'black', 'Pair' = 'darkblue')) + 
    scale_colour_manual(values = c('Any NSO' = 'black', 'Pair' = 'darkblue')) +
    geom_rug(data = boTotalRaw, mapping = aes(x = total), inherit.aes = FALSE) +
    theme(panel.background = element_rect(fill = 'transparent'),
          axis.line = element_line(),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 16),
          legend.title = element_blank(),
          #legend.position = 'none',
          legend.background = element_rect(fill='transparent')
    )   
  bb
  #don't include this figure; use one from "no_area" model instead
  
  
## -----------------------------------------------------------------------------
  
  