## multistate occupancy models with 2018 owl data

## evaluate most-supported model from submodel selection:
## p(NOISE + EFFORT) δ(.) ψ(.) R(.)

## CONTENTS:
## 1) Calculate derived parameters and compare with naive estimates
## 2) Explore covariate effects on detection (calculate odds-ratios, etc.)
## 3) Calculate predicted/fitted values
## 4) Create marginal plots
## 5)

library(dplyr)
library(magrittr)
library(knitr)
library(jagsUI)
library(MCMCvis)
library(kableExtra)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)


## -----------------------------------------------------------------------------
## 1) Calculate derived parameters and compare with naive estimates  ####
 
## Import model
final_model <- readRDS('results/04_model/model_output.rds')

  ## FOR TABLES 3 AND 4
  ## Model summary with parameter estimates and derived parameters tracked in model:
  final_model$summary %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')

## Import naive values
nso_dh <- fread('data2018/nso_dh_2018.csv')
  
  #quick exploration
  maxState <- apply(nso_dh[,c(3:42)], 1, max, na.rm = TRUE)
  (naive_psi <- length(maxState[maxState > 0]) / length(maxState))               # naive psi = 16.4%
  (naive_R <- length(maxState[maxState == 2]) / length(maxState[maxState > 0]))  # naive R = 17.6%
  (naive_occPair <- length(maxState[maxState == 2]) / length(maxState))
  (naive_occWithout <- length(maxState[maxState == 1]) / length(maxState))
  (naive_unOcc <- length(maxState[maxState == 0])/ length(maxState))
  
  naive_df <- data.frame('naive' = c('"psi" (34 out of 207 hexagons had NSO)', 
                                     '"R"   (6 of those 34 had an NSO pair)', 
                                     'occupied without pair (28 out of 270 hexagons)',
                                     'occupied with pair (6 out of 270 hexagons had an NSO pair)',
                                     'unoccupied (173 out of 207 hexagons)'), 
                         'estimate' = c(naive_psi, naive_R, naive_occWithout, 
                                        naive_occPair, naive_unOcc))
  
  ## TO REPORT IN TEXT
  naive_df %>% 
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')


## -----------------------------------------------------------------------------
## 2) Explore covariate effects on detection  ####

## Plot posterior means and CI (50% = thick lines) (95% = thin lines):
  MCMCplot(final_model, params = 'gamma', main = '"p" parameters', ref_ovl = TRUE)
  #[1]: intercept, [2]: intercept adjustment for pairs, [3]: noise, [4]: effort
  
## Calculate means, SDs, and odds ratios:
  
  #intercepts
  p_intMean  = mean(final_model$sims.list$gamma[,1])
  p_intSD    = sd(final_model$sims.list$gamma[,1])
  p_intQuant = quantile(final_model$sims.list$gamma[,1], probs = c(0.025, 0.975))
  
  p_intPairMean  = mean(final_model$sims.list$gamma[,2])
  p_intPairSD    = sd(final_model$sims.list$gamma[,2])
  p_intPairQuant = quantile(final_model$sims.list$gamma[,2], probs = c(0.025, 0.975))
  
  #noise
  p_noiseSlope = mean(final_model$sims.list$gamma[,3])
  p_noiseSD    = sd(final_model$sims.list$gamma[,3])
  p_noiseQuant = quantile(final_model$sims.list$gamma[,3], probs = c(0.025, 0.975))
  p_noiseOR    = median(exp(final_model$sims.list$gamma[,3])) #median
  p_noiseOR_LCI = quantile(exp(final_model$sims.list$gamma[,3]), probs = c(0.025))
  p_noiseOR_UCI = quantile(exp(final_model$sims.list$gamma[,3]), probs = c(0.975))
  
  #effort
  p_effortSlope = mean(final_model$sims.list$gamma[,4])
  p_effortSD    = sd(final_model$sims.list$gamma[,4])
  p_effortQuant = quantile(final_model$sims.list$gamma[,4], probs = c(0.025, 0.975))
  p_effortOR    = median(exp(final_model$sims.list$gamma[,4])) #median
  p_effortOR_LCI = quantile(exp(final_model$sims.list$gamma[,4]), probs = c(0.025))
  p_effortOR_UCI = quantile(exp(final_model$sims.list$gamma[,4]), probs = c(0.975))
  
  p_params <- data.frame('parameter' = rep('p', 4),
                         'covariate' = c('Intercept','Intercept for pairs','Background noise',
                                         'Recording effort'),
                         'mean' = c(p_intMean, p_intPairMean, p_noiseSlope, p_effortSlope),
                         'SD' = c(p_intSD, p_intPairSD, p_noiseSD, p_effortSD),
                         'LCI_95' = c(p_intQuant[[1]], p_intPairQuant[[1]], 
                                      p_noiseQuant[[1]], p_effortQuant[[1]]),
                         'UCI_95' = c(p_intQuant[[2]], p_intPairQuant[[2]], 
                                      p_noiseQuant[[2]], p_effortQuant[[2]]),
                         'F_score' = final_model$f$gamma,
                         'Odds_ratio' = c(NA, NA, p_noiseOR, p_effortOR),
                         'OR_LCI' = c(NA, NA, p_noiseOR_LCI, p_effortOR_LCI),
                         'OR_UCI' = c(NA, NA, p_noiseOR_UCI, p_effortOR_UCI))

  ## FOR TABLE 4  
  p_params %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')
  
  
## -----------------------------------------------------------------------------
## 3) Calculate predicted/fitted values ####
  
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
  simsNoise <- final_model$sims.list$gamma[,c(1:3)]
  
  #set up matrices
  f1 <- matrix(nrow = nrow(simsNoise), ncol = length(noiseScaled))
  f2 <- matrix(nrow = nrow(simsNoise), ncol = length(noiseScaled))
  
  #multiply intercept by 1; slope by value of the covariate
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
  noisePlot <- data.frame(x = rep(noiseBT, 4),
                          y = c(apply(f1bt, 2, mean), apply(f2bt, 2, mean)), 
                          lo = c(apply(f1bt, 2, quantile, probs = 0.025), 
                                 apply(f2bt, 2, quantile, probs = 0.025)),
                          hi = c(apply(f1bt, 2, quantile, probs = 0.975), 
                                 apply(f2bt, 2, quantile, probs = 0.975)),
                          grp = c(rep('Non-pairs', length(noiseBT)), rep('Pairs', length(noiseBT))),
                          season = c(rep('Early', length(noiseBT)), rep('Late', length(noiseBT))))
  
  
## EFFORT ####  

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
  simsEffort <- final_model$sims.list$gamma[,c(1:2,4)]
  
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
  effortPlot <- data.frame(x = rep(effortBT, 4),
                           y = c(apply(g1bt, 2, mean), apply(g2bt, 2, mean)),
                           lo = c(apply(g1bt, 2, quantile, probs = 0.025), 
                                  apply(g2bt, 2, quantile, probs = 0.025)),
                           hi = c(apply(g1bt, 2, quantile, probs = 0.975), 
                                  apply(g2bt, 2, quantile, probs = 0.975)),
                           grp = c(rep('Non-pairs', length(effortBT)), rep('Pairs', length(effortBT))),
                           season = c(rep('Early', length(noiseBT)), rep('Late', length(noiseBT))))
  
  
## -----------------------------------------------------------------------------
## 4) Create marginal plots ####
  
## NOISE
  nn <- ggplot(noisePlot, aes(x, y)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = grp), alpha = 0.3) +
    geom_line(size = 1.5, aes(color = grp, linetype = grp)) +
    ylab('Weekly detection probability (p) \u00B1 95% CI') + xlab('Noise (dBFS)') +
    #scale_x_continuous(breaks = seq(0, 8000, 2000), limits = c(0, 5000)) +
    ylim(c(0,1)) +
    geom_rug(data = noiseRaw, mapping = aes(x = noise), inherit.aes = FALSE) +
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
  nn
  
## EFFORT
  ee <- ggplot(effortPlot, aes(x, y)) +
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
  ee
  

## Combine and export
  combinedFig <- ggarrange(nn, ee, labels = c("a", "b"), ncol = 2, nrow = 1, 
                           common.legend = TRUE, legend = 'bottom', vjust = -1, hjust = -4.5, 
                           font.label = list(size = 24)) +
    theme(plot.margin = margin(2,0.1,0.1,0.1, "cm"))
  
  ## FIGURE 3
  # ggexport(combinedFig, width = 850, height = 450,
  #          filename = 'figures/marginals_noise_effort_model04.png')

  combinedFig2 <- ggarrange(nn + rremove("ylab"), ee + rremove("ylab"), labels = c("a", "b"), ncol = 1, nrow = 2, 
                           common.legend = TRUE, legend = 'bottom',  vjust = -1, hjust = -4.5, 
                           font.label = list(size = 24)) +
    theme(plot.margin = margin(2,0.1,0.1,0.1, "cm"))
  combinedFig2_shared_axis <- annotate_figure(combinedFig2, left = textGrob('Weekly detection probability (p) \u00B1 95% CI', rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
  
  tiff(filename = 'figures/test.tif', height = 5600, width = 5200, units = 'px', 
       res = '800', compression = 'lzw')
  print(combinedFig2_shared_axis)
  dev.off()  
  
  ## give it more room on the left 

  