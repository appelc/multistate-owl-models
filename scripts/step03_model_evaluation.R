## multistate occupancy models with 2018 owl data

## evaluate most-supported model from submodel selection for prediction:
## p(NOISE + EFFORT) δ(.) ψ(.) R(.)

## CONTENTS:
## 1) Calculate derived parameters and compare with naive estimates
## 2) Explore covariate effects on detection (calculate odds-ratios, etc.)
## 3) Calculate psi_cond and psi_hat

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
library(grid)


## -----------------------------------------------------------------------------
## 1) Calculate derived parameters and compare with naive estimates  ####
 
## Import model
final_model <- readRDS('results/04_model/model_output.rds')
  #temporarily:
  final_model <- readRDS('C:/Users/caral/Documents/_RESEARCH/Models/results2018_Jan2022/17_model_final_norm1priors/model_output.rds')

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

  ## FOR TABLE 3   
  p_params %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')
  
  
  
## -----------------------------------------------------------------------------
## 3) Calculate psi_cond and psi_hat
  
## psi_cond: probability that a unit is occupied given owls NOT detected there
## (after MacKenzie et al. 2018 page 138)
  psiD = mean(final_model$sims.list$psiD)
  pOcc = mean(final_model$sims.list$pOcc) #use p1 for simplicity
  k = 40  #survey occasions
  
  ## TO REPORT IN-TEXT  
  (psiCond = (psiD * (1 - pOcc)^k) / (1 - psiD * (1 - (1 - pOcc)^k))) 

  
## psi_hat: proportion of sites occupied
## (after MacKenzie et al. 2018 section 6.1)
  s  = 207  #number of sites surveyed 
  sd = 34   #number of sites with detections
  
  ## TO REPORT IN-TEXT
  (psiHat = (sd + (s - sd)*psiCond) / s) 
  

  