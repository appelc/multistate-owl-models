## multistate occupancy models with 2018 owl data

## sequential-by-submodel selection 

library(dplyr)
library(magrittr)
library(knitr)
library(jagsUI)
library(MCMCvis)
library(kableExtra)
library(data.table)
library(ggplot2)
library(gridExtra)


## -----------------------------------------------------------------------------
## "p" SUBMODEL: p(NOISE + EFFORT + SEASON + BO_WEEKLY + NR_200m) δ(.) ψ(.) R(.) ####

p_submodel <- readRDS('results/01b_submodel_p_spikeSlab/model_output.rds')

## Look at indicator variable means/SDs:
  pCovariates <- c('int_adjustment','NOISE','EFFORT','SEASON','BO_WEEKLY','NR_200m')
  p_indicators <- NULL
  for(cc in 1:ncol(p_submodel$sims.list$w)){
    wMean <- mean(p_submodel$sims.list$w[,cc])
    wSD   <- sd(p_submodel$sims.list$w[,cc])
    wDF   <- data.frame('indicator' = paste('w[', cc, ']', sep = ''), 'covariate' = pCovariates[cc], 
                        'posterior_mean' = wMean, 'posterior_SD' = wSD)
    p_indicators <- rbind(p_indicators, wDF)
  }
  
  p_indicators %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')
  
## Plot parameter estimates with CI (50% = thick lines) (95% = thin lines):
  MCMCplot(p_submodel, params = 'gamma', ref_ovl = TRUE)  
  #[1]: intercept, [2]: intercept adjustment for pairs, [3]: NOISE, [4]: EFFORT, 
  #[5]: SEASON, [6]: BO_WEEKLY, [7]: NR_200m
  
## Summarize indicator variables across model iterations
  pmod <- p_submodel$sims.list$w
  pmod <- paste(pmod[,1],pmod[,2],pmod[,3],pmod[,4],pmod[,5],pmod[,6],sep="")
  pCombos <- sort(round(table(pmod)/length(pmod), 3), decreasing = TRUE)

  ## FOR TABLE 2  (order is: intercept adjustment for pairs, NOISE, EFFORT, SEASON, BO_WEEKLY, NR_200m)
  pCombos %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')
  
  
## -----------------------------------------------------------------------------
## no "delta" submodel  ####
  
## -----------------------------------------------------------------------------
## "psi" SUBMODEL: p(NOISE + EFFORT) δ(.) ψ(NR500 + STUDY_AREA + BO_TOTAL) R(.) ####
  
psi_submodel <- readRDS('results/02b_submodel_psi_spikeSlab/model_output.rds')

## Look at indicator variable means/SDs:
  psiCovariates <- c('NR_500m','AREA','BO_TOTAL')
  psi_indicators <- NULL
  for(dd in 1:ncol(psi_submodel$sims.list$w)){
    wMean <- mean(psi_submodel$sims.list$w[,dd])
    wSD   <- sd(psi_submodel$sims.list$w[,dd])
    wDF   <- data.frame('indicator' = paste('w[', dd, ']', sep = ''), 'covariate' = psiCovariates[dd],
                        'posterior_mean' = wMean, 'posterior_SD' = wSD)
    psi_indicators <- rbind(psi_indicators, wDF)
  }
  
  psi_indicators %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')

## Plot parameter estimates with CI (50% = thick lines) (95% = thin lines):
  MCMCplot(psi_submodel, params = 'beta', ref_ovl = TRUE)  
  #[1]: intercept, [2]: NR_500m, [3]: STUDY_AREA, [4]: BO_TOTAL

## Summarize indicator variables across model iterations
  psimod <- psi_submodel$sims.list$w
  psimod <- paste(psimod[,1],psimod[,2],psimod[,3],sep="")
  psiCombos <- sort(round(table(psimod)/length(psimod), 3),decreasing = TRUE)

  ## FOR TABLE 2 (order is: NR_500m, STUDY_AREA, BO_TOTAL)  
  psiCombos %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')


## --------------------------------------------------------------------------
## "R" SUBMODEL: p(NOISE + EFFORT) δ(.) ψ(.) R(NR500 + STUDY_AREA + BO_TOTAL) ####
  
r_submodel <- readRDS('results/03b_submodel_R_spikeSlab/model_output.rds')

## Look at indicator variable means/SDs:
  rCovariates <- c('NR_500m','AREA','BO_TOTAL')
  r_indicators <- NULL
  for(ee in 1:ncol(r_submodel$sims.list$w)){
    wMean <- mean(r_submodel$sims.list$w[,ee])
    wSD   <- sd(r_submodel$sims.list$w[,ee])
    wDF   <- data.frame('indicator' = paste('w[', ee, ']', sep = ''), 'covariate' = rCovariates[ee],
                        'posterior_mean' = wMean, 'posterior_SD' = wSD)
    r_indicators <- rbind(r_indicators, wDF)
  }
  
  r_indicators %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')

## Plot parameter estimates with CI (50% = thick lines) (95% = thin lines):
  MCMCplot(r_submodel, params = 'beta2', ref_ovl = TRUE)  
  #[1]: intercept, [2]: nr_500m, [3]: area, [4]: BO_total
  
## Summarize indicator variables across model iterations
  rmod <- r_submodel$sims.list$w
  rmod <- paste(rmod[,1],rmod[,2],rmod[,3],sep="")
  rCombos <- sort(round(table(rmod)/length(rmod), 3),decreasing = TRUE)
  
  ## FOR TABLE 2 (order is: NR_500m, STUDY_AREA, BO_TOTAL)
  rCombos %>%
    kbl(digits = 3) %>%
    kable_styling(bootstrap_options = 'striped', font_size = 14, full_width = FALSE, position = 'left')
  
  