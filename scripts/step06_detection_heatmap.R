## multistate occupancy models with 2018 owl data

## create detection probability heatmap from "posthoc" model
## p(NOISE + EFFORT) δ(.) ψ(.) R(.)

## CONTENTS:
## 1) Import model and extract detection probabilities
## 2) Extrapolate for 1-5 stations

library(data.table)
library(ggplot2)
library(reshape2)
library(MCMCvis)

## -----------------------------------------------------------------------------
## 1) Import model and extract detection probabilities ####

final_model <- readRDS('results/04_model/model_output.rds')
  #temporarily:
  # final_model <- readRDS('C:/Users/caral/Documents/_RESEARCH/Models/results2018_Jan2022/17_model_final_norm1priors/model_output.rds')

  (p1 <- mean(final_model$sims.list$pOcc))
  (p2 <- mean(final_model$sims.list$pPair))
  (delta <- mean(final_model$sims.list$condPair))
  (pPair <- mean(final_model$sims.list$truePair))


## -----------------------------------------------------------------------------
## 2) Extrapolate for 1-5 stations ####

#detecting any owl at a site where true state is paired
  n = 1:5
  (p_hexagon = 1 - (1 - p2)^n)
  plot(n, p_hexagon, type = "b", ylim = c(0,1), xlab = c('stations'), 
       main = 'hexagon-level p \n (site has a pair)')
  abline(h = 0.95, lty = 3)
  
  m = 1:8 #and for 1-8 weeks
  (p_week = 1 - (1 - p_hexagon[5])^m)
  plot(m, p_week, type = "b", ylim = c(0,1), xlab = c('weeks'), 
       main = 'weekly p for a hexagon with 5 stations \n (site has a pair)')
  abline(h = 0.95, lty = 3)

#detecting any owl at a site where true state is NOT paired
  o = 1:5
  (p_hexagon = 1 - (1 - p1)^o)
  plot(n, p_hexagon, type = "b", ylim = c(0,1), xlab = c('stations'), 
       main = 'hexagon-level p \n (site has a single)')
  abline(h = 0.95, lty = 3)
  
  q = 1:8 #and for 1-8 weeks
  (p_week = 1 - (1 - p_hexagon[5])^q)
  plot(q, p_week, type = "b", ylim = c(0,1), xlab = c('weeks'), 
       main = 'weekly p for a hexagon with 5 stations \n (site has a single)')
  abline(h = 0.95, lty = 3)

#Pr detect an owl AND PROPERLY IDENTIFY IT AS A PAIR at a site where there is a pair
  r = 1:5
  (p_hexagon = 1 - (1 - p2*delta)^r)
  plot(r, p_hexagon, type = "b", ylim = c(0,1), xlab = 'stations', 
       main = 'hexagon-level p \n (detect an owl and properly ID as pair)')
  abline(h = 0.95, lty = 3)
  
  s = 1:8 #and for 1-8 weeks
  (p_week = 1 - (1 - p_hexagon[5])^s)
  plot(s, p_week, type = 'b', ylim = c(0,1), xlab = 'weeks', 
       main = 'hexagon-level p \n (detect an owl and properly ID as pair)')  
  abline(h = 0.95, lty = 3)

  
## -----------------------------------------------------------------------------
## 2) Create heatmap ####
  
nStn = 1:5
nWk  = 1:8

#P(detect pair | site with pair)
  p_hex_truePair  <- 1 - (1 - p2*delta)^nStn
  truePair_df <- data.frame()
  for (ss in nStn){
    truePair_df <- rbind(truePair_df, 1 - (1- p_hex_truePair[ss])^nWk)
  }
  colnames(truePair_df) <- nWk
  truePair_df$stn <- nStn
  
  truePair_df_long <- melt(truePair_df, id.vars = 'stn')  #convert to long
  colnames(truePair_df_long) <- c('stn','wk','p')
  truePair_df_long$grp <- 'c'

#P(detect use | site with non-pair)
  p_hex_noPair  <- 1 - (1 - p1)^nStn
  nonPair_df <- data.frame()
  for (tt in nStn){
    nonPair_df <- rbind(nonPair_df, 1 - (1- p_hex_noPair[tt])^nWk)
  }
  colnames(nonPair_df) <- nWk
  nonPair_df$stn <- nStn
  
  nonPair_df_long <- melt(nonPair_df, id.vars = 'stn')    #convert to long
  colnames(nonPair_df_long) <- c('stn','wk','p')
  nonPair_df_long$grp <- 'b'

#P(detect use | site with pair)
  p_hex_Pair  <- 1 - (1 - p2)^nStn
  pair_df <- data.frame()
  for (uu in nStn){
    pair_df <- rbind(pair_df, 1 - (1- p_hex_Pair[uu])^nWk)
  }
  colnames(pair_df) <- nWk
  pair_df$stn <- nStn
  
  pair_df_long <- melt(pair_df, id.vars = 'stn')    #convert to long
  colnames(pair_df_long) <- c('stn','wk','p')
  pair_df_long$grp <- 'a'


#merge and plot
  comboDF <- rbind(truePair_df_long, nonPair_df_long)
  comboDF <- rbind(comboDF, pair_df_long)    

  hm <- ggplot(comboDF, aes(x = stn, wk)) +
    geom_tile(aes(fill = p)) +
    geom_text(aes(label = round(p, 2))) +
    facet_grid(~grp) +
    scale_fill_distiller(palette = 'Spectral') +        
    xlab('Stations') + ylab('Weeks') +
    theme(panel.background = element_rect(fill = 'transparent'),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5),
          # legend.text = element_text(size = 14),
          # legend.title = element_text(size = 16),
          # legend.background = element_rect(fill='transparent'),
          legend.position = 'none')
  hm
  #export 1000 x 450

#export
  tiff(filename = 'figures/s3_fig1.tif', height = 3500, width = 7000, units = 'px',
       res = '800', compression = 'lzw')
  print(hm)
  dev.off()  
  
  