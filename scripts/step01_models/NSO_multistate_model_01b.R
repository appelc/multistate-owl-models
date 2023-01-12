## multistate occupancy models with 2018 owl data

## (01b) 'p' submodel: p(NOISE + EFFORT + SEASON + BO_WEEKLY) delta(.) psi(.) R(.)
## with spike-and-slab priors

library(data.table)
library(jagsUI)

## clear R environment
ls()
rm(list = ls())
ls()

## -----------------------------------------------------------------------------
## set model name here

  modname = '01b_submodel_p_spikeSlab'

## -----------------------------------------------------------------------------
## create a folder to save outputs

  filedir = file.path('results//', modname)
  if (!dir.exists(filedir)) {dir.create(filedir, recursive = TRUE)
    print(paste('Created directory: ', modname, sep = ''))
  } else {stop('Directory already exists! Provide unique run_name')}


## -----------------------------------------------------------------------------
## import detection history ####
  nso_dh <- fread('data2018/nso_dh_2018.csv')
  nso_dh <- nso_dh[,-1]

#turn into matrix with 3s instead of 0s (better for JAGS)
  maxState <- apply(nso_dh[,c(2:41)], 1, max, na.rm = TRUE)
  nso_matrix <- unname(as.matrix(nso_dh[,-1]))
  det_matrix <- ifelse(nso_matrix == 0, 3, nso_matrix)


## -----------------------------------------------------------------------------
## import site-level covariates (already standardized) ####

#none for this submodel


## -----------------------------------------------------------------------------
## import survey-level covariates (already standardized) ####

#background noise
  noise <- fread('data2018/covariates/noise_2018.csv')
  noise_covar <- unname(as.matrix(noise[,-c(1:2)]))
  
#effort
  effort <- fread('data2018/covariates/effort_2018.csv')
  effort_covar <- unname(as.matrix(effort[,-c(1:2)]))
  
#season as a factor
  season <- fread('data2018/covariates/season_2018.csv')
  season_covar <- unname(as.matrix(season[,-c(1,2)]))

#barred owl calling (weekly)    
  bo_weekly <- fread('data2018/covariates/bo_weekly_2018.csv')
  bo_weekly_covar <- unname(as.matrix(bo_weekly[,-c(1:2)]))

#forest NR suitability (200m) -- station-scale
  nr200 <- fread('data2018/covariates/nr200_2018.csv')
  nr200_mean <- unname(as.matrix(nr200[,-c(1:2)]))


## -----------------------------------------------------------------------------
## set up data input ####

#load means/SDs from initial models to use for spike-and-slab priors
  spike_mu_gamma <- readRDS('results/01a_submodel_p_initial/spike_gamma_mu.RDS')
  spike_SD_gamma <- readRDS('results/01a_submodel_p_initial/spike_gamma_SD.RDS')
  
jdata <- list(nhex = as.numeric(nrow(det_matrix)), nVisits = as.numeric(ncol(det_matrix)), 
              y = det_matrix,
              noise = noise_covar, effort = effort_covar, season = season_covar, 
              bo_weekly = bo_weekly_covar, nr200 = nr200_mean,
              SpikeMuGamma = c(NA, spike_mu_gamma[2:length(spike_mu_gamma)]), #replace intercept with NA
              SpikeSDGamma = c(NA, spike_SD_gamma[2:length(spike_SD_gamma)])) #replace intercept with NA

str(jdata) #check dimensions (1x207 and 1x40)


## -----------------------------------------------------------------------------
## set up tracked parameters ####

params<-c('beta','beta2','gamma','iota','w')


## -----------------------------------------------------------------------------
## set initial values ####

#get max state of each site to use for initial values 
  maxState[maxState == 0] <- 2  #change 0s to 2s
  y2 <- maxState

#input number of covariates on each parameter
  ncov_psi   = 0
  ncov_r     = 0
  ncov_p     = 5
  ncov_delta = 0

#set initial values (don't provide inits for gamma):
inits<-function(){list(Occ.hex=y2,
                       beta=rep(0,ncov_psi+1), beta2=0, 
                       iota=rep(0,ncov_delta+1), w=rep(1,ncov_p+1))}

  
## -----------------------------------------------------------------------------
## run JAGS ####
  start <- Sys.time()
  output <- jagsUI(data = jdata, inits = inits, parameters = params, store.data = TRUE,
                   model = 'bugs/01b_submodel_p_spikeSlab.bugs',
                   n.chains=3, n.thin=1, n.iter=50000, n.burnin=10000, n.adapt=10000, parallel=TRUE)
  (end = Sys.time() - start)

  output
  #plot(output)

  
## -----------------------------------------------------------------------------
## save results

  saveRDS(output, file = paste(filedir, 'model_output.rds', sep = '/'))
  write.table(round(output$summary, 3), file = paste(filedir, 'model_summary.txt', sep = '/'), 
              sep = '\t', col.names = NA, row.names = TRUE)


## -----------------------------------------------------------------------------
## check indicator variable combinations from spike-and-slab run:

  mod <- output$sims.list$w
  mod <- paste(mod[,1],mod[,2],mod[,3],mod[,4],mod[,5],sep = '')
  table(mod)
  sort(round(table(mod)/length(mod), 3)) #top combo for "p"
  
  