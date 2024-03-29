## multistate occupancy models with 2018 owl data

## (02a) 'psi' submodel: p(NOISE + EFFORT) delta(.) psi(STUDY_AREA + NR_SUIT500 + BO_TOTAL) R(.)
## with initial 'normal1' priors: dnorm(0,0.368)

library(data.table)
library(jagsUI)

## clear R environment
ls()
rm(list = ls())
ls()

## -----------------------------------------------------------------------------
## set model name here

  modname = '02a_submodel_psi_initial'

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

#study area
  area <- fread('data2018/covariates/area_2018.csv')
  area_covar <- area$AREA
  area_covar[area_covar == 2] <- 0  #0 = OLY, 1 = COA
  
#barred owl calling (total)
  bo_total <- fread('data2018/covariates/bo_total_2018.csv')
  bo_total_covar <- bo_total$total_std
  
#forest nesting suitability (500m)
  nr500 <- fread('data2018/covariates/nr500_2018.csv')  
  nr500_mean <- as.numeric(nr500$MEAN_STD)


## -----------------------------------------------------------------------------
## import survey-level covariates (already standardized) ####

#background noise
  noise <- fread('data2018/covariates/noise_2018.csv')
  noise_covar <- unname(as.matrix(noise[,-c(1:2)]))
  
#effort
  effort <- fread('data2018/covariates/effort_2018.csv')
  effort_covar <- unname(as.matrix(effort[,-c(1:2)]))
  

## -----------------------------------------------------------------------------
## set up data input ####

jdata <- list(nhex = as.numeric(nrow(det_matrix)), nVisits = as.numeric(ncol(det_matrix)), 
              y = det_matrix,
              noise = noise_covar, effort = effort_covar,
              area = area_covar, nr500 = nr500_mean, bo_total = bo_total_covar)

str(jdata) #check dimensions (1x207 and 1x40)


## -----------------------------------------------------------------------------
## set up tracked parameters ####

params<-c('beta','beta2','gamma','iota')


## -----------------------------------------------------------------------------
## set initial values ####

#get max state of each site to use for initial values 
  maxState[maxState == 0] <- 2  #change 0s to 2s
  y2 <- maxState

#input number of covariates on each parameter:
  ncov_psi   = 3
  ncov_r     = 0
  ncov_p     = 2
  ncov_delta = 0

#set initial values:
inits<-function(){list(Occ.hex=y2,
                       beta=rep(0,ncov_psi+1), beta2=0, 
                       gamma=rep(0,ncov_p+2), iota=rep(0,ncov_delta+1))}
  

## -----------------------------------------------------------------------------
## run JAGS ####

  start <- Sys.time()
  output <- jagsUI(data = jdata, inits = inits, parameters = params, store.data = TRUE,
                   model = 'bugs/02a_submodel_psi_initial.bugs',
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
## save the mean and SD of initial run:
  
  saveRDS(output$mean$beta, file = paste(filedir, 'spike_gamma_mu.RDS', sep = '/'))
  saveRDS(output$sd$beta, file = paste(filedir, 'spike_gamma_SD.RDS', sep = '/'))

