## compare parameters from model runs with different prior sets

library(MCMCvis)

## import model(s) -------------------------------------------------------------

model05a <- readRDS('results/05a_model_posthoc_norm1priors/model_output.rds')
  MCMCtrace(model05a, ind = TRUE)
  (model05a_means <- model05a$summary[,1])

model05b <- readRDS('results/05b_model_posthoc_norm2priors/model_output.rds')
  MCMCtrace(model05b, ind = TRUE)
  (model05b_means <- model05b$summary[,1])
  
model05c <- readRDS('results/05c_model_posthoc_unifpriors/model_output.rds')
  MCMCtrace(model05c, ind = TRUE)
  (model05c_means <- model05c$summary[,1])


## now compare them ------------------------------------------------------------

par(mfrow=c(1,3))

## 05a vs 05b
plot(mod05a_means[names(mod05a_means) %in% names(mod05b_means)], mod05b_means[names(mod05b_means) %in% names(mod05a_means)],
     col = 'lightblue', pch = 19, cex = 2, xlim = c(-2,3), ylim = c(-2,3), xlab = 'normal1', ylab = 'normal2') ; abline(0,1)
text(mod05a_means[names(mod05a_means) %in% names(mod05b_means)] ~ mod05b_means[names(mod05b_means) %in% names(mod05a_means)],
     labels = names(mod05b_means), cex = 0.9, font = 1, adj = 1)

## 05a vs 05c
plot(mod05a_means[names(mod05a_means) %in% names(mod05c_means)], mod05c_means[names(mod05c_means) %in% names(mod05a_means)],
     col = 'lightblue', pch = 19, cex = 2, xlim = c(-2,3), ylim = c(-2,3), xlab = 'normal1', ylab = 'uniform') ; abline(0,1)
text(mod05a_means[names(mod05a_means) %in% names(mod05c_means)] ~ mod05c_means[names(mod05c_means) %in% names(mod05a_means)],
     labels = names(mod05c_means), cex = 0.9, font = 1, adj = 1)

## 05b vs 05c
plot(mod05b_means[names(mod05b_means) %in% names(mod05c_means)], mod05c_means[names(mod05c_means) %in% names(mod05b_means)],
     col = 'lightblue', pch = 19, cex = 2, xlim = c(-2,3), ylim = c(-2,3), xlab = 'normal2', ylab = 'uniform') ; abline(0,1)
text(mod05b_means[names(mod05b_means) %in% names(mod05c_means)] ~ mod05c_means[names(mod05c_means) %in% names(mod05b_means)],
     labels = names(mod05c_means), cex = 0.9, font = 1, adj = 1)

