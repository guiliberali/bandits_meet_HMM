#######################################################################################################################
# Author: 
# Purpose:computes the mu & rho values, used for subsequent simulation
#          More specifically, puts out 'mu_2Sterminal.csv' and 'rho_2Sterminal.csv'

#  
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
#  - After that, set working directory to Study1/2. HMM estimates
#  - Before running this file, run 'stan_computing_posterior.R' in the same directory
#
#  Overview:
#     A) loads results from 'stan_computing_posterior.R'
#     B) computes mu parameters
#     C) computes rho parameters
#
#######################################################################################################################


#########################
# A) loads results from 'stan_computing_posterior.R'
#########################

PATH_RESULTS = 'Estimates for 2 states'

# load HMM M2/M1 results 
load("input_hmm_estimates/hmm2S_M2_terminal.RData")
load("input_hmm_estimates/hmm2S_M1_terminal.RData")

##  HMM  PARAMETERS DRAWS ----
stfit_m1=hmm2S_M1_terminal
stfit_m2=hmm2S_M2_terminal

#########################
# B) computes mu parameters
#########################

# Summarize 95% confidence interval for Mu parameters
mu_m1 <- summary(stfit_m1, pars = c("mu"), probs=c(0.025, 0.975), digits=3)$summary
mu_m2 <- summary(stfit_m2, pars = c("mu"), probs=c(0.025, 0.975), digits=3)$summary
mu_table <- cbind(mu_m1[, c(1, 4, 5)], mu_m2 [, c(1, 4, 5)])
stargazer(mu_table, summary = F, align=T)

# Take mean, write out MU result for simulations
mu_m1 <- colMeans(as.matrix(stfit_m1, pars = c("mu")))
mu_m2 <- colMeans(as.matrix(stfit_m2, pars = c("mu")))
mu=c(mu_m1, mu_m2)
mu_2Sterminal <- data.frame(cond = c('mu1mu11', 'mu2m11'), value = mu)

# write out estimated param
write.csv(mu_2Sterminal, file=paste0(PATH_RESULTS, "/mu_2Sterminal.csv"))

#########################
# B) computes rho parameters
#########################

# Summarize 95% confidence interval for Rho parameters
rho_m1 <- summary(stfit_m1, pars = c("rho"), probs=c(0.025, 0.975), digits=3)$summary
rho_m2 <- summary(stfit_m2, pars = c("rho"), probs=c(0.025, 0.975), digits=3)$summary
rho_table <- cbind(rho_m1[, c(1, 4, 5)],rho_m2[, c(1, 4, 5)])
stargazer(rho_table, summary = F, align=T)

# Take mean, write out Rho result for simulations
rho_m1 <- matrix(colMeans(as.matrix(stfit_m1, pars = c("rho"))), nrow=1)
rho_m2 <- matrix(colMeans(as.matrix(stfit_m2, pars = c("rho"))), nrow=1)
rho= c(rho_m1 , rho_m2)

# put out final results
rho_2Sterminal <- data.frame(cond =c('m1rho11', 'm1rho12', 'm1rho13', 'm1rho14', 'm2rho11', 'm2rho12', 'm2rho13', 'm2rho14'),
                             coef = rho)

write.csv(rho_2Sterminal, file=paste0(PATH_RESULTS, "/rho_2Sterminal.csv"))
