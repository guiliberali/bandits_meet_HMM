
PATH_RESULTS = 'Estimates for 2 states'

# load HMM M2/M1 results 
load("hmm2S_M2_terminal.RData")
load("hmm2S_M1_terminal.RData")

##  HMM  PARAMETERS DRAWS ----
stfit_m1=hmm2S_M1_terminal
stfit_m2=hmm2S_M2_terminal

# Summarize 95% confidence interval for Mu parameters
mu_m1 <- summary(stfit_m1, pars = c("mu"), probs=c(0.025, 0.975), digits=3)$summary
mu_m2 <- summary(stfit_m2, pars = c("mu"), probs=c(0.025, 0.975), digits=3)$summary
stargazer(cbind(mu_m1[, c(1, 4, 5)], mu_m2 [, c(1, 4, 5)]), summary = F, align=T)

# Take mean, write out MU result for simulations
mu_m1 <- colMeans(as.matrix(stfit_m1, pars = c("mu")))
mu_m2 <- colMeans(as.matrix(stfit_m2, pars = c("mu")))
mu=c(mu_m1, mu_m2)
write.csv(mu, file=paste0(PATH_RESULTS, "mu_2Sterminal.csv"))

# Summarize 95% confidence interval for Rho parameters
rho_m1 <- summary(stfit_m1, pars = c("rho"), probs=c(0.025, 0.975), digits=3)$summary
rho_m2 <- summary(stfit_m2, pars = c("rho"), probs=c(0.025, 0.975), digits=3)$summary
stargazer(cbind(rho_m1[, c(1, 4, 5)],rho_m2[, c(1, 4, 5)]), summary = F, align=T)

# Take mean, write out Rho result for simulations
rho_m1 <- matrix(colMeans(as.matrix(stfit_m1, pars = c("rho"))), nrow=1)
rho_m2 <- matrix(colMeans(as.matrix(stfit_m2, pars = c("rho"))), nrow=1)
rho= c(rho_m1 , rho_m2)
write.csv(rho, file=paste0(PATH_RESULTS, "rho_2Sterminal.csv"))
