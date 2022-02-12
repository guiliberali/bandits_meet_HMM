#######################################################################################################################
# Author: 
# Purpose: computes the bounce probability per state
#          More specifically, puts out 'psm_2Sterminal.csv' and 'pmf_bounce_user_2Sterminal.csv'
# 
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
#  - After that, set working directory to Study1/2. HMM estimates
#  - Before running this file, run 'stan_computing_posterior.R' in the same directory
#
#  Overview: 
#     A) loads data, and augment data to prepare for estimation 
#     B) estimate bounce probabilities 
#     C) Calculate the PSM
#
#######################################################################################################################


##########################
#  A) loads data, and augment data to prepare for estimation 
########################## 


# Get final dataset, as well as data per users from m1 and m2 
Data_final=read.csv(paste0('input_hmm_estimates',"/AB_data_05_June_2018.csv"), header = T)
user_sample_m1=read.csv(paste0('input_hmm_estimates',"/user_sample_m1.csv"), header=T)
user_sample_m2=read.csv(paste0('input_hmm_estimates',"/user_sample_m2.csv"), header=T)


# Gather posterior results
load("input_hmm_estimates/hmm2S_M2_terminal.RData")
load("input_hmm_estimates/hmm2S_M1_terminal.RData")
stfit_m1=hmm2S_M1_terminal
stfit_m2=hmm2S_M2_terminal

# Viterbi draws for states - m1
zstar_m1 <- rstan::extract(stfit_m1, pars = "zstar_t")[[1]]
zdraw_m1 = apply(zstar_m1, 2, median)

# Viterbi draws for states - m2 
zstar_m2 <- rstan::extract(stfit_m2, pars = "zstar_t")[[1]]
zdraw_m2 = apply(zstar_m2, 2, median)

## attach to stacked dataset Data_final - m1
hmm_data_morph1=Data_final %>% 
  dplyr::select(user, auxvar.1, auxvar.2,timestamp_Clickstream, time, c3, c15, c17, link_depth2, morph, 
                Outcome_DB, Outcome_CV, Outcome_CU, Outcome_RE, Outcome_EA, visit, 
                Cumul_Int_Outcomes, Cumul_Outcome_DB, Cumul_Outcome_CV, Cumul_Outcome_CU,  Cumul_Outcome_RE) %>%
  filter(morph==1) %>%
  arrange(user, time) %>% 
  na.omit() %>% 
  filter(user %in% user_sample_m1$U_m1)%>% mutate(group_id = group_indices(., user)) %>% 
  group_by(user) %>% 
  mutate(T.length=row_number(), lastT=ifelse(row_number()==n(), 1,0)) %>% 
  ungroup()

# add data from Viterbi draws - m1 
hmm_data_morph1%<>% 
  cbind(zdraw_m1)%>% 
  mutate(z_draw=zdraw_m1) %>% 
  dplyr::select(-c(zdraw_m1))

## attach to stacked dataset Data_final - m2
hmm_data_morph2=Data_final %>% 
  dplyr::select(user, auxvar.1, auxvar.2,timestamp_Clickstream, time,c3, c15, c17, link_depth2,morph,
                Outcome_DB, Outcome_CV, Outcome_CU, Outcome_RE, Outcome_EA, visit, 
                Cumul_Int_Outcomes, Cumul_Outcome_DB, Cumul_Outcome_CV, Cumul_Outcome_CU,  Cumul_Outcome_RE) %>%
  filter(morph==2) %>% 
  arrange(user, time) %>% 
  na.omit() %>% 
  filter(user %in% user_sample_m2$U_m2)%>% mutate(group_id = group_indices(., user)) %>% 
  group_by(user) %>% 
  mutate(T.length=row_number(), lastT=ifelse(row_number()==n(), 1,0)) %>% 
  ungroup()

# add data from Viterbi draws - m2
hmm_data_morph2%<>% 
  cbind(zdraw_m2)%>% 
  mutate(z_draw=zdraw_m2) %>% 
  dplyr::select(-c(zdraw_m2))


##########################
#  B) estimate bounce probabilities 
########################## 

mdata=hmm_data_morph1
clicks_sm=mdata %>%
  group_by(user) %>%
  summarise(total.count=n(), 
            last_state=last(z_draw)) %>% 
  filter(last_state==1) %>% 
  dplyr::select(total.count)

## fit a geometric distribution per state per morph - Morph 1 State 1
fs1m1=MASS::fitdistr(clicks_sm$total.count,"geometric")
pS1M1=fs1m1$estimate
p_bounce_s1_m1=1-(1-pS1M1)^(1:20)

## fit a geometric distribution per state per morph - Morph 1 State 2
mdata=hmm_data_morph1
clicks_sm=mdata %>%
  group_by(user) %>%
  summarise(total.count=n(), 
            last_state=last(z_draw)) %>% 
  filter(last_state==2) %>% 
  dplyr::select(total.count)

fs2m1=MASS::fitdistr(clicks_sm$total.count,"geometric")
pS2M1=fs2m1$estimate
p_bounce_s2_m1=1-(1-pS2M1)^(1:20)

## fit a geometric distribution per state per morph - Morph 2 State 1
mdata=hmm_data_morph2
clicks_sm=mdata %>%
  group_by(user) %>%
  summarise(total.count=n(), 
            last_state=last(z_draw)) %>% 
  filter(last_state==1) %>% 
  dplyr::select(total.count)
fs1m2=MASS::fitdistr(clicks_sm$total.count,"geometric")
pS1M2=fs1m2$estimate

## fit a geometric distribution per state per morph - Morph 2 State 2
mdata=hmm_data_morph2
clicks_sm=mdata %>%
  group_by(user) %>%
  summarise(total.count=n(), 
            last_state=last(z_draw)) %>% 
  filter(last_state==2) %>% 
  dplyr::select(total.count)
fs2m2=MASS::fitdistr(clicks_sm$total.count,"geometric")
pS2M2=fs2m2$estimate

# write out results
p_bounce_click=c(pS1M1, pS1M2, pS2M1, pS2M2)
pmf_bounce_user_2Sterminal <- data.frame(Condition = c('S1M1', 'S1M2', 'S2M1', 'S2m2'), coef = p_bounce_click)
write.csv(paste0('Estimates for 2 states/', pmf_bounce_user_2Sterminal), 'pmf_bounce_user_2Sterminal.csv')



##########################
#  C) Calculate the PSM 
########################## 
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))



### Morph 1

# matrices to be filled
zstar_draw <- as.matrix(stfit_m1, pars = c("zstar_t"))
beta_draw <- as.matrix(stfit_m1, pars = c("beta"))
theta_draw<- as.matrix(stfit_m1, pars = c("theta"))

# define n observations and iterations
N=nrow(hmm_data_morph1)
niter=nrow(beta_draw)

# adapt data
hmm_data_morph1 %<>%
  mutate(logCIO=log(Cumul_Int_Outcomes+1),) %>%
  group_by(user) %>% 
  mutate(lagCIO=lag(logCIO, 1)) %>% 
  mutate(lag1CIO=ifelse(is.na(lagCIO), 0,lagCIO)) %>% ungroup()



## Overall averages
phi_tk <- matrix(rep(0, N*niter), N, niter)
for (i in  1:niter)
{
  data = hmm_data_morph1 %>% 
    mutate(iota=rep(1, length(lag1CIO)), beta=beta_draw[i, zstar_draw[i,]], theta=theta_draw[i, zstar_draw[i,]]) %>% 
    rowwise() %>% 
    mutate(phi_tk=inv_logit(c(iota,lag1CIO) %*% c(beta, theta)) )
  
  phi_tk[ ,i]=data$phi_tk
}
phi_m1=apply(phi_tk, 1, mean)

# get final results for morph 1, with confidence interval
emission_m1=rbind(
  c(mean(phi_morph[zdraw_morph==1]), quantile(phi_morph[zdraw_morph==1], probs=c(0.025, 0.975))), 
  c(mean(phi_morph[zdraw_morph==2]), quantile(phi_morph[zdraw_morph==2], probs=c(0.025, 0.975))))



## Morph 2

# matrices to be filled
zstar_draw <- as.matrix(stfit_m2, pars = c("zstar_t"))
beta_draw <- as.matrix(stfit_m2, pars = c("beta"))
theta_draw<- as.matrix(stfit_m2, pars = c("theta"))

# define the N of observations and iterations
N=nrow(dataset_m2)
niter=nrow(beta_draw)

# adapt data
hmm_data_morph2 %<>%
  mutate(logCIO=log(Cumul_Int_Outcomes+1)) %>%
  group_by(user) %>% 
  mutate(lagCIO=lag(logCIO, 1)) %>% 
  mutate(lag1CIO=ifelse(is.na(lagCIO), 0,lagCIO)) %>% 
  ungroup()

# overall averages
phi_tk <- matrix(rep(0, N*niter), N, niter)
for (i in  1:niter)
{
  data = hmm_data_morph2 %>% 
    mutate(iota=rep(1, length(lag1CIO)), beta=beta_draw[i, zstar_draw[i,]], theta=theta_draw[i, zstar_draw[i,]]) %>% 
    rowwise() %>% 
    mutate(phi_tk=inv_logit(c(iota,lag1CIO) %*% c(beta, theta)) )
  
  phi_tk[ ,i]=data$phi_tk
}
phi_m2=apply(phi_tk, 1, mean)

# get final results for morph 2, with confidence interval
emission_m2=rbind(
  c(mean(phi_m2[zdraw_m2==1]), quantile(phi_m2[zdraw_m2==1], probs=c(0.025, 0.975))), 
  c(mean(phi_m2[zdraw_m2==2]), quantile(phi_m2[zdraw_m2==2], probs=c(0.025, 0.975))))

# create dataframe with df_psm 
df_psm <- data.frame(cond = c('m1s1', 'm1s2', 'm2s1', 'm2s2'), value = emission_m1_m2[,1])
