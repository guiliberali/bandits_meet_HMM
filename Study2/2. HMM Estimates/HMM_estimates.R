#######################################################################################################################
# 
# Author: 
# Purpose: Get HMM estimates for several parameters used in the simulations, per condition
# 
# Note: 
#  -Run 'Configuration.R' before running this file
#  -After that, ensure the working directory is Replication_Morphing_HMM/Study2/2. HMM Estimates
# 
# Overview:
#     A) Read in dataset, initial transformations, set parameters
#     B) Make HMM estimates with Stan
#     C) Generate several parameters for simulations: 
#           - rho_c14_2ST_RCT_April2021.csv
#           - mu_c14_2ST_RCT_April2021.csv
#     D) Get state each user was in at each click
#     E) Generate the probability of being in a state across clicks
#     F) Generate several parameters for simulations: 
#           -  psm_condition_user_level_RCT_April2021_prePost.csv
#           -  pmf_bounce_geometric_RCT_April2021.csv 
#           - geometric_emission_probs_RCT_April2021.csv
#    
#######################################################################################################################

###########
# A) Read in dataset, initial transformations, set parameters
############

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Source functions ----
source("HMM_functions_Application2.R" ) 
## Compile the HMM model
m = stan_model(model_code = hmm_terminal_geometric_ind_cov)
devtools::session_info("rstan")


## Import data sets ---- 
## use Tlength below for clicksstream length 
hmm_dataset = read_csv("Application2_Calibration_unstacked_data_April21.csv")
participant_stacked_data = read_csv("Application2_Calibration_stacked_data_April2021.csv" ) 


### Build data set ----
##* Create "add-to-comparison" dummy variable ----
hmm_dataset_T=hmm_dataset %>%  
  mutate(outcome_addTo=case_when(grepl(".*/add-to-comparison", destination_url) ~ 1,
                                 TRUE ~ 0)) %>% 
  group_by(user_id) %>% 
  mutate(Cumul_outcome_addTo = cumsum(outcome_addTo), Tlength=row_number()) %>% 
  mutate(total_clicks=last(Tlength)) %>% 
  ungroup() %>% 
  ## transform cumsum into a dummy vriable
  mutate(Cumul_outcome_addTo= if_else(Cumul_outcome_addTo > 0, 1, 0) )


# summarize sample size per condition
##* sample size per condition----
stacked_dataset_with_Tlength=hmm_dataset_T %>% group_by(user_id) %>% 
  slice(n()) %>% ungroup %>% 
  mutate(stayed_7clicks=if_else(Tlength>6, 1, 0))


#### Regressions for Section 6.2.3 self-selection -----
funnel_reg=glm(formula = bounce ~ as.factor(condition), 
               family="binomial",
               data = participant_stacked_data %>% filter(total_clicks<7))
summary(funnel_reg)

funnel_reg=glm(formula = total_clicks ~ as.factor(condition), 
               family="poisson",
               data = participant_stacked_data %>% filter(total_clicks<7))
summary(funnel_reg)


## prob=.229 chance of clicking on "add-to-comparison" by click 7 (including).
## used in the simulation study
prob_cov_click7=hmm_dataset_T %>% filter(Tlength<8) %>% 
  group_by(user_id) %>% 
  summarise(addTo_cov=sum(Cumul_outcome_addTo)) %>% 
  mutate(addTo_cov_user= if_else(addTo_cov > 0, 1, 0) ) %>% 
  summarise(rate_click=mean(addTo_cov_user))


###########
# B) Make HMM estimates with Stan
###########

HMM_estimate_condition <- function(hmm_dataset_T, condition){
  
  
  
  ## no cov on emission prob, "add-to-comparison" covariate on transition probs
  #* Construct condition-level hmm dataset ----
  ## Change condition_reassigned to 2, 3, and 4 to estimate the HMM models for each condition. 
  hmm_data_condition=hmm_dataset_T[hmm_dataset_T$condition_reassigned== condition,]
  
  if(condition==1){
    hmm_data_condition <- hmm_data_condition %>%filter(user_id!= "GVy6y1x89W8LKRRHpFAAEyam587pz")
    
  }
  if(condition==4){
    hmm_data_condition <- hmm_data_condition %>%filter(user_id!= "0wMm9bCTVEAlsXnszDlSJxLO63JtdeB")
    
  }
  
  
  hmm_data_condition <- hmm_data_condition %>% 
    filter(total_clicks!=1) %>% 
    filter(!grepl("finished-checking-out", destination_url)) %>% ## remove end-page with explanation in the clickstream
    group_by(user_id) %>% 
    mutate(IND = cur_group_id(),
           lastT=ifelse(row_number()==n(), 1, 0) ) %>% 
    ungroup() 
  
  XA=as.matrix(hmm_data_condition  %>%
                 dplyr::select(Cumul_outcome_addTo))
  
  
  ## number of clicks per product, per individual
  N=nrow(hmm_data_condition)
  Nind=length(unique(hmm_data_condition$user_id))
  
  
  
  K_EST=2
  stan.data = list(
    N=N,                   # number of observations
    Nind=Nind,             # number of individuals
    K=K_EST,                # number of hidden states  
    A=ncol(XA),
    IND=pull(hmm_data_condition, IND),
    T=pull(hmm_data_condition, Tlength), 
    lastT=pull(hmm_data_condition, lastT),  
    XA= XA,  
    z=1:K_EST, 
    lastQij=c(rep(0, K_EST-1), 1),
    pi0=c(1, rep(0, K_EST-1))   #start probability, everyone starts in the early state
  )
  
  hmm2ST_geometric_ind_cov_compOnly <- stan(model_code = hmm_terminal_geometric_ind_cov, data = stan.data,
                                                verbose = TRUE, chains = 2, seed = 123456789, 
                                                iter = 50, warmup = 40, 
                                                control = list(adapt_delta=0.99, max_treedepth=15) )
  
  return(hmm2ST_geometric_ind_cov_compOnly)
  
  
  
}

# get HMM estimates per condition
hmm2ST_geometric_ind_c1r_cov_compOnly <- HMM_estimate_condition(hmm_dataset_T = hmm_dataset_T, 1)
hmm2ST_geometric_ind_c2r_cov_compOnly <- HMM_estimate_condition(hmm_dataset_T = hmm_dataset_T, 2)
hmm2ST_geometric_ind_c3r_cov_compOnly <- HMM_estimate_condition(hmm_dataset_T = hmm_dataset_T, 3)
hmm2ST_geometric_ind_c4r_cov_compOnly <- HMM_estimate_condition(hmm_dataset_T = hmm_dataset_T, 4)

hmm_dataset_T %>%
  filter(condition ==3)

###################
# C) Generate several parameters for simulations: 
#             - rho_c14_2ST_RCT_April2021.csv
#             - mu_c14_2ST_RCT_April2021.csv
###################

##* tables with estimates for Appendix C ----
mu_m1 <- summary(hmm2ST_geometric_ind_c1r_cov_compOnly, pars = c("mu_bar"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
mu_m2 <- summary(hmm2ST_geometric_ind_c2r_cov_compOnly, pars = c("mu_bar"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
mu_m3 <- summary(hmm2ST_geometric_ind_c3r_cov_compOnly, pars = c("mu_bar"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
mu_m4 <- summary(hmm2ST_geometric_ind_c4r_cov_compOnly, pars = c("mu_bar"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]

rho_m1 <- summary(hmm2ST_geometric_ind_c1r_cov_compOnly, pars = c("rho"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
rho_m2 <- summary(hmm2ST_geometric_ind_c2r_cov_compOnly, pars = c("rho"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
rho_m3 <- summary(hmm2ST_geometric_ind_c3r_cov_compOnly, pars = c("rho"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
rho_m4 <- summary(hmm2ST_geometric_ind_c4r_cov_compOnly, pars = c("rho"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]

beta_m1 <- summary(hmm2ST_geometric_ind_c1r_cov_compOnly, pars = c("beta_bar"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
beta_m1[2,]=beta_m1[1,]+exp(beta_m1[2,])
beta_m2 <- summary(hmm2ST_geometric_ind_c2r_cov_compOnly, pars = c("beta_bar"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
beta_m2[2,]=beta_m2[1,]+exp(beta_m2[2,])
beta_m3 <- summary(hmm2ST_geometric_ind_c3r_cov_compOnly, pars = c("beta_bar"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
beta_m3[2,]=beta_m3[1,]+exp(beta_m3[2,])
beta_m4 <- summary(hmm2ST_geometric_ind_c4r_cov_compOnly, pars = c("beta_bar"), probs=c(0.025, 0.975), digits=3)$summary[,c(1, 4, 5)]
beta_m4[2,]=beta_m4[1,]+exp(beta_m4[2,])


pars=rbind(c(mu_m1, mu_m2), c(rho_m1, rho_m2), c(beta_m1[1,], beta_m2[1,]), c(beta_m1[2,], beta_m2[2,]),
           c(mu_m3, mu_m4), c(rho_m3, rho_m4), c(beta_m3[1,], beta_m4[1,]), c(beta_m3[2,], beta_m4[2,]))
pars

stargazer(pars, summary = F, align=T)


##########
#     D) Get state each user was in at each click
##########


##** CONDITION 1----
hmm_data_condition=hmm_dataset_T %>% 
  filter(condition_reassigned==1) %>% 
  filter(user_id!= "GVy6y1x89W8LKRRHpFAAEyam587pz") %>% # for condition 1 "GVy6y1x89W8LKRRHpFAAEyam587pz") %>% #for condition 4 #"0wMm9bCTVEAlsXnszDlSJxLO63JtdeB") %>% 
  filter(total_clicks!=1) %>% 
  filter(!grepl("finished-checking-out", destination_url)) %>%
  group_by(user_id) %>% 
  mutate(IND = cur_group_id(),
         lastT=ifelse(row_number()==n(), 1, 0) ) %>% 
  ungroup()

alpha_tk <- rstan::extract(hmm2ST_geometric_ind_c1r_cov_compOnly, pars = 'alpha_tk')[[1]]
prob_state = apply(alpha_tk, c(2, 3), median)
estimated = apply(apply(alpha_tk, c(2, 3), median), 1, which.is.max)
hmm_data_with_states_c1 = hmm_data_condition %>%  
  bind_cols(est_state_2ST = estimated, 
            prob_state1_2ST = round(prob_state[,1], 2), 
            prob_state2_2ST = round(prob_state[,2],2))



##** CONDITION 2----
hmm_data_condition=hmm_dataset_T %>% 
  filter(condition_reassigned==2) %>% 
  #filter(user_id!= "GVy6y1x89W8LKRRHpFAAEyam587pz") %>% # for condition 1 "GVy6y1x89W8LKRRHpFAAEyam587pz") %>% #for condition 4 #"0wMm9bCTVEAlsXnszDlSJxLO63JtdeB") %>% 
  filter(total_clicks!=1) %>% 
  filter(!grepl("finished-checking-out", destination_url)) %>%
  group_by(user_id) %>% 
  mutate(IND = cur_group_id(),
         lastT=ifelse(row_number()==n(), 1, 0) ) %>% 
  ungroup() %>%  
  mutate(trolley_dummy=if_else(grepl("/trolley", destination_url), 1, 0))

alpha_tk <- rstan::extract(hmm2ST_geometric_ind_c2r_cov_compOnly, pars = 'alpha_tk')[[1]]
prob_state = apply(alpha_tk, c(2, 3), median)
estimated = apply(apply(alpha_tk, c(2, 3), median), 1, which.is.max)
hmm_data_with_states_c2 = hmm_data_condition %>%  
  bind_cols(est_state_2ST = estimated, 
            prob_state1_2ST = round(prob_state[,1], 2), 
            prob_state2_2ST = round(prob_state[,2],2))

##* CONDITION 3 ----
#** State recovery ----
## 2 states terminal
hmm_data_condition=hmm_dataset_T %>% 
  filter(condition_reassigned==3) %>% 
  #filter(user_id!= "GVy6y1x89W8LKRRHpFAAEyam587pz") %>% # for condition 1 "GVy6y1x89W8LKRRHpFAAEyam587pz") %>% #for condition 4 #"0wMm9bCTVEAlsXnszDlSJxLO63JtdeB") %>% 
  filter(total_clicks!=1) %>% 
  filter(!grepl("finished-checking-out", destination_url)) %>%
  group_by(user_id) %>% 
  mutate(IND = cur_group_id(),
         lastT=ifelse(row_number()==n(), 1, 0) ) %>% 
  ungroup() %>%  
  mutate(trolley_dummy=if_else(grepl("/trolley", destination_url), 1, 0))

alpha_tk <- rstan::extract(hmm2ST_geometric_ind_c3r_cov_compOnly, pars = 'alpha_tk')[[1]]
prob_state = apply(alpha_tk, c(2, 3), median)
estimated = apply(apply(alpha_tk, c(2, 3), median), 1, which.is.max)
hmm_data_with_states_c3 = hmm_data_condition %>%  
  bind_cols(est_state_2ST = estimated, 
            prob_state1_2ST = round(prob_state[,1], 2), 
            prob_state2_2ST = round(prob_state[,2],2))

##* CONDITION 4 ----
#** State recovery ----
## 2 states terminal
hmm_data_condition=hmm_dataset_T %>% 
  filter(condition_reassigned==4) %>% 
  filter(user_id!= "0wMm9bCTVEAlsXnszDlSJxLO63JtdeB") %>% # for condition 1 "GVy6y1x89W8LKRRHpFAAEyam587pz") %>% #for condition 4 #"0wMm9bCTVEAlsXnszDlSJxLO63JtdeB") %>% 
  filter(total_clicks!=1) %>% 
  filter(!grepl("finished-checking-out", destination_url)) %>%
  group_by(user_id) %>% 
  mutate(IND = cur_group_id(),
         lastT=ifelse(row_number()==n(), 1, 0) ) %>% 
  ungroup() %>%  
  mutate(trolley_dummy=if_else(grepl("/trolley", destination_url), 1, 0))

alpha_tk <- rstan::extract(hmm2ST_geometric_ind_c4r_cov_compOnly, pars = 'alpha_tk')[[1]]
prob_state = apply(alpha_tk, c(2, 3), median)
estimated = apply(apply(alpha_tk, c(2, 3), median), 1, which.is.max)
hmm_data_with_states_c4 = hmm_data_condition %>%  
  bind_cols(est_state_2ST = estimated, 
            prob_state1_2ST = round(prob_state[,1], 2), 
            prob_state2_2ST = round(prob_state[,2],2))

##########
# E) Generate the probability of being in a state across clicks
##########

hmm_data_with_states_cond=hmm_data_with_states_c1 %>% 
  bind_rows(hmm_data_with_states_c2) %>% 
  bind_rows(hmm_data_with_states_c3) %>% 
  bind_rows(hmm_data_with_states_c4)

a=hmm_data_with_states_cond %>% 
  dplyr::select(Tlength, prob_state1_2ST, prob_state2_2ST, outcome_addTo, condition_reassigned) %>% 
  group_by(Tlength, outcome_addTo, condition_reassigned) %>% 
  summarise(pb1=mean(prob_state1_2ST), pb2=mean(prob_state2_2ST)) %>% 
  ungroup() 
colnames(a)=c("Click", "add_to_comp", "Condition", "Early state", "Late state")

a=a %>%  
  mutate(add_to_comp=if_else(add_to_comp==0, "Baseline", "Click on covariate" ), 
         Condition = case_when(Condition==1~"m1m1",
                               Condition==2~"m1m2",
                               Condition==3~"m2m1",
                               TRUE~"m2m2")) %>% 
  pivot_longer(-c(Click, add_to_comp, Condition), names_to = "State", values_to = "Prob")
a$add_to_comp=factor(a$add_to_comp, levels = c("Baseline", "Click on covariate"))

alphaTR=a %>% 
  filter(Click<15) %>% 
  ggplot(aes(x=Click, y=Prob, fill=State))+
  geom_bar(stat="identity")+
  labs(y="State membership probability", fill="")+
  scale_x_continuous(breaks=c(1, 4, 7, 10, 14))+
  scale_fill_manual(values=c("#C1CDCD", "#454545"))+
  theme_light()+theme(legend.position = "bottom") +
  facet_grid(add_to_comp~Condition)

##################
#     E) Generate several parameters for simulations: 
#           -  psm_condition_user_level_RCT_April2021_prePost.csv
#           -  pmf_bounce_geometric_RCT_April2021.csv 
#           - geometric_emission_probs_RCT_April2021.csv
###################

## success probs for those who stayed for less than 7 clicks
success_prob_before_c7=stacked_dataset_with_Tlength %>% 
  filter(stayed_7clicks==0) %>% 
  group_by(condition) %>% ## because this includes only those who saw M1 only or M2 only, in conditions 1 or 4
  summarise(success_prob=mean(outcome_var))

success_prob_after_c7=stacked_dataset_with_Tlength %>% 
  filter(stayed_7clicks==1) %>% 
  group_by(condition_reassigned) %>% 
  summarise(success_prob=mean(outcome_var))

## used as success probabilities in the simulations - (psm_condition_user_level_RCT_April2021_prePost.csv) 
psm_condition_user_level_RCT_April2021_prePost <- success_prob_after_c7
psm_condition_user_level_RCT_April2021_prePost

###* Bounce rate per conditon, necessary for the sims (pmf_bounce_geometric_RCT_April2021.csv)
last_click_distr = hmm_data_with_states_cond %>% 
  group_by(user_id) %>% 
  slice(n()) %>% 
  ungroup() %>% 
  dplyr::select(user_id, Tlength) %>% 
  mutate(last_click_number=Tlength) %>% 
  dplyr::select(-Tlength) %>% 
  left_join(participant_stacked_data %>% 
              dplyr::select(user_id, outcome_var, condition_reassigned)) 

## fit a geometric distribution to evaluate the bounce rate at every click, conditions 1 and 4
geom_distr=MASS::fitdistr(last_click_distr$last_click_number[last_click_distr$condition_reassigned==1],"geometric")
geom_distr$estimate
geom_distr=MASS::fitdistr(last_click_distr$last_click_number[last_click_distr$condition_reassigned==4],"geometric")
geom_distr$estimate


setwd('/Users/flori/OneDrive/Documents/Github/Replication_Morphing_HMM/Study2/2. HMM Estimates')
