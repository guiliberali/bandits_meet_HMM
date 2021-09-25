#######################################################################################################################
# Author: 
# Purpose: Replication file for table 5, page 29
# 
#
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
# - run the following files in 'study1/3. Simulation Code'
#     - Simulation results for 2 states with fixed states
#     - Simulation results for 2 states with HMM Bandit
# - After that, set working directory to Study2/4. Figures & Tables
# - Then, run this file to generate table 7
#
# Overview
#  A) load the results for each version
#  B) Create table 5
# 
#######################################################################################################################

##########
# A) load the results for each version
##########

delta_sm_average_fixed = colMeans(sim_1k_reps_2S_fixed[,1:(2*TOT_MORPHS)])
delta_sm_average_bandit = colMeans(sim_1k_reps_2S_bandit[,1:(2*TOT_MORPHS)])

##########
# B) Create table 5
##########


df_sm_average_fixed <- data.frame(t(matrix(delta_sm_average_fixed, ncol=2)))
colnames(df_sm_average_fixed) <- c('early_state_fixed_states','late_state_fixed_states')
df_sm_average_fixed$total_fixed_states <- rowSums(df_sm_average_fixed)

df_sm_average_bandit <- data.frame(t(matrix(delta_sm_average_bandit, ncol=2)))
colnames(df_sm_average_bandit) <- c('early_state_HMM_Bandit','late_state_HMM_Banit')
df_sm_average_bandit$total_HMM_Bandit <- rowSums(df_sm_average_bandit)

table5 <- cbind(df_sm_average_fixed, df_sm_average_bandit)
improvement <- (sum(table5$total_fixed_states) - sum(table5$total_HMM_Bandit))/sum(table5$total_HMM_Bandit)
