
## for Table 5
delta_sm_average_fixed = colMeans(sim_1k_reps_2S_fixed[,1:(2*TOT_MORPHS)])
delta_sm_average_bandit = colMeans(sim_1k_reps_2S_bandit[,1:(2*TOT_MORPHS)])

delta_sm_average
df_sm_average_fixed <- data.frame(t(matrix(delta_sm_average_fixed, ncol=2)))
colnames(df_sm_average_fixed) <- c('early_state','late_state')
df_sm_average_fixed$total <- rowSums(df_sm_average_fixed)

df_sm_average_bandit <- data.frame(t(matrix(delta_sm_average_bandit, ncol=2)))
colnames(df_sm_average_bandit) <- c('early_state','late_state')
df_sm_average_bandit$total <- rowSums(df_sm_average_bandit)


############
# 
# # load("App1_sim_Random_2ST.RData")
# delta_replicates_baseline=rowSums(sim_1k_reps[,1:(TOT_STATES*TOT_MORPHS)])
# c(mean(delta_replicates_baseline), 
#   t.test(delta_replicates_baseline)$conf.int)
# 
# ## efficiency measures
# (mean(delta_replicates_mab_fixed_state)-mean(delta_replicates_baseline))/mean(delta_replicates_baseline)
# percentage_difference_mab_fixed_state=(delta_replicates_mab_fixed_state-delta_replicates_baseline)/delta_replicates_baseline
# c(mean(percentage_difference_mab_fixed_state), 
#   t.test(percentage_difference_mab_fixed_state)$conf.int)
# 
# (mean(delta_replicates_hmm_bandit)-mean(delta_replicates_baseline))/mean(delta_replicates_baseline)
# percentage_difference_hmmBandit=(delta_replicates_hmm_bandit-delta_replicates_baseline)/delta_replicates_baseline
# c(mean(percentage_difference_hmmBandit), 
#   t.test(percentage_difference_hmmBandit)$conf.int)
