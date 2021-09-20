
## 4. Summary statistics ----
## for Table 4
# load("App1_sim_hmmBandit_2ST.RData")


delta_replicates_hmm_bandit
delta_replicates_hmm_bandit_2S_random = rowSums(sim_1k_reps_2S_random[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_2S_fixed = rowSums(sim_1k_reps_2S_fixed[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_2S_bandit = rowSums(sim_1k_reps_2S_bandit[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_random = rowSums(sim_1k_reps_3S_random[,1:(3*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_fixed = rowSums(sim_1k_reps_3S_fixed[,1:(3*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_bandit = rowSums(sim_1k_reps_3S_bandit[,1:(3*TOT_MORPHS)])


calc_mean_CI_applications <- function(delta_replicates){
  
  mean_applications <- mean(delta_replicates)
  CI_applications <- t.test(delta_replicates)$conf.int
  
  C_applications_formatted <- paste0(round(CI_applications[1], 4), '-', round(CI_applications[2], 4))

}

avg_applications_2S_random <- mean(delta_replicates_hmm_bandit_2S_random)
CI_applications_2S_random <- t.test(delta_replicates_hmm_bandit_2S_random)$conf.int

avg_applications_2S_fixed <- mean(delta_replicates_hmm_bandit_2S_fixed)
CI_applications_2S_fixed <- t.test(delta_replicates_hmm_bandit_2S_fixed)$conf.int



c(mean(delta_replicates_hmm_bandit), 
  t.test(delta_replicates_hmm_bandit)$conf.int)

## for Table 5
delta_sm_average = colMeans(sim_1k_reps[,1:(TOT_STATES*TOT_MORPHS)])

# load("App1_sim_MabFixedState_2ST.RData")
delta_replicates_mab_fixed_state=rowSums(sim_1k_reps[,1:(TOT_STATES*TOT_MORPHS)])
c(mean(delta_replicates_mab_fixed_state), 
  t.test(delta_replicates_mab_fixed_state)$conf.int)
