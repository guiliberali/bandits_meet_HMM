#######################################################################################################################
# Author: 
# Purpose: Replication file for Figure 3, page 36 of the paper 
# 
# Note: the results for this table are created with the following steps:
#       1) run all the files in 'study1/simulation codes' folder
#       2) run this file
#
# Overview:
#     Loads in the results from the various simulation files, and determines the average number of applications (with confidence interval)
#
#
#######################################################################################################################


delta_replicates_hmm_bandit_2S_random = rowSums(sim_1k_reps_2S_random[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_2S_fixed = rowSums(sim_1k_reps_2S_fixed[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_2S_bandit = rowSums(sim_1k_reps_2S_bandit[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_random = rowSums(sim_1k_reps_3S_random[,1:(3*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_fixed = rowSums(sim_1k_reps_3S_fixed[,1:(3*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_bandit = rowSums(sim_1k_reps_3S_bandit[,1:(3*TOT_MORPHS)])


calc_mean_CI_applications <- function(delta_replicates){
  
  mean_applications <- mean(delta_replicates)
  CI_applications <- t.test(delta_replicates)$conf.int
  
  C_applications_formatted <- paste0(round(CI_applications[1], 1), '-', round(CI_applications[2], 1))
  
  return(c(mean_applications,C_applications_formatted))

}

n_apps_result_2S_random <- calc_mean_CI_applications(delta_replicates_hmm_bandit_2S_random)
n_apps_result_2S_fixed <- calc_mean_CI_applications(delta_replicates_hmm_bandit_2S_fixed)
n_apps_result_2S_bandit <- calc_mean_CI_applications(delta_replicates_hmm_bandit_2S_bandit)
n_apps_result_3S_random <- calc_mean_CI_applications(delta_replicates_hmm_bandit_3S_random)
n_apps_result_3S_fixed <- calc_mean_CI_applications(delta_replicates_hmm_bandit_3S_fixed)
n_apps_result_3S_bandit <- calc_mean_CI_applications(delta_replicates_hmm_bandit_3S_bandit)


table4 <- data.frame(n_applications_2S = c(n_apps_result_2S_random, n_apps_result_2S_fixed, n_apps_result_2S_bandit),
                     n_applications_3S = c(n_apps_result_3S_random, n_apps_result_3S_fixed, n_apps_result_3S_bandit))