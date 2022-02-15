#######################################################################################################################
#  
# Purpose: Table 4, page 28
# 
#
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
#  -After that, run each of the following files from '3. Simulation Code':
#     - Simulation Results for One Morphing Opportunity.R
#     - Simulation Results for Two Morphing Opportunities.R
#     - Simulation Results for Random Morphing.R
# - After that, set working directory to Study2/4. Figures & Tables
# - Then, run this file to generate table 7
#
# Overview
#  A) load the results for each version
#  B) Get the average numer of applications and  confidence interval for each version
#  C) Create table 4
# 
#######################################################################################################################

############
# A) load the results for each version
############

delta_replicates_hmm_bandit_2S_random = rowSums(sim_1k_reps_2S_random[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_2S_fixed = rowSums(sim_1k_reps_2S_fixed[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_2S_bandit = rowSums(sim_1k_reps_2S_bandit[,1:(2*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_random = rowSums(sim_1k_reps_3S_random[,1:(3*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_fixed = rowSums(sim_1k_reps_3S_fixed[,1:(3*TOT_MORPHS)])
delta_replicates_hmm_bandit_3S_bandit = rowSums(sim_1k_reps_3S_bandit[,1:(3*TOT_MORPHS)])

############
# B) Get the average numer of applications and  confidence interval for each version
############

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

############
# C) Create table 4
############

table4 <- data.frame(n_applications_2S = c(n_apps_result_2S_random, n_apps_result_2S_fixed, n_apps_result_2S_bandit),
                     n_applications_3S = c(n_apps_result_3S_random, n_apps_result_3S_fixed, n_apps_result_3S_bandit))
