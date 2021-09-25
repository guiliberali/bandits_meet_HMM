
#######################################################################################################################
# Author: 
# Purpose: Generates table 7 from page 42
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
#  A) set the perfect information benchmark
#  B) Get the average purchase rate and confidence interval for each version
#  C) Create table 7
#
#
#
#######################################################################################################################


############
# A) set the perfect information benchmark
###########

## perfect information benchmark
load('../1. Raw Data/HMM_Bandit_Store_single_run.RData')
last_click_dummy=ifelse(last_click > 6, 1, 0)
pr_two_morph_exposure=sum(last_click_dummy)/TOT_VISITORS
perfect_info_benchmark=(1-pr_two_morph_exposure)*P_mt_pre[2] + pr_two_morph_exposure*P_mt_post[3] 

############
# B) Get the average purchase rate and confidence interval for each version
###########

# Average Purchase Rate and Confidence Interval - Random Morphing (Baseline)
purchase_rates_baseline <- rowSums(sim_1krep_random_morph [,-1])/TOT_VISITORS
avg_purchase_rate_baseline <- mean(purchase_rates_baseline)
CI_purchase_rate_basline <- t.test(purchase_rates_baseline)$conf.int

# Average Purchase Rate and Confidence Interval - One Morphing Opportunity
purchase_rates_one_morph <- rowSums(sim_1krep_one_morph [,-1])/TOT_VISITORS
avg_purchase_rate_one_morph <- mean(purchase_rates_one_morph)
CI_purchase_rate_one_morph <- t.test(purchase_rates_one_morph)$conf.int

# Average Purchase Rate and Confidence Interval - Two Morphing Opportunities
purchase_rates_two_morph <- rowSums(sim_1krep_two_morph[,-1])/TOT_VISITORS
avg_purchase_rate_two_morph <- mean(purchase_rates_two_morph)
CI_purchase_rate_two_morph <- t.test(purchase_rates_two_morph)$conf.int

############
# C) Create table 7
###########

# Calculate Efficiency of one and two morphing opportunities
efficiency_one_morph <- (avg_purchase_rate_one_morph-avg_purchase_rate_baseline)/(perfect_info_benchmark-avg_purchase_rate_baseline)
efficiency_two_morph <- (avg_purchase_rate_two_morph-avg_purchase_rate_baseline)/(perfect_info_benchmark-avg_purchase_rate_baseline)

# Format the confidence intervals for the table
CI_purchase_rate_random_morph_format <- paste0(round(CI_purchase_rate_random_morph[1],4), '-',round(CI_purchase_rate_random_morph[2],4))
CI_purchase_rate_one_morph_format <- paste0(round(CI_purchase_rate_one_morph[1],4), '-',round(CI_purchase_rate_one_morph[2],4))
CI_purchase_rate_two_morph_format <- paste0(round(CI_purchase_rate_two_morph[1],4), '-',round(CI_purchase_rate_two_morph[2],4))

# put all of the information together in the table
table7 <- data.frame(Avg_purchase_rate = c(avg_purchase_rate_baseline,
                                           avg_purchase_rate_one_morph, 
                                           avg_purchase_rate_two_morph,
                                           perfect_info_benchmark ),
                     CI_purchase_rate = c(CI_purchase_rate_random_morph_format,
                                          CI_purchase_rate_one_morph_format,
                                          CI_purchase_rate_two_morph_format,
                                          NA),
                     Efficiency = c(0, efficiency_one_morph,efficiency_two_morph, 1 ))

table7