#######################################################################################################################
# Author: 
# Purpose: Generates table 6 from page 40
# 
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
#  -After that, set working directory to Study2/4. Figures & Tables
# - Then, run this file to generate table 6
#
# Overview
#   A) Load the data on the experiment
#   B) Create table with number of bounces per morph
#
#######################################################################################################################



###############
# A) Load the data on the experiment
###############

## Import data sets ---- 
## use Tlength below for clickstream length 
hmm_dataset = read_csv("../2. HMM Estimates/Application2_Calibration_unstacked_data_April21.csv")


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

###############
# B) Create table with number of bounces per morph
###############

# summarize sample size per condition
##* sample size per condition----
stacked_dataset_with_Tlength=hmm_dataset_T %>% group_by(user_id) %>% 
  slice(n()) %>% ungroup %>% 
  mutate(stayed_7clicks=if_else(Tlength>6, 1, 0))

# get number per condition
n_per_condition_7clicks <- stacked_dataset_with_Tlength %>% 
  group_by(stayed_7clicks, condition_reassigned) %>% 
  summarise(n=n())

n_per_condition_stayed_7clicks <- n_per_condition_7clicks %>% filter(stayed_7clicks == 1)
n_per_condition_before_7clicks <- n_per_condition_7clicks %>% filter(stayed_7clicks == 0)
n_per_second_morph_seen_stayed_7clicks <-matrix(t(n_per_condition_stayed_7clicks$n),nrow=2)
n_per_first_morph_seen_before_7clicks <- matrix(colSums(matrix(n_per_condition_before_7clicks$n,nrow=2)))

## Define Table 6
table6 <- cbind(n_per_first_morph_seen_before_7clicks, n_per_second_morph_seen_stayed_7clicks)
colnames(table6) <- c('Bounced before 7th click', 'Second Morph: 1', 'Second Morph: 2')
rownames(table6) <- c('First Morph: 1', 'First Morph: 2')
table6
