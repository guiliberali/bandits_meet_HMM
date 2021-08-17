
#######################################################################################################################
# Author: Floris Holstege
# Purpose: Gather data on alpha + beta values and Gittens index from http://mcd.hypermorphingtechnologies.com/HMT/report/22/2021/05/07   
# 
#
# Overview:
#     A) Load necessary packages, and define several functions to download the raw data from the website, and clean it
#     B) Functions to extract alpha + beta values, calculate gittens index
#     C) Create plot with GIttens index
#
#
#######################################################################################################################

################
# A) Load in data (processed and raw), and packages
################
rm(list=ls()) # clean up R environment

# Packes required for subsequent functions. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               stringr,
               qdapRegex, 
               anytime,
               stringr,
               writexl,
               reshape2) 
library(tidyverse)
# get helper functions to download data

#start and end dates to analyse
start <- "2021/05/17"
end <- "2021/05/17"


# # sort the data on ms
# df_cleanData_sorted <- df_cleanData %>% arrange(ms_timestamp)
DAY   <- "17"
MONTH <- "05"
DATE     <- paste("2021_",MONTH,"_" , DAY, sep="")

RAW_DATA="../Raw Data/" ## change as necessary

### SURVEY DATA ----
survey_data <-read_csv(paste0(RAW_DATA, "/survey_data_", DATE, ".csv"))
survey_data=survey_data %>% 
  mutate(survey_user_id=gsub("[^[:alnum:]]", "", survey_user_id) ) %>% 
  filter(grepl("es|ES", att_check_6_TEXT ) )


### CLICKS DATA ----
processed_clicks_data <- list.files(path=paste0(RAW_DATA, "/processed_clicks"), full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

## filter by participants, then arrange by userid and timestamp
processed_clicks_filtered = processed_clicks_data %>% 
  mutate(user_id=gsub("[^[:alnum:]]", "", user_id)) %>% 
  filter(user_id %in% survey_data$survey_user_id) %>% 
  arrange(ms_timestamp, user_id)
dim(processed_clicks_filtered)
length(unique(processed_clicks_filtered$user_id))

winner=processed_clicks_filtered %>% filter(grepl("IgV0dcGD2WKZSpAJoEcRHTbqGQbY5", user_id))


## check last alpha/ beta
test_abUpdate=processed_clicks_data %>% arrange(ms_timestamp) %>% 
  dplyr::select(aux2) %>% 
  tail()

test = processed_clicks_data %>% group_by(user_id) %>% 
  slice(n()) %>% 
  # summarise(a=sum(cell), last_morph=last(morph_id), purchase=sum(outcome0)) %>% 
  group_by(morph_id, cell) %>% 
  summarise(purchase_rate=mean(outcome0), count=n())

processed_clicks_data %>% filter(outcome3 == 1) %>% 
  group_by(cell, morph_id) %>% 
  summarise(cell1_prop=mean(cell), purchase_rate=mean(outcome0), count_last_m=n())

## add condiitons variable
clicks_data_with_conditions= processed_clicks_filtered %>%
  filter(aux0!=(-1)) %>%  ## removes clicks on instructions 
  group_by(user_id) %>% 
  mutate(total_clicks=last(aux0+1), first_morph=first(morph_id), 
         last_morph=last(morph_id) ) %>% 
  ungroup() %>% 
  dplyr::select(-c('...1', device, source, starts_with("segment"), outcome4, outcome5, outcome6, outcome7, page_id, link_id)) %>% 
  mutate(condition=as.factor(case_when(first_morph==1 & last_morph==1 ~ 1, 
                                       first_morph==1 & last_morph==2 ~ 2,
                                       first_morph==2 & last_morph==1 ~ 3,
                                       TRUE ~ 4) ))

## performance of morph combo
stacked_data = clicks_data_with_conditions %>% 
  group_by(user_id) %>% 
  summarise(condition=first(condition), first_morph=first(first_morph), last_morph=first(last_morph), 
            cell=first(cell), purchase=last(outcome0), total_clicks=last(aux0+1),
            ms_timestamp=last(ms_timestamp)) %>%
  mutate(stayed_for_7clicks=if_else(total_clicks<7, 0, 1)) %>% 
  ungroup() %>% 
  arrange(ms_timestamp) %>% 
  group_by(cell) %>% 
  mutate(Visitor=1:n()) %>% 
  ungroup()

#### total clicks -----
funnel_reg=glm(formula = total_clicks ~ as.factor(condition), 
               family="poisson",
               data = stacked_data)
summary(funnel_reg)

##** success rates and people in teratment and control----
stacked_data %>% group_by(condition) %>% 
  filter(cell==0 & total_clicks>6) %>% 
  summarise(n())

stacked_data %>% group_by(cell) %>% 
  summarise(n(), mean(purchase))

## purchase rates ----
conditions_purchase_rate = stacked_data %>%  
  group_by(cell, condition) %>%
  summarise(purchase_rate=mean(purchase), count_n=n(), mean_nr_clicks=mean(total_clicks), sd_clicks=sd(total_clicks))
conditions_purchase_rate   

purchase_reg=glm(formula = purchase ~ as.factor(cell), 
               family="binomial",
               data = stacked_data )

conditions_purchase_rate = stacked_data %>%  
  group_by(cell, condition, first_morph, last_morph) %>%
  summarise(purchase_rate=mean(purchase), sd=sd(purchase), count_n=n(), mean_nr_clicks=mean(total_clicks))
   

## purchase rate and burnin ----
conditions_purchase_rate = stacked_data %>%  
  filter(cell==1) %>% 
  filter(Visitor %in% 1:2000) %>% #(n()/2):n() ) ) %>% 
  summarise(purchase_rate=mean(purchase), count_n=n(), mean_nr_clicks=mean(total_clicks))
 

conditions_purchase_rate = stacked_data %>%  
  filter(!(cell==1 & Visitor %in% 1:500)) %>% 
  group_by(cell) %>%  summarise(purchase_rate=mean(purchase))


conditions_purchase_rate = stacked_data %>%  
  filter(cell==1 & Visitor %in% 1:50) %>% 
  group_by(cell, condition) %>%  
  summarise(purchase_rate=mean(purchase), mean_nr_clicks=mean(total_clicks), 
            count_n=n())


## create running average of purchase for Figure XYZ ---- 
stacked_data_purchase_RA = stacked_data %>% 
  arrange(ms_timestamp) %>% 
  mutate(Obs=rep(1, n())) %>% 
  group_by(cell) %>% 
  mutate(Visitor=1:n(), success_prob_real=cumsum(purchase)/cumsum(Obs)) %>% 
  ungroup() %>% 
  mutate(Treatment=if_else(cell==0,"Random", "HMM_Bandit_Webstore") )

plot(stacked_data_purchase_RA$success_prob_real[stacked_data_purchase_RA$cell==1], type="l")

### PLOT OF RUNNING AV OF SUCCESSES WITH SIMS -----
load(paste0(RAW_DATA, "/Random_single_run.RData") )

### PLOTS RANDOM vs. various treatments based on success rates ----
monitor_Success_Prob=tibble(success_prob=success_prob_runAverage, 
                            Treatment="Random")

load(paste0(RAW_DATA, "/HMM_Bandit_Store_single_run.RData") )

TOT_VISITORS = 100000
monitor_Success_Prob=monitor_Success_Prob %>% 
  bind_rows(tibble(success_prob=success_prob_runAverage, 
                   Treatment="HMM_Bandit_Webstore") ) %>%
  mutate(Visitor=rep(1:TOT_VISITORS, 2) )


monitor_Success_Prob %>% 
  filter(Visitor %in% 5000:20000) %>% 
  ggplot(aes(x=Visitor, y=success_prob, color=Treatment))+
  geom_line()

blended_data=monitor_Success_Prob %>% 
  left_join(stacked_data_purchase_RA) %>% 
  mutate(Data_type=if_else(is.na(success_prob_real), "Simulated", "Observed") ) %>% 
  mutate(Success_rate=if_else(is.na(success_prob_real), success_prob, success_prob_real))
  
SR_plot=blended_data %>% 
  filter(Visitor %in% 25:5000) %>% 
  ggplot(aes(x=Visitor, y=Success_rate, linetype=Treatment, color=Data_type))+
  geom_line()+
  geom_vline(xintercept = c(275, 1066), size=.4, color="gray60")+
  scale_linetype_manual(values = c(1, 3), labels=c("HMM Bandit (Webstore)", "Random") ) +
  scale_color_manual(values=c( "black", "red") ) +
  labs(y='Purchase rate (running average)',  linetype= "Treatment", color="Data type")+
  theme_classic()+theme(legend.position = "bottom") 

SR_plot


plot_scale <- 3.5
plot_aspect <- 2.5
save_plot <- purrr::partial(ggsave, width = plot_aspect * plot_scale, height = 1 * plot_scale)
PATH_PLOTS ="~/Documents/Algo_study_May2021/plots"

save_plot(paste0(PATH_PLOTS, '/Purchase_rate_running_average.pdf') )


### PLOTS RANDOM vs. various treatments based on success dummies----
## based on success (0, 1)
load(paste0(RAW_DATA, "/Random_single_run.RData") )

monitor_Success_Prob=tibble(success_per_visitor=success_per_visitor, 
                            Treatment="Random")

TOT_VISITORS = 100000

load(paste0(RAW_DATA, "/HMM_Bandit_Store_single_run.RData") )
monitor_Success_Prob=monitor_Success_Prob %>% 
  bind_rows(tibble(success_per_visitor=success_per_visitor, 
                   Treatment="HMM_Bandit_Webstore") )%>% 
  mutate(Visitor=rep(1:TOT_VISITORS, 2) )

monitor_Success_Prob %>% 
  filter(Visitor %in% 5000:20000) %>% 
  ggplot(aes(x=Visitor, y=success_per_visitor, color=Treatment))+
  geom_line()

blended_data=monitor_Success_Prob %>% 
  left_join(stacked_data_purchase_RA) %>% 
  mutate(Data_type=if_else(is.na(purchase), "Simulated", "Observed") ) %>% 
  mutate(Success_dummy=if_else(is.na(purchase), success_per_visitor, purchase)) 

SR_data = blended_data %>%  group_by(Treatment) %>% 
  mutate(Success_rate=cumsum(Success_dummy)/Visitor) %>% 
  ungroup() 


