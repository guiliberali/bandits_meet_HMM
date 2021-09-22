
#######################################################################################################################
# Author: 
# Purpose: Calculates Purchasing Rates 
# 
# Note: 
#  -Make sure to have run 'Configuration.R' and 'FUnctions.R' before running this file
#  -After that, ensure the working directory is Replication_Morphing_HMM/Study2/3. Simulation Code
#
# Overview:
#     A) Load in the Survey and clicks data 
#     B) Get performance of each morph (clicks, purchase rate, etc.)
#     C) Compare performance of morph vs random
#
#
#######################################################################################################################

################
# A) Load in the Survey and clicks data 
################

# Day to analyse 
DAY   <- "17"
MONTH <- "05"
DATE     <- paste("2021_",MONTH,"_" , DAY, sep="")

# Location of the raw data
PATH_RAW="../1. Raw Data/" 

### SURVEY DATA ----
survey_data <-read_csv(paste0(PATH_RAW, "/survey_data_", DATE, ".csv"))
survey_data=survey_data %>% 
  mutate(survey_user_id=gsub("[^[:alnum:]]", "", survey_user_id) ) %>% 
  filter(grepl("es|ES", att_check_6_TEXT ) )


### CLICKS DATA ----
processed_clicks_data <- list.files(path=paste0(PATH_RAW, "/processed_clicks"), full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

## filter by participants, then arrange by userid and timestamp
processed_clicks_filtered = processed_clicks_data %>% 
  mutate(user_id=gsub("[^[:alnum:]]", "", user_id)) %>% 
  filter(user_id %in% survey_data$survey_user_id) %>% 
  arrange(ms_timestamp, user_id)

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

################
# B) Get performance of each morph (clicks, purchase rate, etc.)
################

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

##** success rates and people in treatment and control----
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


## create running average of purchase for Figure 
stacked_data_purchase_RA = stacked_data %>% 
  arrange(ms_timestamp) %>% 
  mutate(Obs=rep(1, n())) %>% 
  group_by(cell) %>% 
  mutate(Visitor=1:n(), success_prob_real=cumsum(purchase)/cumsum(Obs)) %>% 
  ungroup() %>% 
  mutate(Treatment=if_else(cell==0,"Random", "HMM_Bandit_Webstore") )

stacked_data_purchase_RA %>% filter(Treatment == 'HMM_Bandit_Webstore')

################
# C) Compare performance of morph vs random
################

getwd()


load(paste0(PATH_RAW, "Random_single_run.RData") )


monitor_Success_Prob=tibble(success_per_visitor=success_per_visitor, 
                            Treatment="Random")

load(paste0(PATH_RAW, "HMM_Bandit_Store_single_run.RData") )
TOT_VISITORS<-100000

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

read.csv(SR_data, "Real_time_app_success_rate_data.csv")



SR_plot=SR_data %>% 
  #filter(Treatment=="Random") %>% 
  filter(Visitor %in% 100:15000) %>% 
  ggplot(aes(x=Visitor, y=Success_rate, linetype=Treatment, color=Data_type))+
  geom_line()+
  #geom_vline(xintercept = c(275, 1066), size=.4, color="gray60")+
  scale_linetype_manual(values = c(1, 3), labels=c("HMM Bandit", "Baseline") ) +
  scale_color_manual(values=c( "black", "gray60") ) +
  labs(y='Purchase rate (running average)',  linetype= "Treatment", color="Data type")+
  theme_classic()+theme(legend.position = "bottom") 
SR_plot




###############
### PLOTS RANDOM vs. various treatments based on success rates ----
monitor_Success_Prob=tibble(success_prob=success_prob_runAverage, 
                            Treatment="Random")

load(paste0( PATH_RAW, "HMM_Bandit_Store_single_run.RData") )
monitor_Success_Prob=monitor_Success_Prob %>% 
  bind_rows(tibble(success_prob=success_prob_runAverage, 
                   Treatment="HMM_Bandit_Webstore") )%>% 
  mutate(Visitor=rep(1:TOT_VISITORS, 5) )


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



## Load performance of random single run
load(paste0(PATH_RAW, "/Random_single_run.RData") )

monitor_Success_Prob=tibble(success_per_visitor=success_per_visitor, 
                            Treatment="Random")

TOT_VISITORS = 100000
getwd()
success_per_visitor
load(paste0(PATH_RAW, "/HMM_Bandit_Store_single_run.RData") )
monitor_Success_Prob=monitor_Success_Prob %>% 
  bind_rows(tibble(success_per_visitor=success_per_visitor, 
                   Treatment="HMM_Bandit_Webstore") )%>% 
  mutate(Visitor=rep(1:TOT_VISITORS, 3) )


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


