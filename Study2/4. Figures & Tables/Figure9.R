

#######################################################################################################################
# Author: 
# Purpose: Generates Figure 9 from page 40 of the paper
# 
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
#
#  Overview:
#     A) Plot for 15000 visitors
#     B) Plot for 100000 visitors
#
#######################################################################################################################


SR_data <- read.csv('../1. Raw Data/Real_time_app_success_rate_data.csv')

############
# A) Plot for 15000 visitors
############

SR_plot_15000=SR_data %>% 
  filter(Visitor %in% 100:15000) %>% 
  ggplot(aes(x=Visitor, y=Success_rate, linetype=Treatment, color=Data_type))+
  geom_line()+
  geom_vline(xintercept = c(575), size=.4, color="red", linetype="longdash")+
  geom_vline(xintercept = c(2204), size=.4, color="red")+
  scale_linetype_manual(values = c(1, 3), labels=c("HMM Bandit", "Baseline") ) +
  scale_color_manual(values=c( "black", "gray60") ) +
  labs(y='Purchase rate (running average)',  linetype= "Treatment", color="Data type")+
  theme_classic()+theme(legend.position = "bottom")
SR_plot_15000


############
# A) Plot for 100000 visitors
############

SR_plot_100000=SR_data %>% 
  filter(Visitor %in% 100:100000) %>% 
  ggplot(aes(x=Visitor, y=Success_rate, color=Treatment))+
  geom_line()+
  geom_vline(xintercept = c(575), size=.4, color="red", linetype="longdash")+
  geom_vline(xintercept = c(2204), size=.4, color="red")+
  scale_color_manual(values=c( "black", "gray60"), labels=c("HMM Bandit", "Baseline") ) +
  labs(y='Purchase rate',  color= "Benchmark")+
  theme_classic()+theme(legend.position = "bottom") 

SR_plot_100000

