#######################################################################################################################
#   
# Purpose: Generate Figure 3, page 36  
# 
# Note: Ensure the working directory is Replication_Morphing_HMM/Study1/Figures & Tables
#
# Overview:
#     A) Load in the data, calculate the total and mean number of visits.
#     B) Determine the clickthrough rate, as well as the visitors bouncing off the landing page
#     C) Create the plot
#
#
#######################################################################################################################

###################
# A) Load in the data, calculate the total and mean number of visits.
###################



## Import dataset / run summary stats ##########
Data_final=read.csv("../1. Raw Data/AB_data_05_June_2018.csv", header = T)
## Total number of clicks
nrow(Data_final)
## Total number of visitors
length(unique(Data_final$user))

## Total and mean number of visits  
visits_clicks=Data_final %>%
  group_by(user, morph, visit) %>%
  summarise(total.count=n())
nrow(visits_clicks)
a=visits_clicks %>% group_by(user) %>% summarise(last_visit=last(visit))
mean(a$last_visit)


###################
# B) Determine the clickthrough rate, as well as the visitors bouncing off the landing page
###################

## Overall clickthrough rate
Data_final %>%summarize(rateEA=mean(Outcome_EA))

### Visitors bouncing off the landing page and after the first click
bouncing_clicks=Data_final %>%
  group_by(user, visit, morph) %>% 
  summarise(total.count=n(), bounce_page=last(auxvar.2)) 
bouncing_Off_landing_page=nrow(bouncing_clicks %>% 
                                 filter(bounce_page=='/mba/international-full-time-mba/why-the-rotterdam-mba/'))/nrow(bouncing_clicks)
bouncing_Off_landing_page
bouncing_after_first_click=nrow(bouncing_clicks %>% 
                                  filter(total.count==1))/nrow(bouncing_clicks)
bouncing_after_first_click

## CTR on apply now per click & bounce rate (Figure 3)
## overall ctr user level 
CTR_OutcomeEA_perclick=Data_final %>%  
  group_by(user) %>% 
  mutate(Obs=1:n()) %>% 
  ungroup() %>% group_by(Obs) %>% 
  summarise(rateEA=mean(Outcome_EA))

### overall bounce rate user/visit level
bounce_userlevel=Data_final %>% 
  group_by(morph, user, visit) %>% 
  summarise(total.count=n())

BR=data.frame(click=sort(unique(bounce_userlevel$total.count)), 
              freq=c(with(bounce_userlevel, table(total.count)))) %>% 
  mutate(bounce_rate=freq/sum(freq)) %>% mutate(prob_stay=1-cumsum(bounce_rate))
 
###################
# C) Create the plot
###################

## dataset for plot
EA_BR_df=CTR_OutcomeEA_perclick %>% 
  filter(Obs<100) %>% 
  full_join(BR, by=c("Obs"="click"))

p <- EA_BR_df %>% filter(Obs<100) %>% 
  ggplot(aes(x=Obs))+
  geom_point(aes(y=prob_stay), shape=9)+
  geom_line(aes(y=prob_stay, colour="Probability of continuation"))

p <- p + 
  geom_point(aes(y=rateEA*5))+
  geom_smooth(aes(y=rateEA*5,  colour="CTR on Apply now"), method = "lm") + 
  scale_y_continuous(limits=c(0,1),"Probability of continuing one's visit (1-bounce rate)", sec.axis = sec_axis(~./5, name = "Clickthrough rate on Apply now"))+ 
  scale_colour_manual(values = c("Probability of continuation"="#C1CDCD", "CTR on Apply now"="#525252"))+
  labs(x="Clicks", color="") + theme_light()+
  theme(legend.position = "bottom", legend.direction = "horizontal")
p
