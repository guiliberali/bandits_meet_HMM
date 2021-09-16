
### Replication File for Figure 9, page 40 of the paper 
### Run the 'Calculate Average Purchase Rates.R' File before running any code here.


# For 15000 visitors
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

# For 100000 visitors

##** Figure 9 reported in the manuscript 3rd submission v June 1st 2021----
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

