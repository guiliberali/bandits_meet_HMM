

#######################################################################################################################
# Author: 
# Purpose: Generates Figure 8 from page 38 of the paper
# 
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
#  -After that, run '3. Simulation Code/Results for One MorphingOpportunity.R'
#
#  Overview:
#     A) Monitor Gittens Index for morph 1 vs morph 2 - POST
#     B) Plot the Gittens Index
#
#
#######################################################################################################################



########
# A) Monitor Gittens Index for morph 1 vs morph 2 - POST
########

colnames(monitor_G_m1_s_POST)<- colnames(monitor_G_m2_s_POST)<-c("State 1", "State 2")
monitor_G=as_tibble(monitor_G_m1_s_POST) %>% 
  mutate(Visitor=1:TOT_VISITORS) %>% 
  pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
  mutate(Morph=c("Morph 1")) %>% 
  bind_rows(as_tibble(monitor_G_m2_s_POST) %>% 
              mutate(Visitor=1:TOT_VISITORS) %>% 
              pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
              mutate(Morph=c("Morph 2")) ) %>% 
  mutate(state_morph=case_when(State== 'State 1' & Morph== "Morph 1" ~ "S1M1",
                               State== 'State 2' & Morph== "Morph 1" ~ "S2M1",
                               State== 'State 1' & Morph== "Morph 2" ~ "S1M2",
                               TRUE ~ "S2M2"))


########
# B) PLot the Gittens Index
########

GI_plot=monitor_G%>% 
  filter(Visitor %in% 1:15000) %>% 
  ggplot(aes(x=Visitor, y=Index, color=state_morph))+
  geom_line( )+
  scale_color_manual(labels=c( "State 1 / Morph 1", "State 1 / Morph 2", "State 2 / Morph 1", "State 2 / Morph 2"), 
                     values=c( "darkred", "firebrick2", "darkblue", "blue1") ) +
  labs(y='Gittins Index', color="State/Morph")+
  theme_classic()+theme(legend.position = "bottom") 

GI_plot
