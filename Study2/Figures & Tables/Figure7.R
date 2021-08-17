

### Replication File for Figure 7, page 36 of the paper 


RAW_DATA="../Raw Data"

#* State membership probs plots ----
###### membership change over clicks
hmm_data_with_states_c14=read_csv(paste0(RAW_DATA,"/hmm_data_with_states_c14.csv"))
hmm_data_with_states_c23=read_csv(paste0(RAW_DATA,"/hmm_data_with_states_c23.csv"))
hmm_data_with_states_cond=hmm_data_with_states_c14 %>% 
  bind_rows(hmm_data_with_states_c23)

a=hmm_data_with_states_cond %>% 
  dplyr::select(Tlength, prob_state1_2ST, prob_state2_2ST, outcome_addTo, condition_reassigned) %>% 
  group_by(Tlength, outcome_addTo, condition_reassigned) %>% 
  summarise(pb1=mean(prob_state1_2ST), pb2=mean(prob_state2_2ST)) %>% 
  ungroup() 
colnames(a)=c("Click", "add_to_comp", "Condition", "Early state", "Late state")

# Changes per condition
a=a %>%  
  mutate(add_to_comp=if_else(add_to_comp==0, "Baseline", "Click on covariate" ), 
         Condition = case_when(Condition==1~"m1m1",
                               Condition==2~"m1m2",
                               Condition==3~"m2m1",
                               TRUE~"m2m2")) %>% 
  pivot_longer(-c(Click, add_to_comp, Condition), names_to = "State", values_to = "Prob")
a$add_to_comp=factor(a$add_to_comp, levels = c("Baseline", "Click on covariate"))

# Create plot
alphaTR=a %>% 
  filter(Click<15) %>% 
  ggplot(aes(x=Click, y=Prob, fill=State))+
  geom_bar(stat="identity")+
  labs(y="State membership probability", fill="")+
  scale_x_continuous(breaks=c(1, 4, 7, 10, 14))+
  scale_fill_manual(values=c("#C1CDCD", "#454545"))+
  theme_light()+theme(legend.position = "bottom") +
  facet_grid(add_to_comp~Condition)
alphaTR
