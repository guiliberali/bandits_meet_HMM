#######################################################################################################################
# Author: 
# Purpose: Generates Figure 8 from page 38 of the paper
# 
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
#
#  Overview:
#     A) Load in Survey Data, the Gmatrix, and calculate the Gittens index
#     B) Plot the Gittens Index
#
#
#######################################################################################################################

getwd()
setwd('C:/Users/flori/OneDrive/Documents/Github/Replication_Morphing_HMM/study2/4. Figures & Tables')
##################
# A) Load in Survey Data, the Gmatrix, and calculate the Gittens index
##################
RAW_DATA = '../1. Raw Data'

### CLICKS DATA ----
processed_clicks_filtered <- read_csv(paste0( RAW_DATA, "/Application2_Production_study_unstacked_data.csv"))

## add conditions variable
clicks_data_with_conditions = processed_clicks_filtered %>%
  filter(aux0!=(-1)) %>%  ## removes clicks on instructions 
  group_by(user_id) %>% 
  mutate(first_morph=first(morph_id), 
         last_morph=last(morph_id) ) %>% 
  ungroup() %>% 
  mutate(condition=as.factor(case_when(first_morph==1 & last_morph==1 ~ 1, 
                                       first_morph==1 & last_morph==2 ~ 2,
                                       first_morph==2 & last_morph==1 ~ 3,
                                       TRUE ~ 4) ))

## load Gmatrix
Gmatrix <- read.table(DROPBOX_GMATRIX_LINK)  
m_Gmatrix=Gmatrix

# create dataframe wih alpha, beta, gittins - per click
df_cleanData_alpha_beta <- create_df_with_alpha_beta(clicks_data_with_conditions) 
df_cleanData_Gittins <- calc_G_indeces(df_cleanData_alpha_beta)

df_Gittins_at_assignment <- df_cleanData_Gittins %>% 
  filter(grepl("Morph recalculated using G index", aux3)) %>% 
  arrange(ms_timestamp)

df_users_Gittins=df_cleanData_Gittins %>% group_by(user_id) %>% 
  slice(n()) %>% 
  ungroup() %>% filter(cell==1) %>% 
  arrange(timestamp_unix) %>% 
  dplyr::select(-c(source_url, destination_url, context, ms_timestamp, cell)) %>% 
  mutate(Obs=1:n())

##################
# B)Plot the Gittens Index
##################

df_Gittins_plot<- df_users_Gittins %>%
  dplyr::select(starts_with("Gitt")) %>%
  mutate(Visitor=1:n()) %>% 
  pivot_longer(-Visitor, names_to="State_Morph", values_to="Index")

GI_plot=df_Gittins_plot%>% 
  ggplot(aes(x=Visitor, y=Index, color=State_Morph))+
  geom_line( size=1)+
  scale_color_manual(labels=c( "State 1 / Morph 1", "State 1 / Morph 2", "State 2 / Morph 1", "State 2 / Morph 2"), 
                     values=c( "darkred", "firebrick2", "darkblue", "blue1") ) +
  labs(y='Gittins Index', color="State/Morph")+
  theme_classic()+theme(legend.position = "bottom") 

GI_plot
