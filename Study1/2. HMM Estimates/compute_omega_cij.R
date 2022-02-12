#######################################################################################################################
# Author: 
# Purpose:computes the omega & and cij values, used for subsequent simulation
#  
# Note: 
#  -If you have not ran 'Configuration.R' and 'Functions.R' before running this file, make sure to do so
#  - After that, set working directory to Study1/2. HMM estimates
#
#  Overview:
#     A) Load data and perform initial transformations
#     B) Estimate Omega parameters
#     C) Create c_ij10.csv for simulation 
#
#######################################################################################################################

##############
# A) Load data and perform initial transformations
##############

PATH_RESULTS = 'Estimates for 2 states'
PATH_INPUT = 'input_hmm_estimates'

## Load clicks dataset ----
Hmm_dataset=read.csv(paste0(PATH_INPUT, "/Hmm_dataset_2Sterminal.csv"), header = T)


## Create state dummy
Hmm_dataset%<>% 
  group_by(user) %>% 
  mutate(state2_fixed= first(z_draw)) %>% 
  ## computes arrival state omega for the MAB_FIXED_STATE benchmark
  ungroup() %>% 
  mutate(State=ifelse(state2_fixed==1,-1, 1)) 


##* for 2-state terminal----
# Hmm_dataset%<>% 
#   mutate(State=ifelse(z_draw==1,-1, 1)) 

prepare_CLdataset <- function(Hmm_dataset, var_codes){
  
  ### join the clicks data with the variable description data
  CLdataset=Hmm_dataset %>% 
    filter(!grepl("#",auxvar.2)) %>% 
    filter(as.character(auxvar.1)!=as.character(auxvar.2)) %>% 
    dplyr::select(user, timestamp_Clickstream, auxvar.1, auxvar.2, State, link_depth2) %>% 
    mutate(y=rep(0, n())) %>% 
    left_join(var_codes,by=c("auxvar.1"="page_URL")) %>% 
    na.omit(page_id)%>% 
    dplyr::select(user, timestamp_Clickstream, auxvar.1, auxvar.2, page_id, URL, State, link_depth2, linkdepth, c3, c15, c17,y) 
  
  ## create an index for page times user interaction
  CLdataset %<>% 
    mutate(last=rep(0, n())) %>% 
    group_by(user,timestamp_Clickstream, page_id)%>%
    mutate(last=replace(last, length(last), 1)) 
  
  # define choice id
  choice_id=c(1,(cumsum(CLdataset$last)+1))
  choice_id=choice_id[1:nrow(CLdataset)]
  
  # add to CLdataset
  CLdataset=cbind(CLdataset,choice_id=choice_id)
  
  # get distance
  CLdataset %<>% 
    group_by(choice_id)%>%
    mutate(dist=levenshteinSim(as.character(auxvar.2), as.character(URL))) %>% 
    mutate(y=replace(y,dist==max(dist), 1))
  
  # gather double counted data
  double_counted_data=CLdataset %>% 
    group_by(choice_id) %>% 
    summarize(sumy=sum(y))
  double_counting=double_counted_data$choice_id[double_counted_data$sumy>=2]
  ## remove observations that are double counted
  CLdataset=filter(CLdataset, !(choice_id %in% double_counting ))
  
  
  ############# Arrival
  
  ## create interaction variables
  CLdataset <- CLdataset %<>% 
    mutate(linkdepth_s2=linkdepth*State, c3_s2=c3*State, c15_s2=c15*State, c17_s2=c17*State)
  
  return(CLdataset)
  
  
}


## Load dataset with the variable coding at the link level
var_codes=read.csv(paste0(PATH_INPUT, "/variable_codebook.csv"), header=T, sep=",", row.names=NULL)
var_codes=var_codes %>% 
  mutate(linkdepth=Hmm_dataset$link_depth2[match(URL, Hmm_dataset$auxvar.2)]) %>% 
  group_by(page_id) %>% 
  mutate(page_URL=rep(URL[1], length(URL))) %>% 
  na.omit()

CLdataset_arrival <- prepare_CLdataset(Hmm_dataset, var_codes)

Hmm_dataset%<>% 
   mutate(State=ifelse(z_draw==1,-1, 1)) 

CLdataset_actual <- prepare_CLdataset(Hmm_dataset, var_codes)


##############
# B) Estimate Omega parameters
##############

### arrival

# estimate the fit
fit.CLOMEGA_arrival<-clogit(y~linkdepth+c3+c15+c17+
                      linkdepth_s2+c3_s2+c15_s2+c17_s2+
                      strata(choice_id), method="exact",
                    data=CLdataset_arrival)      

# define object for writing out to csv
omega_2Sterminal_arrival <- data.frame(var = names(fit.CLOMEGA_arrival$coefficients), coef = fit.CLOMEGA_arrival$coefficients)
rownames(omega_2Sterminal_arrival) <- 1:8


#### Actual

# estimate the fit
fit.CLOMEGA_actual<-clogit(y~linkdepth+c3+c15+c17+
                      linkdepth_s2+c3_s2+c15_s2+c17_s2+
                      strata(choice_id), method="exact",
                    data=CLdataset_actual)      

# define object for writing out to csv
omega_2Sterminal_actual <- data.frame(var = names(fit.CLOMEGA_actual$coefficients), coef = fit.CLOMEGA_actual$coefficients)
rownames(omega_2Sterminal_actual) <- 1:8

# write out result 
write.csv(omega_2Sterminal_arrival, file=paste0(PATH_RESULTS, "/ArrivalState_omega_2Sterminal.csv"))
write.csv(omega_2Sterminal_actual, file=paste0(PATH_RESULTS, "/omega_2Sterminal.csv"))

##############
# C) Create c_ij10.csv for simulation 
##############

## take page ids and link ids for the 10 most frequently clicked on pages (used in simulations)
c_ij=var_codes %>% 
  dplyr::select(page_id, link_id, linkdepth, c3, c15, c17) %>% 
  filter(page_id %in% c(1, 5, 2, 38, 3, 34,33, 4, 11,9)) %>% 
  mutate(page_id=as.factor(page_id))

# write out c_ij10.csv for simulation
c_ij[order(factor(c_ij$page_id, levels=c(1, 5, 2, 38, 3, 34,33, 4, 11,9))),]
write.csv(c_ij, file=paste0( "c_ij10.csv"))

