
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


## Load dataset with the variable coding at the link level
var_codes=read.csv(paste0(PATH_DATA, "/variable_codebook.csv"), header=T, sep=",", row.names=NULL)
var_codes=var_codes %>% 
  mutate(linkdepth=Hmm_dataset$link_depth2[match(URL, Hmm_dataset$auxvar.2)]) %>% 
  group_by(page_id) %>% 
  mutate(page_URL=rep(URL[1], length(URL))) %>% 
  na.omit()

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

## create interaction variables
CLdataset <- CLdataset %<>% 
  mutate(linkdepth_s2=linkdepth*State, c3_s2=c3*State, c15_s2=c15*State, c17_s2=c17*State)

# estimate the fit
fit.CLOMEGA<-clogit(y~linkdepth+c3+c15+c17+
                      linkdepth_s2+c3_s2+c15_s2+c17_s2+
                      strata(choice_id), method="exact",
                    data=CLdataset)      

# write out result 
write.csv(fit.CLOMEGA$coefficients, file=paste0(PATH_RESULTS, "/omega_2Sterminal.csv"))

## take page ids and link ids for the 10 most frequently clicked on pages (used in simulations)
c_ij=var_codes %>% 
  dplyr::select(page_id, link_id, linkdepth, c3, c15, c17) %>% 
  filter(page_id %in% c(1, 5, 2, 38, 3, 34,33, 4, 11,9)) %>% 
  mutate(page_id=as.factor(page_id))

# write out c_ij10.csv for simulation
c_ij[order(factor(c_ij$page_id, levels=c(1, 5, 2, 38, 3, 34,33, 4, 11,9))),]
write.csv(c_ij, file=paste0( "c_ij10.csv"))

