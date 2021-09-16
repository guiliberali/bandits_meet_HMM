

### Replication File for Table, page 40 of the paper 
### Run the '2. HMM Estimates/HMM_estimates.R' File before running any code here.

## Table 6
stacked_dataset_with_Tlength %>% 
  group_by(stayed_7clicks, condition_reassigned) %>% 
  summarise(n())