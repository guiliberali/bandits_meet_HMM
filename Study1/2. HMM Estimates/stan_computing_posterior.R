
set.seed(9000)

## Load required packages ----

rstan_options(auto_write = TRUE)

## hmm terminal model written in stan ------
## the code is used for the estimation of 2-, 3-, and 4-state terminal HMM models
## the 1-state model is a binary logit model for the emission probabilities
hmm_terminal=" 
functions {
vector normalize(vector x) {
return x / sum(x);
}
}

data {
int<lower=1> N;                   // number of observations
int<lower=1> IND;                  // number of indivivuals
int<lower=1> K;                   // number of hidden states  
int<lower=1> A; 
int<lower=1> L; 
int<lower=1> maxTN[N];            // N for the max number of obs. per individual; 
int<lower=1> T[N];               // vector for the time series per individual
int<lower=0, upper=1> lastT[N];  // dummy for last observation per individual
int<lower=1> ID[N];              // vector of individual index
int<lower=0, upper=1> y[N];      // observations for the choice eq.
row_vector[L] XL[N];            // covariates of the emission model, intercept not included
row_vector[A] XA[N];            // covariates of the transition model, intercept not included
int<lower=0> z[K];
vector[K] lastQij;
}


parameters {
ordered[K-1] mu[K-1];
ordered[K] beta;
vector[L] theta[K];
vector[A] rho[K-1];
simplex[K] pi0;
}

transformed parameters {
vector[K] unalpha_tk[N];
vector[K] unalpha_Tkn[IND];
vector[K] Q_ij[K];
for (n in 1:N)
{
real accumulator[K];

if(T[n]==1)
{
  for (j in 1:K) {
  unalpha_tk[n][j] = log(pi0[j]) + bernoulli_logit_lpmf(y[n]|beta[j]+XL[n]*theta[j]);
  }
}else{
  for (j in 1:K) { // j = current (t)
    for (i in 1:K) { // i = previous (t-1)
      if (i==K) Q_ij[i,j] = log(lastQij[j]);
         else Q_ij[i,j]=ordered_logistic_lpmf(z[j]|XA[n]*rho[i], mu[i]);
        accumulator[i] = unalpha_tk[n-1, i] + Q_ij[i, j]+bernoulli_logit_lpmf(y[n]|beta[j]+XL[n]*theta[j]);
  }
  unalpha_tk[n, j] = log_sum_exp(accumulator);
  }
}
if(lastT[n]==1) 
{unalpha_Tkn[ID[n]]=unalpha_tk[n];}

} // Forward
}

model {

for (k in 1:(K-1))
mu[k]~normal(0,1);

beta~normal(0,1);

for(k in 1:K)
theta[k]~normal(0, 1);

for(k in 1:(K-1))
rho[k]~normal(0,1);

for (i in 1:IND)
target += log_sum_exp(unalpha_Tkn[i]); // Update based only on last unalpha_tk

}

generated quantities{
vector<lower=0, upper=1>[K] phi_k[N];  // choice probability: pb of 1 or 0, given state, at each obs.
matrix[K, K] Qest_ij[N];               // transition probabilities
vector[K] unbeta_tk[N];
vector[K] ungamma_tk[N];
vector[K] alpha_tk[N];
vector[K] beta_tk[N];
vector[K] gamma_tk[N];

int<lower=1, upper=K> zstar_t[N];
real log_lik[N];
real logp_zstar_t;
int a_tk[N, K];                 // backpointer to the most likely previous state on the most probable path
real delta_tk[N, K];            // max prob for the seq up to t


// generate phi_k

for (n in 1:N)
{ // Forwards algorithm log p(z_t = j | x_{1:t})

for (j in 1:K) { // j = current (t)
  phi_k[n,j]=(y[n]==1 ? inv_logit(beta[j]+XL[n]*theta[j]) : (1-inv_logit(beta[j]+XL[n]*theta[j])));
  for (i in 1:K){
    if(i==K) Qest_ij[n,i,j] = log(lastQij[j]);
    else Qest_ij[n,i,j]=ordered_logistic_lpmf(z[j]|XA[n]*rho[i], mu[i]);
  }
}
}

// DRAW STATES

{ // Forward algortihm
for (n in 1:N)
alpha_tk[n] = softmax(unalpha_tk[n]);
} // Forward


{ // Backward algorithm log p(x_{t+1:T} | z_t = j)
real accumulator[K];

for (n in 1:N)
{
  if(lastT[n]==1)
  {for (j in 1:K)
  {unbeta_tk[n,j] = 1;}
  }
}

for (n in 1:N){
vector[K] Qij[K];
if(lastT[n]==1){
for (j in 1:K)
unbeta_tk[n,j] = 1;
}else{
int t;
t = maxTN[n] - (T[n]-1);

for(j in 1:K){
for (i in 1:K){
if (i==K) Qij[i,j] = log(lastQij[j]);
       else Qij[i,j]=ordered_logistic_lpmf(z[j]|XA[t]*rho[i], mu[i]);
}
}

for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// backwards t    + transition prob + local evIDence at t
accumulator[i] = unbeta_tk[t, i] + Qij[j, i]+bernoulli_logit_lpmf(y[t]|beta[i]+XL[t]*theta[i]);
}
unbeta_tk[t-1, j] = log_sum_exp(accumulator);
}
}
}

for (n in 1:N)
beta_tk[n] = softmax(unbeta_tk[n]);
} // Backward

{ // Forwards-backwards algorithm log p(z_t = j | x_{1:T})
for(n in 1:N) {
ungamma_tk[n] = alpha_tk[n] .* beta_tk[n];
}

for(n in 1:N)
gamma_tk[n] = normalize(ungamma_tk[n]);
} // Forwards-backwards

// Viterbi algorithm
// with final output from state k for time t

for (n in 1:N)
{
vector[K] Qij[K];
for (j in 1:K){
for (i in 1:K){

  if(i==K) Qij[i,j] = log(lastQij[j]);
  else Qij[i,j]=ordered_logistic_lpmf(z[j]|XA[n]*rho[i], mu[i]);
}
}

  if(T[n]==1)
  {for (j in 1:K)
  delta_tk[n, j] =  bernoulli_logit_lpmf(y[n]|beta[j]+XL[n]*theta[j]);// log(phi_k[j, x[n]]);
  }else{
  for (j in 1:K) { // j = current (t)
  delta_tk[n, j] = negative_infinity();
  for (i in 1:K) { // i = previous (t-1)
  real logp;
  logp = delta_tk[n-1, i] +Qij[i, j]+ bernoulli_logit_lpmf(y[n]|beta[j]+XL[n]*theta[j]);
  if (logp > delta_tk[n, j]) {
  a_tk[n, j] = i;
  delta_tk[n, j] = logp;
  }
  }
  }
  }
}

for (n in 1:N)
{
  if(lastT[n]==1)
  {logp_zstar_t = max(delta_tk[n]);
  
  for (j in 1:K)
  if (delta_tk[n, j] == logp_zstar_t)  zstar_t[n] = j;
  }
  log_lik[n]=logp_zstar_t;
  
}


for (n in 1:N) {
if(lastT[n]==1)
{zstar_t[n]=zstar_t[n];
}else{
int t;
t = maxTN[n] - T[n];
zstar_t[t] = a_tk[t + 1, zstar_t[t+1]];
log_lik[t]=delta_tk[t, zstar_t[t]];
}
}

}
"

m = stan_model(model_code = hmm_terminal)

## Import data ----

PATH_INPUT = 'input_hmm_estimates'
Data_final=read.csv(paste0(PATH_INPUT,"/AB_data_05_June_2018.csv"), header = T)

## link choice variables
## c3 - alumni, c15 - FAQ, c17 - programme info
## navigation variabe: link_depth2
## Intermediate outcomes: DB, CU, CV, RE
## create morph-specific data sets ----
hmm_data_morph1=Data_final %>% 
  select(user, timestamp_Clickstream, c3, c15, c17, link_depth2,morph, 
         Cumul_Outcome_DB, Cumul_Outcome_CV, Cumul_Outcome_CU, Cumul_Outcome_RE, Cumul_Int_Outcomes,
         Outcome_EA,  Outcome_DB, Outcome_CV, Outcome_CU, Outcome_RE, visit) %>%
  filter(morph==1) %>% 
  arrange(user, timestamp_Clickstream) %>% 
  na.omit() 

hmm_data_morph2=Data_final %>% 
  select(user, timestamp_Clickstream,c3, c15, c17, link_depth2,morph, 
         Cumul_Outcome_DB, Cumul_Outcome_CV, Cumul_Outcome_CU, Cumul_Outcome_RE, Cumul_Int_Outcomes,
         Outcome_EA,  Outcome_DB, Outcome_CV, Outcome_CU, Outcome_RE, visit) %>%
  filter(morph==2) %>% 
  arrange(user, timestamp_Clickstream) %>% 
  na.omit() 

####### SAMPLE SELECTION -----
## subset the data
## Morph 1 (abstract morph)
user_sample_m1=read_csv(paste0(PATH_INPUT, "/user_sample_m1.csv"))
hmm_data_morph1 %<>% filter(user %in% user_sample_m1$U_m1 ) 

## Morph 2 (concrete morph)
user_sample_m2=read_csv(paste0(PATH_INPUT, "/user_sample_m2.csv"))
hmm_data_morph2 %<>% filter(user %in% user_sample_m2$U_m2) 

############### RUN THE HMM ################################################
## run the HMM for the abstract morph (morph 1)


hmm_estimation <- function(hmm_data_morph){
  
  hmm_data=hmm_data_morph ## change to hmm_data_morph2 for the concrete morph (morph 2)
  hmm_data %<>% 
    group_by(user) %>% 
    mutate(T.length=row_number(), lastT=ifelse(row_number()==n(), 1,0)) %>% 
    ungroup()
  
  
  hmm_data %<>% 
    mutate(group_id = group_indices(., user), Obs=1:n())
  
  hmm_data %<>% 
    group_by(user) %>%
    mutate(maxTN = max(Obs)) %>% 
    ungroup()
  
  ## create logged lagged cumulative intermediate outcome variable
  hmm_data %<>%
    mutate(logCIO=log(Cumul_Int_Outcomes+1)) %>%
    group_by(user) %>% 
    mutate(lagCIO=lag(logCIO, 1)) %>% 
    mutate(lag1CIO=ifelse(is.na(lagCIO), 0,lagCIO)) %>% 
    ungroup()
  
  XL=hmm_data %>% select(lag1CIO)
  L=ncol(XL)  
  
  XA=as.matrix(hmm_data %>%
                 select(link_depth2, c3,c15, c17)%>% 
                 mutate(l_ld2= log(link_depth2))%>% 
                 select(l_ld2, c3, c15, c17))
  
  A=ncol(XA)
  K=2 ## change to K=3 for the 3-state model
  if (K==2) {
    z=c(1, 2)
    lastQij=c(0, 1)
  }else{ z=c(1, 2, 3) 
  lastQij=c(0, 0, 1)
  }
  
  maxTN= pull(hmm_data, maxTN)
  ID=pull(hmm_data, group_id)
  lastT=pull(hmm_data, lastT)
  Outcome_EA=pull(hmm_data, Outcome_EA)
  T.length=pull(hmm_data, T.length)
  
  ## Create stan data ----
  stan.data = list(
    maxTN=maxTN, 
    N=length(lastT),                   # number of observations
    IND=length(unique(ID)),             # number of visitors
    K=K,                               # number of hidden states      
    L=L,
    A=A,
    T=T.length,                       # vector for the time series per visitor
    lastT=lastT,                      # vector for the time series per visitor
    ID=ID,                            # vector of visitor index
    y= Outcome_EA,                   # observations for the choice eq.
    XL= XL, # predictors for the choice equation;
    XA=XA,                  # predictors for the state equation;
    z=z,
    lastQij=lastQij
  )
  
  ### Run the model ------
  hmm2S_Morph_terminal <- stan(model_code = hmm_terminal, data = stan.data,
                               verbose = TRUE, chains = 2,
                               seed = 9000, iter = 2000, warmup = 1500, control = list(adapt_delta=0.99, max_treedepth=15))
  
  return(hmm2S_Morph_terminal)
  
}

### Run the model ------
hmm2S_M1_terminal <-hmm_estimation(hmm_data_morph1)
hmm2S_M2_terminal <-hmm_estimation(hmm_data_morph2)

## save stan object for posterior inferences
save(hmm2S_M1_terminal, file="hmm2S_M1_terminal.RData")
save(hmm2S_M2_terminal, file="hmm2S_M2_terminal.RData")
