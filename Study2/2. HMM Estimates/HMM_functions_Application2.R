# Replication code for "Morphing for Consumer Dynamics: Bandits Meet HMM"
# Author: Gui Liberali and Alina Ferecatu 
# Date: October 2021
# Analysis for Application 2 - Calibration study
# Functions necessary for the HMM estimation and the Gittins plots in Application 2

## Hmm terminal model for Application 2, calibration data ----
hmm_terminal_geometric_ind_cov = "
data {
int<lower=1> N;                   // number of observations
int<lower=1> Nind;                // number of indivivuals
int<lower=1> K;                   // number of hidden states    
int<lower=1> A;                
int<lower=1> T[N];                // vector for the time series per individual
int<lower=0, upper=1> lastT[N];   // vector for the time series per individual
int<lower=1> IND[N];              // vector of individual index
row_vector[A] XA[N];            // transition prob variables, intercept NOT included
int<lower=0> z[K];
vector[K] lastQij;
vector[K] pi0;                 // everyone starts in the early morph in the experiment;
}

parameters {
vector[K] beta_bar;
vector[K-1] mu_bar[K-1];
vector<lower=0>[K] sigma_beta;
vector[K] alpha_beta[Nind];
vector[K-1] alpha_mu[Nind, K-1]; 
vector<lower=0>[K-1] sigma_mu[K-1]; 
vector[A] rho[K-1];
}

transformed parameters {
ordered[K] beta[Nind];
ordered[K-1] mu[Nind, K-1]; 

vector<lower=0, upper=1>[K] empb;
vector[K] unalpha_tk[N];
vector[K] unalpha_Tkn[Nind];

for (n in 1:Nind){

for (i in 1:(K-1)){ //for terminal;
  mu[n, i, 1]=mu_bar[i, 1]+sigma_mu[i, 1] * alpha_mu[n, i, 1];
  for (j in 2:(K-1)){ // for the ordered logit probs
  mu[n, i, j]=mu[n, i, (j-1)]+exp(mu_bar[i, j]+sigma_mu[i, j] * alpha_mu[n, i, j]);
  }
}

beta[n, 1]=beta_bar[1]+sigma_beta[1]*alpha_beta[n, 1];
  for (j in 2:K)
    beta[n, j]=beta[n, (j-1)]+exp(beta_bar[j]+sigma_beta[j]*alpha_beta[n, j]);
}


for (t in 1:N)
{ // Forwards algorithm log p(z_t = j | x_{1:t})
real accumulator[K];
// at T=1    
if(T[t]==1){
  for (j in 1:K) {
  unalpha_tk[t][j] = log(pi0[j]) + log(1-Phi_approx(beta[IND[t], j]));
  }
}else{
  if(lastT[t]==1) empb=Phi_approx(beta[IND[t]]);
    else empb=(1-Phi_approx(beta[IND[t]]));
    
  for (j in 1:K) { // j = current (t)
  for (i in 1:K) { // i = previous (t-1)
  if(i==K) accumulator[i] = unalpha_tk[t-1, i] + log(lastQij[j])+log(empb[j]);
       else accumulator[i] = unalpha_tk[t-1, i] + ordered_logistic_lpmf(z[j]|XA[t]*rho[i], mu[IND[t], i])+log(empb[j]);
  }
  unalpha_tk[t, j] = log_sum_exp(accumulator);
  }
}
if(lastT[t]==1) 
{unalpha_Tkn[IND[t]]=unalpha_tk[t];}

} // Forward
}

model {
beta_bar~normal(0,1);
sigma_beta~exponential(1);

for (j in 1:(K-1)) 
{mu_bar[j] ~ normal(0,1);
sigma_mu[j] ~ exponential(1);
}

for (n in 1:Nind){
alpha_beta[n] ~ normal(0,1);
for (i in 1:(K-1))
alpha_mu[n, i] ~ normal(0,1);
}

for(k in 1:(K-1))
rho[k]~normal(0,1);

for (n in 1:Nind)
target += log_sum_exp(unalpha_Tkn[n]); // Note: update based only on last unalpha_tk
}

generated quantities{
vector[K] alpha_tk[N];
real log_lik[Nind]; 

{ // Forward algortihm
for (t in 1:N)
alpha_tk[t] = softmax(unalpha_tk[t]);
} // Forward

for (n in 1:Nind)
log_lik[n]=log_sum_exp(unalpha_Tkn[n]);
}
"
