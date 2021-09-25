##############  Morphing Consumer Dynamics - Repository of functions #####                                                  
## Terminal state - Compute individual specific lambda given previous morph
## Computes transition probabilities for previous state 1 to state K-1/ in our case we only have p11 and p12. 
##  Because if one starts in previous state 2, she stays there. So p21=0 and p22=1, this makes Q0_K=c(0,1)


############### Basic functions: logsumexp, softmax, inv_logit ##############  
logsumexp <- function(x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}
softmax <- function(x) {
  exp(x - logsumexp(x))
}

inv_logit <- function(logit) {
  exp(logit) / (1 + exp(logit))
}

############### Transition_Matrix ##############################
Transition_Matrix <- function(XA, rho, mu,Q0_K, DIM_NON_TERMINAL_STATES, TERMINAL)
{  
  if(TRANSITION_COVARIATES){
    
    ## initialize a (K-1)*K transition matrix per morph. 
    Aij=array( ,dim=c(DIM_NON_TERMINAL_STATES, TOT_STATES, TOT_MORPHS))
    for(m in 1: TOT_MORPHS)
    { 
      ## for each morph, compute the probability of staying in state 1
      if(TOT_STATES==2){
        if (TERMINAL) {Aij[,1,m]=1-inv_logit(-mu[,,m]+XA%*%rho[,,m])} else
        { Aij[,1,m]=1-inv_logit(-mu[,,m]+XA%*%t(rho[,,m]))}
        Aij[,2,m]=1-Aij[,1,m]
      }else{
        Aij[ ,1,m]=1-inv_logit(c(-mu[ ,1,m])+c(XA%*%t(rho[, ,m])))
        Aij[ ,2,m]=inv_logit(c(-mu[ ,1,m])+c(XA%*%t(rho[, ,m])))-inv_logit(c(-mu[ ,2,m])+c(XA%*%t(rho[, ,m])))
        Aij[ ,3,m]=inv_logit(c(-mu[ ,2,m])+c(XA%*%t(rho[,,m])))
      }
    }
    
    ## Adds the last row with zero-prob transitions out of terminal state
    if (TERMINAL) {Aij   <- array(c(rbind(Aij[,,1], Q0_K[,,1]), rbind(Aij[,,2], Q0_K[,,2])), dim=c(TOT_STATES, TOT_STATES, TOT_MORPHS))} 
  } else {
    ### NO COVARIATES IN THE TRANSITION MATRIX
    ## initialize a (K-1)*K transition matrix per morph. 
    Aij=array( ,dim=c(DIM_NON_TERMINAL_STATES, TOT_STATES, TOT_MORPHS))
    for(m in 1: TOT_MORPHS)
    { 
      ## for each morph, compute the probability of staying in state 1
      if(TOT_STATES==2){
        if (TERMINAL) {Aij[,1,m]=1-inv_logit(-mu[,,m]+1)} else
        { Aij[,1,m]=1-inv_logit(-mu[,,m]+1)}
        Aij[,2,m]=1-Aij[,1,m]
      }else{
        Aij[ ,1,m]=1-inv_logit(c(-mu[ ,1,m])+1)
        Aij[ ,2,m]=inv_logit(c(-mu[ ,1,m])+1)-inv_logit(c(-mu[ ,2,m])+1)
        Aij[ ,3,m]=inv_logit(c(-mu[ ,2,m])+1)
      }
    }
    
    ## Adds the last row with zero-prob transitions out of terminal state
    if (TERMINAL) {Aij   <- array(c(rbind(Aij[,,1], Q0_K[,,1]), rbind(Aij[,,2], Q0_K[,,2])), dim=c(TOT_STATES, TOT_STATES, TOT_MORPHS))} 
  }
  ## for each morph, computes the probability of moving to state 2
  ## in the main simulation, add the extra line to the transition matrix, for the terminal state
  return(Aij)
}

Transition_Matrix_study1 <- function(XA, rho, mu, DIM_STATES, TERMINAL)
{  
  ## initialize a (K-1)*K transition matrix per morph. 
  Aij=array( ,dim=c(DIM_STATES, TOT_STATES, TOT_MORPHS))
  for(m in 1: TOT_MORPHS)
  { 
    ## for each morph, compute the probability of staying in state 1
    if(TOT_STATES==2){
      if (TERMINAL) {Aij[,1,m]=1-inv_logit(-mu[,,m]+XA%*%rho[,,m])} else
      {Aij[,1,m]=1-inv_logit(-mu[,,m]+XA%*%t(rho[,,m]))}
      Aij[,2,m]=1-Aij[,1,m]
    }else{
      Aij[ ,1,m]=1-inv_logit(c(-mu[ ,1,m])+c(XA%*%t(rho[, ,m])))
      Aij[ ,2,m]=inv_logit(c(-mu[ ,1,m])+c(XA%*%t(rho[, ,m])))-inv_logit(c(-mu[ ,2,m])+c(XA%*%t(rho[, ,m])))
      Aij[ ,3,m]=inv_logit(c(-mu[ ,2,m])+c(XA%*%t(rho[,,m])))
    }
  }
  return(Aij)
}


############### Estimate_future_period_qs_general(Q0)  ############################
# qs=QS_ARRIVAL
# tr_matrix=Q0
estimate_future_period_qs_general<-function(tr_matrix)
{
  qs_expected_general    <-array(, dim=c(INIT_STATES, TOT_STATES, TOT_MORPHS, TOT_PERIODS))
  
  for (period in 1:TOT_PERIODS)
  {
    for(m in 1:TOT_MORPHS)
      qs_expected_general[ , , m,period] <- (tr_matrix[, ,m]%^% (period))
  }
  
  return(qs_expected_general)
}

############### Estimate_future_period_qs(Q0, QS_ARRIVAL) ############
estimate_future_period_qs<-function(tr_matrix, qs)
{
  qs_expected    <-array( , dim=c(INIT_STATES, TOT_STATES, TOT_MORPHS, TOT_PERIODS))
  
  # morph 1
  qs_expected[,,1,1]<-matrix(rep(qs, 2), nrow=2,byrow=T)
  
  # morph 2
  qs_expected[,,2,1]<-matrix(rep(qs, 2), nrow=2,byrow=T)
  
  for (period in 2:TOT_PERIODS)
  {
    period_minus_1 <- period-1
    
    ## morph 1
    ## initial stage 1
    qs_expected[1,1,1,period] <-tr_matrix[1,1,1]*qs_expected[1,1,1,period_minus_1]
    qs_expected[1,2,1,period] <- 1-qs_expected[1,1,1,period]
    
    ## initial stage 2
    qs_expected[2,1,1,period] <- tr_matrix[2,1,1]*qs_expected[2,1, 1,period_minus_1]
    qs_expected[2,2,1,period] <- 1-qs_expected[2,1,1,period]
    
    ## morph 2
    ## initial stage 1
    qs_expected[1,1,2,period] <- tr_matrix[1,1,2]*qs_expected[1,1,2,period_minus_1]
    qs_expected[1,2,2,period] <- 1-qs_expected[1,1,2,period]
    
    ## initial stage 2
    qs_expected[2,1,2,period] <- tr_matrix[2,1,2]*qs_expected[2,1, 2,period_minus_1]
    qs_expected[2,2,2,period] <- 1-qs_expected[2,1,2,period]
  }
  return(qs_expected)
}
#click=1
############### Within-person  DP #################
# version control
#  - Gui added the m- and s- specific bouncing probability instead of the discount rate
#  - Note that Alina suggested separating number of clicks from number of periods locally in this function only.
#    Gui rejected that change because (1) in this paper we are not doing this separation (2) if and when we decide to do this,
#    we will have to do this globally, not in a single function. Even if our initial intuition is that this is the only place 
#    that this matters, this needs to be thought trhoughly instead of "ad-hoc"-ally :) otherwise we simply can create new 
#    errors for no reason

# with seven periods ahead ----
DP = function(K, qs_future, qs, G_current_sm_POST, G_current_sm_PRE, p_bounce) # , delta_dp)   
{
  V_immediate    <- V_continuation <- array(0, c(TOT_STATES, TOT_MORPHS, TOT_PERIODS) )  
  V_total_cond   <- array(0,c(TOT_STATES, TOT_MORPHS, TOT_PERIODS+1) )  
  V_total        <- array(0,c(TOT_MORPHS, TOT_PERIODS) ) 
  V_max          <- rep(0,TOT_PERIODS)  
  best_morph     <- rep(0,TOT_PERIODS)
  V_total_cond[,,K_FULL+1] <-  0    # F_FULL PERIOD
  
  for (period in K_FULL:K)
  {			
    next_period <- period+1
    for (morph in 1: TOT_MORPHS) 
    {
      for (stage in 1:TOT_STATES) 
      { 
        # exponent on p_bounce should be (period - K), not (next_period-k) 
        pmf_bounce_ms=((1-p_bounce[stage, morph])^(period-K))*p_bounce[stage, morph]
        
        # IMMEDIATE REWARD  
        if (period < 7) { V_immediate[stage,morph, period]  <- crossprod(qs_future[stage, ,morph ,(next_period-K)], G_current_sm_PRE[,morph])  }
        if (period >=7) { V_immediate[stage,morph, period]  <- crossprod(qs_future[stage, ,morph ,(next_period-K)], G_current_sm_POST[,morph]) }
        # CONTINUATION REWARD    
        V_continuation[stage, morph,period] <-pmf_bounce_ms * V_total_cond[ stage, morph, next_period]   
        
        # CONDITIONAL TOTAL REWARD  
        V_total_cond[ stage, morph, period] <- V_immediate[ stage, morph, period]  + V_continuation[stage, morph, period]          
        
      } # close loop over STATES
      
      # UNCONDITIONAL TOTAL 
      V_total[morph, period] <-crossprod(qs, V_total_cond[ ,morph, period] )
      
    } # close loop over morphs
    
    #  VMAX: BEST MORPH GIVEN PREVIOUS     
    V_max[period] <-	max(V_total[, period])   
    best_morph[period] <- which.is.max(V_total[, period]) 
  }  #  close loop over period
  
  return(best_morph)
} 


## with four periods ahead ----
DP_4periods = function(K, qs_future, qs, G_current_sm, p_bounce) # , delta_dp)   
{
  V_immediate    <- V_continuation <- array(0, c(TOT_STATES, TOT_MORPHS, TOT_CONSIDERED_PERIODS) )  
  V_total_cond   <- array(0,c(TOT_STATES, TOT_MORPHS, TOT_CONSIDERED_PERIODS+1) )  
  V_total        <- array(0,c(TOT_MORPHS, TOT_CONSIDERED_PERIODS) ) 
  V_max          <- rep(0,TOT_CONSIDERED_PERIODS)  
  best_morph     <- rep(0,TOT_CONSIDERED_PERIODS)
  V_total_cond[,,TOT_CONSIDERED_PERIODS+1] <-  0    # F_FULL PERIOD
  
  Future_period=c(7, 9, 12, 14)
  
  for (i in TOT_CONSIDERED_PERIODS:1)
  {			
    period=Future_period[i]
    next_period <- period+1
    
    for (morph in 1:TOT_MORPHS) 
    {
      for (stage in 1:TOT_STATES) 
      { 
        # exponent on p_bounce should be (period - K), not (next_period-k) 
        pmf_bounce_ms=((1-p_bounce[stage, morph])^(period-K))*p_bounce[stage, morph]
        
        # IMMEDIATE REWARD   
        V_immediate[stage,morph, i]  <- crossprod(qs_future[stage, ,morph ,(next_period-K)], G_current_sm[,morph])  
        
        # CONTINUATION REWARD    
        V_continuation[stage, morph,i] <-pmf_bounce_ms * V_total_cond[ stage, morph, (i+1)]   
        
        # CONDITIONAL TOTAL REWARD  
        V_total_cond[ stage, morph, i] <- V_immediate[ stage, morph, i]  + V_continuation[stage, morph, i]          
        
      } # close loop over STATES
      
      # UNCONDITIONAL TOTAL 
      V_total[morph, i] <-crossprod(qs, V_total_cond[ ,morph, i] )
      
    } # close loop over morphs
    
    #  VMAX: BEST MORPH GIVEN PREVIOUS     
    V_max[i] <-	max(V_total[, i])   
    best_morph[i] <- which.is.max(V_total[, i]) 
  }  #  close loop over period
  
  return(best_morph)
} 

## with four periods ahead ----
DP_webstore_implementation = function(G_current_sm) # , delta_dp)   
{
  
  Value_M1=ifelse(XA==0, (0.23564*G_current_sm[1, 1] + 0.80561 * G_current_sm[2, 1]),
                  ( 0.4116289*G_current_sm[1, 1]+ 0.63883*G_current_sm[2, 1]))
  Value_M2=ifelse(XA==0, (0.099 *G_current_sm[1, 2] + 0.96086 * G_current_sm[2, 2]),
                  ( 0.23884 * G_current_sm[1, 2] + 0.82523*G_current_sm[2, 2])) 
  
  V_total        <- c(Value_M1, Value_M2)
  
  best_morph <- which.is.max(V_total)
  
  return(best_morph)
} 

############### Compute_click_prob_conditionalLogit ################### 
# computes clicks probabilities per morph per stage following the conditional logit model
# set.seed(123456789)
Compute_click_prob_conditionalLogit<-function(K_FULL)
{
  
  ## computes the denominator per page_id
  den<-matrix(, nrow=K_FULL, ncol=TOT_STATES)
  for (k in 1:K_FULL)   {for(s in 1:TOT_STATES)
  {den[k, s] = sum(exp(Ckjn[Ckjn[,1]==k,2:ncol(Ckjn)] %*% omega %*% all_s[s,]))   }    }
  
  ## computes the probability per link_id
  temp <- array(,c(TOT_STATES,TOT_LINKS, TOT_MORPHS))
  for (m in 1:TOT_MORPHS)  {for(s in 1:TOT_STATES) {for (j in 1:TOT_LINKS) 
  { temp[s, j, m] = exp(Ckjn[j,2:ncol(Ckjn)] %*% omega %*% all_s[s,]) / den[Ckjn[j,1], s] }}  }
  
  ## appends the page_id, to be used in the Y_vec
  temp1 <- array(,c((TOT_STATES+1),TOT_LINKS, TOT_MORPHS))
  for (m in 1:TOT_MORPHS) {temp1[,,m]=rbind(Ckjn[,1], temp[,,m])}
  return(temp1)
}

############### Compute_click_prob_restrictedMultinomialLogit ################### 
## computes clicks probabilities per morph per stage
## following the restricted multinomial logit model
## Input - variables are lagged
## c1, c5, c7: characteristics of the page previously chosen
## c8-c12: characteristics of the source page 
## (the source page = the page participants were on when deciding to choose the previously chosen page)
Compute_click_prob_restrictedMultinomialLogit <- function(D_in, Omega_in)
{  
  prob = array( , dim=c(nrow(index_tr_matrix), nrow(index_tr_matrix), TOT_STATES, D_in))
  ## This prob says: given that participant chose page "previous_pagechosen" starting from "previous_source" page in the previous click
  ## and they're in state S, they will now have the choice prob for current page d in {1....D_in}.
  ## prob has dimensions [previous_source (D_in-1) * previous_pagechosen (TOT_PAGES-1) * state S * current choice D_in ]
  ## previous_source has dimension (D_in-1) because they will not start on the last page 6 (that's the WTP page)
  ## In the simulation, keep track of previous_source (or index_previous_source) and previous_pagechosen (or index_previous_pagechosen)
  ## to recover choice prob for the current page, which is prob[index_previous_source, index_previous_pagechosen, state, ],
  ## and draw a click give this D_in-dimensional probability vector
  
  ## source page at the previous click
  for (index_previous_source in 1:(D_in)) # -1 removed   by Gui to handle c_kjm without wtp
  {
    ## destination page at the previous click  
    for (index_previous_pagechosen in 1:(D_in))  # -1 removed   by Gui to handle c_kjm without wtp
    {
      Ckjn_previous_source=as.matrix(Ckjn[index_previous_source, 10:14])
      Ckjn_previous_pagechosen=c(as.matrix(Ckjn[index_previous_pagechosen, c(3, 7, 9)]), Ckjn_previous_source)
      omega_restr = Omega_in
      for (d in 1:D_in)
      {     
        if (index_tr_matrix[index_previous_pagechosen, d] == 0) {omega_restr[1, d]= -100} 
        # ensures X[,1] stays positive, to get -100
        for (state in 1:TOT_STATES)
        {
          temp = rep(0, D_in)
          for(d in 1:D_in) {temp[d]=(t(c(Ckjn_previous_pagechosen %*% t(all_s[state, ]))) %*% omega_restr[, d] )   }
          prob[index_previous_source, index_previous_pagechosen, state, ] = softmax(temp)
        } # close state loop 
        
      }  # close d loop 
      
    }  # close index_previous_pagechosen
  }    # close index_previous_pagechosen loop
  return(prob)
}

#####################
# Used for study 1

Compute_click_prob<-function(K_FULL)
{
  ## computes the denominator per page_id
  den<-matrix(, nrow=K_FULL, ncol=TOT_STATES)
  for (k in 1:K_FULL)   {for(s in 1:TOT_STATES)
  {den[k, s] = sum(exp(Ckjn[Ckjn[,1]==k,2:ncol(Ckjn)] %*% omega %*% all_s[s,]))   }    }
  
  ## computes the probability per link_id
  temp <- array(,c(TOT_STATES,TOT_LINKS, TOT_MORPHS))
  for (m in 1:TOT_MORPHS)  {for(s in 1:TOT_STATES) {for (j in 1:TOT_LINKS) 
  { temp[s, j, m] = exp(Ckjn[j,2:ncol(Ckjn)] %*% omega %*% all_s[s,]) / den[Ckjn[j,1], s] }}  }
  
  ## appends the page_id, to be used in the Y_vec
  temp1 <- array(,c((TOT_STATES+1),TOT_LINKS, TOT_MORPHS))
  for (m in 1:TOT_MORPHS) {temp1[,,m]=rbind(Ckjn[,1], temp[,,m])}
  return(temp1)
}


############### G_interpolation_fast  ####### 
# retrives the gittins index via interpolation, because we don't have integer alphas and betas 
G_interpolation_fast<-function(alpha,beta)
{return(mean(Gmatrix[round(beta), round(alpha)], Gmatrix[round(beta), trunc(alpha)],  Gmatrix[trunc(beta), round(alpha)], Gmatrix[trunc(beta), trunc(alpha)]))}

############### Draw_state ####################################
# Choose true_state - multinomial with prob pi0 
#set.seed(123456789)
Draw_state<-function(STATES, prob) { return(which.is.max(rmultinom(STATES, 1, prob)))}

############### Get_starting_morph  ###################################
# Provides the initial morph as specified by INITIAL_MORPH
Get_starting_morph =  function (G_current_sm) {
  current_m <- rep(0,TOT_MORPHS)
  for (m in 1:TOT_MORPHS) {current_m[m] <-pi0  %*% G_current_sm[,m] }
  return (which.is.max(current_m) ) 
}

############### Bayesian_updater   #####################################
## bayesian updater for q_r, with br_updated as prior
Bayesian_updater<-function(K, M, navigation_prior, omega, Ckjn,all_s, Y_mat )
{
  num=num_jointclicks=matrix(rep(0, K*TOT_STATES), ncol=TOT_STATES)
  den_jointclicks=rep(0, K)
  for (k in 1:K)
  {
    for (r in 1:TOT_STATES)
    {
      num[k,r]=prod((exp(Ckjn[Ckjn[,1]==k,2:ncol(Ckjn)]%*%omega%*%all_s[r,])/sum(exp(Ckjn[Ckjn[,1]==k,2:ncol(Ckjn)]%*%omega%*%all_s[r,])))^(Y_mat[Y_mat[,1]==k,2]))
    }
    if(k==1){
      num_jointclicks[k,]=num[k,]*navigation_prior
      den_jointclicks[k]=sum(num_jointclicks[k,])
    }else{
      num_jointclicks[k,]=(apply(num[1:k,], 2, prod))*navigation_prior
      den_jointclicks[k]=sum(num_jointclicks[k,])
    }
  }
  return(num_jointclicks/den_jointclicks)
}
############### Bayesian updating - Montoya et al 2010 ###############
##  Value:  prob of state s' at t+1 given state s at t, organized in a K*K matrix, where on the rows you have the current state, and on the columns you have the previous states.  
##  Intermediate steps:
##    numerator: probability of ending up in state s' given all possible starting points s
##    denominator: sum over all probabilities of being in a state s', starting from state s and summed acrross all previous states s
# Tr_matrix is [TOT_STATES, TOT_STATES, TOT_MORPHS]
# qs_ is [K_FULL+1, TOT_STATES]
# psm_ is [TOT_MORPHS, TOT_STATES]
## Next 5 lines are for testing purposes - comment out in production
# k=click
# tr_matrix=Q0
# state_probs=qs_HMM
# success_prob_per_state = psm_true
# morph=best_morph_next_t
# CODED BY: ALINA 
#Bayesian_Thm_updater<-function(click, lambda_, qs_, morph, psm_ )
#{
#  numerator   <-  matrix(, TOT_STATES, TOT_STATES)
#  for (current_s in 1:TOT_STATES) {numerator[current_s,] <- lambda_[ , current_s, morph]*qs_[click, ] * psm_[morph, current_s]  }
#  denominator <- rowSums(numerator)
#  return(numerator/denominator)
#}  

############### Bayesian_Thm_updater_element ###############  
Bayesian_Thm_updater_element_cond<-function(lambda_, P_st, morph, geometric_emission_probs)
{
  p_deltas<- as.data.frame(geometric_emission_probs[geometric_emission_probs$Tr_cov_model==TRANSITION_COVARIATES,])
  
  # Notes
  #  1. Now psm depends on the previous morph seen. This is given empirically in the RCT.
  #  2. All that matters is the first morph seen when selecting the morph at t=7 hence 
  #     temporarily we do p_delta[1,] and p_delta[2,] for now
  #p_delta <- matrix(0,nr=TOT_MORPHS, nc=TOT_STATES)
  #if (first_morph ==1) 
  #  {
  #   p_delta[1,] <- p_deltas$p_purchase_cond1[click=7] # cond 1 is m1m1
  #   p_delta[2,] <- p_deltas$p_purchase_cond2[click=7] # cond 2 is m1m2
  #   }
  #if (first_morph ==2) 
  #  {
  #   p_delta[1,] <- p_deltas$p_purchase_cond3[click=7] # cond 1 is m1m1
  #   p_delta[2,] <- p_deltas$p_purchase_cond4[click=7] # cond 2 is m1m2
  #   }
  p_delta=matrix((1-p_deltas$prob), nrow=TOT_MORPHS, ncol=TOT_STATES, byrow = T)
  
  numerator   <-  matrix(, TOT_STATES) 
  if (TOT_STATES ==2) 
  {
    # For clarity I am doing the expectation by hand, explicitly. Later I can generalize it for s =1..TOT_STATES 
    numerator[s=1] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=1, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=1, morph] * p_delta[morph, s_prime=2]      
    
    numerator[s=2] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=2, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=2, morph] * p_delta[morph, s_prime=2]      
    
  } else
  { # TOT_STATES=1
    # For clarity I am doing the expectation by hand, explicitly. Later I can generalize it for s =1.. TOT_STATES 
    numerator[s=1] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=1, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=1, morph] * p_delta[morph, s_prime=2] +
      P_st[s_prime=3] * lambda_[s_prime=3, s=1, morph] * p_delta[morph, s_prime=3] 
    
    
    numerator[s=2] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=2, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=2, morph] * p_delta[morph, s_prime=2] +     
      P_st[s_prime=3] * lambda_[s_prime=3, s=2, morph] * p_delta[morph, s_prime=3]      
    
    numerator[s=3] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=3, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=3, morph] * p_delta[morph, s_prime=2] +     
      P_st[s_prime=3] * lambda_[s_prime=3, s=3, morph] * p_delta[morph, s_prime=3]      
  } 
  return( numerator/sum(numerator)  )
}  



############### Bayesian_Thm_updater_element ###############  
Bayesian_Thm_updater_element<-function(lambda_, P_st, morph, p_delta )
{
  ## for testing
  # lambda_=lambda
  # P_st=qs_HMM[click,]
  # morph=best_morph
  
  numerator   <-  matrix(, TOT_STATES) 
  if (TOT_STATES ==2) 
  {
    # For clarity I am doing the expectation by hand, explicitly. Later I can generalize it for s =1..TOT_STATES 
    numerator[s=1] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=1, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=1, morph] * p_delta[morph, s_prime=2]      
    
    numerator[s=2] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=2, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=2, morph] * p_delta[morph, s_prime=2]      
    
  } else
  { # TOT_STATES=1
    # For clarity I am doing the expectation by hand, explicitly. Later I can generalize it for s =1.. TOT_STATES 
    numerator[s=1] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=1, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=1, morph] * p_delta[morph, s_prime=2] +
      P_st[s_prime=3] * lambda_[s_prime=3, s=1, morph] * p_delta[morph, s_prime=3] 
    
    
    numerator[s=2] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=2, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=2, morph] * p_delta[morph, s_prime=2] +     
      P_st[s_prime=3] * lambda_[s_prime=3, s=2, morph] * p_delta[morph, s_prime=3]      
    
    numerator[s=3] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=3, morph] * p_delta[morph, s_prime=1] + 
      P_st[s_prime=2] * lambda_[s_prime=2, s=3, morph] * p_delta[morph, s_prime=2] +     
      P_st[s_prime=3] * lambda_[s_prime=3, s=3, morph] * p_delta[morph, s_prime=3]      
  } 
  return( numerator/sum(numerator)  )
}  

# state_updater=State_bayesian_updater(click, Q0, qs_HMM, best_morph_next_t, psm_true)
# rowSums(state_updater)
################
# Functions to extract alpha + beta values, calculate gittins index
################

###########################
#text_alpha_beta_to_number
# 
# Input: 
#   text_alpha_beta; character, text that contains alpha and beta values
#
# Output:
#   vector of alpha and beta values (size 8)
############################

text_alpha_beta_to_number <- function(text_alpha_beta){
  
  # grab text with alpha values, split it, and then grab numeric values
  alpha_text <- str_match(text_alpha_beta, "Alpha\\s*(.*?)\\s*]]")[2]
  alpha_split <- str_split(alpha_text,';', simplify = TRUE)
  alpha_m_numbers <- as.numeric(gsub("\\[|\\]", "", alpha_split))
  
  
  # grab text with beta values, split it, and then grab numeric values
  beta_text <- str_match(text_alpha_beta, "Beta\\s*(.*?)\\s*]]")[2]
  beta_split <- str_split(beta_text,';', simplify = TRUE)
  beta_m_numbers <- as.numeric(gsub("\\[|\\]", "", beta_split))
  
  
  return(c(alpha_m_numbers, beta_m_numbers))
  
}

##############################
#create_df_with_alpha_beta
#
#  input: 
#     df_cleanData: dataframe, contains aux2 variable
#
#  output: 
#    returns the dataframe with alpha and beta values added
#############################

create_df_with_alpha_beta <- function(df_cleanData){
  
  # get number of observations, create vectors of this size
  n_obs =  length(df_cleanData$aux2)
  
  # one vector for each value
  v_alpha_s1_m1 <- rep(NA,n_obs)
  v_alpha_s2_m1 <- rep(NA,n_obs)
  v_alpha_s1_m2 <- rep(NA,n_obs)
  v_alpha_s2_m2 <- rep(NA,n_obs)
  v_beta_s1_m1 <- rep(NA,n_obs)
  v_beta_s2_m1 <- rep(NA,n_obs)
  v_beta_s1_m2 <- rep(NA,n_obs)
  v_beta_s2_m2 <- rep(NA,n_obs)
  
  # loop through all values
  for(i in seq(1, n_obs)){
    
    # get the text
    text = df_cleanData$aux2[i]
    
    # get alpha and beta values
    alpha_beta_numbers <- text_alpha_beta_to_number(text)
    
    # add values to vectors
    v_alpha_s1_m1[i] <- alpha_beta_numbers[1]
    v_alpha_s2_m1[i] <- alpha_beta_numbers[2]
    v_alpha_s1_m2[i] <- alpha_beta_numbers[3]
    v_alpha_s2_m2[i] <- alpha_beta_numbers[4]
    v_beta_s1_m1[i] <- alpha_beta_numbers[5]
    v_beta_s2_m1[i] <- alpha_beta_numbers[6]
    v_beta_s1_m2[i] <- alpha_beta_numbers[7]
    v_beta_s2_m2[i] <- alpha_beta_numbers[8]
    
    
  }
  
  # add vectors to dataframe and return
  df_cleanData$alpha_s1_m1 <- v_alpha_s1_m1
  df_cleanData$alpha_s2_m1 <- v_alpha_s2_m1
  df_cleanData$alpha_s1_m2 <- v_alpha_s1_m2
  df_cleanData$alpha_s2_m2 <- v_alpha_s2_m2
  df_cleanData$beta_s1_m1 <- v_beta_s1_m1
  df_cleanData$beta_s2_m1 <- v_beta_s2_m1
  df_cleanData$beta_s1_m2 <- v_beta_s1_m2
  df_cleanData$beta_s2_m2 <- v_beta_s2_m2
  
  return(df_cleanData)
}



##############################
#G_interpolation
#
#  input: 
#     alpha; float, alpha value
#     beta; float, beta value
#     Gmatrix; matrix object, with gittins index values 
#
#  output: 
#    returns the Gittins index score
# #############################
G_interpolation<-function(alpha,beta, Gmatrix){

  # trunc and round the alpha and beta
  trunc_alpha = trunc(alpha)
  round_alpha = round(alpha)
  trunc_beta = trunc(beta)
  round_beta = round(beta)

  # ensure there is a minimum index of 1,1
  if(trunc_alpha == 0){
    trunc_alpha <- 1
  }
  if(round_alpha == 0){
    round_alpha <- 1
  }
  if(trunc_beta == 0){
    trunc_beta <- 1
  }

  # take mean of round/truncated indeces
  return(mean(Gmatrix[round_beta, round_alpha],
              Gmatrix[round_beta,trunc_alpha],
              Gmatrix[trunc_beta, round_alpha],
              Gmatrix[trunc_beta, trunc_alpha]))

}
##############################
#calc_G_indeces
#
#  input: 
#     df_cleanData_alpha_beta: dataframe, contains alpha and beta columns
#
#  output: 
#    returns the dataframe with gittins index value added
#############################

calc_G_indeces <- function(df_cleanData_alpha_beta){
  
  # save gittins index scores in vectors here
  n_rows <- nrow(df_cleanData_alpha_beta)
  vec_Gittins_s1_m1 <- rep(NA, n_rows)
  vec_Gittins_s2_m1 <- rep(NA, n_rows)
  vec_Gittins_s1_m2 <- rep(NA, n_rows)
  vec_Gittins_s2_m2 <- rep(NA, n_rows)
  
  
  # start for loop
  for(i in seq(1, n_rows)){
    
    # get alpha and beta for s1, m1 
    alpha_s1_m1 = df_cleanData_alpha_beta$alpha_s1_m1
    beta_s1_m1 = df_cleanData_alpha_beta$beta_s1_m1
    
    # get alpha and beta for s2, m1
    alpha_s2_m1 = df_cleanData_alpha_beta$alpha_s2_m1
    beta_s2_m1 = df_cleanData_alpha_beta$beta_s2_m1
    
    # get alpha and beta for s1, m2
    alpha_s1_m2 = df_cleanData_alpha_beta$alpha_s1_m2
    beta_s1_m2 = df_cleanData_alpha_beta$beta_s1_m2
    
    # get alpha and beta for s2, m2 
    alpha_s2_m2 = df_cleanData_alpha_beta$alpha_s2_m2
    beta_s2_m2 = df_cleanData_alpha_beta$beta_s2_m2
    
    # calculate the gittins per morph and state
    Gittins_s1_m1 <- G_interpolation(alpha_s1_m1[i], beta_s1_m1[i], m_Gmatrix)
    Gittins_s2_m1 <- G_interpolation(alpha_s2_m1[i], beta_s2_m1[i], m_Gmatrix)
    Gittins_s1_m2 <- G_interpolation(alpha_s1_m2[i], beta_s1_m2[i], m_Gmatrix)
    Gittins_s2_m2 <- G_interpolation(alpha_s2_m2[i], beta_s2_m2[i], m_Gmatrix)
    
    # add gittins to vector
    vec_Gittins_s1_m1[i] <- Gittins_s1_m1
    vec_Gittins_s2_m1[i] <- Gittins_s2_m1
    vec_Gittins_s1_m2[i] <- Gittins_s1_m2
    vec_Gittins_s2_m2[i] <- Gittins_s2_m2
    
  }
  
  
  # add to datafarme and return
  df_cleanData_alpha_beta$Gittins_s1_m1 <- vec_Gittins_s1_m1
  df_cleanData_alpha_beta$Gittins_s2_m1 <- vec_Gittins_s2_m1
  df_cleanData_alpha_beta$Gittins_s1_m2 <- vec_Gittins_s1_m2
  df_cleanData_alpha_beta$Gittins_s2_m2 <- vec_Gittins_s2_m2
  
  return(df_cleanData_alpha_beta)
}
