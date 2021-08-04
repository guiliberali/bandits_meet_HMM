#################################################################################################################
#  Morphing Consumer Dynamics - Repository of functions                                                        
#  March 16, 2019
#
#   Changes made on March 16, 2019 to the version used until then
#       - Relabeled G_current to G_current_sm in the entire file, all functions
## Terminal state - Compute individual specific lambda given previous morph ######
## computes transition probabilities for previous state 1 to state K-1/ in our case we only have p11 and p12. Because if one starts in previous state 2, she stays there. So p21=0 and p22=1, this makes Q0_K=c(0,1)
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))

############################################################################
Transition_Matrix <- function(XA, rho, mu, DIM_STATES, TERMINAL)
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
  ## for each morph, computes the probability of moving to state 2
  ## in the main simulation, add the extra line to the transition matrix, for the terminal state
  return(Aij)
}

############################################################################
#qs=QS_ARRIVAL
#tr_matrix=Q0
#estimate_future_period_qs_general(Q0)
estimate_future_period_qs_general<-function(tr_matrix)
{
  qs_expected_general    <-array(, dim=c(INIT_STATES, TOT_STATES, TOT_MORPHS, TOT_PERIODS))
  
  for (period in 1:TOT_PERIODS)
  {
    for(m in 1:TOT_MORPHS)
      qs_expected_general[ , , m,period] <- (tr_matrix[, ,m] %^% (period))
  }
  
  return(qs_expected_general)
}

############################################################################
#estimate_future_period_qs(Q0, QS_ARRIVAL)
estimate_future_period_qs<-function(tr_matrix, qs)
{
  qs_expected    <-array(, dim=c(INIT_STATES, TOT_STATES, TOT_MORPHS, TOT_PERIODS))
  
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
## WITHIN-PERSON  DP
# version control
#  - Gui added the m- and s- specific bouncing probability instead of the discount rate
#  - Note that Alina suggested separating number of clicks from number of periods locally in this function only.
#    Gui rejected that change because (1) in this paper we are not doing this separation (2) if and when we decide to do this,
#    we will have to do this globally, not in a single function. Even if our initial intuition is that this is the only place 
#    that this matters, this needs to be thought trhoughly instead of "ad-hoc"-ally :) otherwise we simply can create new 
#    errors for no reason
#

DP = function(K, qs_future, qs, G_current_sm, p_bounce) # , delta_dp)   
{
  # for testing
  #K=click
  #qs_future=qs_HMM_sys_future
  #qs=qs_HMM[click,]
  #p_bounce=p_bounce_click
  
  V_immediate    <- V_continuation <- array(0,c(TOT_STATES, TOT_MORPHS, TOT_PERIODS) )  
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
       # cdf_bounce_ms   <-  p_bounce_mt[period, stage, morph ]
        # pmf_bounce_ms   <-  cdf_bounce_ms * (1-cdf_bounce_ms)^(period) 
        pmf_bounce_ms=((1-p_bounce[stage, morph])^(period-K))*p_bounce[stage, morph]
        
        # IMMEDIATE REWARD   
        V_immediate[stage,morph, period]  <- crossprod(qs_future[stage, ,morph ,(next_period-K)], G_current_sm[,morph])  
      
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
 
##set.seed(123456789)
## computes clicks probabilities per morph per stage
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

############################################################################
## retrives the gittins index via interpolation, because we don't have integer alphas and betas 
G_interpolation_fast<-function(alpha,beta)
{return(mean(Gmatrix[round(beta), round(alpha)], Gmatrix[round(beta), trunc(alpha)],  Gmatrix[trunc(beta), round(alpha)], Gmatrix[trunc(beta), trunc(alpha)]))}

############################################################################
# Choose true_state - multinomial with prob pi0 
#set.seed(123456789)
Draw_state<-function(STATES, prob) { return(which.is.max(rmultinom(STATES, 1, prob)))}
 
############################################################################
# Provides the initial morph as specified by INITIAL_MORPH
Get_starting_morph =  function (G_current_sm) {
        current_m <- rep(0,TOT_MORPHS)
        for (m in 1:TOT_MORPHS) {current_m[m] <-pi0  %*% G_current_sm[,m] }
        return (which.is.max(current_m) ) 
        }
  
############################################################################
## bayesian updater for q_r, with br_updated as prior
Bayesian_updater<-function(K, M, navigation_prior,new_omega, Ckjn,all_s, Y_mat )
{
    num=num_jointclicks=matrix(rep(0, K*TOT_STATES), ncol=TOT_STATES)
    den_jointclicks=rep(0, K)
    for (k in 1:K)
    {
        for (r in 1:TOT_STATES)
        {
            num[k,r]=prod((exp(Ckjn[Ckjn[,1]==k,2:ncol(Ckjn)]%*%new_omega%*%all_s[r,])/sum(exp(Ckjn[Ckjn[,1]==k,2:ncol(Ckjn)]%*%new_omega%*%all_s[r,])))^(Y_mat[Y_mat[,1]==k,2]))
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

## Bayesian updating - Montoya et al 2010
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



Bayesian_Thm_updater_element<-function(lambda_, P_st, morph, p_delta )
{
  numerator   <-  matrix(, TOT_STATES) 
  if (TOT_STATES ==2) 
  {
  # For clarity I am doing the expectation by hand, explicitly. Later I can generalize it for s =1..TOT_STATES 
  numerator[s=1] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=1, morph] * p_delta[morph, s_prime=1] + 
                     P_st[s_prime=2] * lambda_[s_prime=2, s=1, morph] * p_delta[morph, s_prime=2]      
       
  numerator[s=2] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=2, morph] * p_delta[morph, s_prime=1] + 
                     P_st[s_prime=2] * lambda_[s_prime=2, s=2, morph] * p_delta[morph, s_prime=2]      

  } else
  {
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

