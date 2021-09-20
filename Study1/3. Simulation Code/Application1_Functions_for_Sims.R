# Replication code for "Morphing for Consumer Dynamics: Bandits Meet HMM"
# Author: Gui Liberali and Alina Ferecatu 
# Date: October 2021
# Analysis for Study 1 - MBA study
# Functions necessary for the empirically-grounded simulations in Application 1

############################################################################
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
  return(Aij)
}

############################################################################
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

## WITHIN-PERSON  DP
DP = function(K, qs_future, qs, G_current_sm, p_bounce)  
{
  V_immediate    <- V_continuation <- array(0,c(TOT_STATES, TOT_MORPHS, TOT_PERIODS) )  
  V_total_cond   <- array(0,c(TOT_STATES, TOT_MORPHS, TOT_PERIODS+1) )  
  V_total        <- array(0,c(TOT_MORPHS, TOT_PERIODS) ) 
  V_max          <- rep(0,TOT_PERIODS)  
  best_morph     <- rep(0,TOT_PERIODS)
  V_total_cond[,,K_FULL+1] <-  0  
  
  for (period in K_FULL:K)
  {			
    next_period <- period+1
    for (morph in 1: TOT_MORPHS) 
    {
      for (stage in 1:TOT_STATES) 
      { 
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
Draw_state<-function(STATES, prob) { return(which.is.max(rmultinom(STATES, 1, prob)))}
 
############################################################################
## Bayesian updater for MAB_FIXED_STATES
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

## Bayesian updating for the HMM
##  Value: prob of state s' at t+1 given state s at t, organized in a K*K matrix, where on the rows you have the current state, and on the columns you have the previous states.  
##  Intermediate steps:
##    numerator: probability of ending up in state s' given all possible starting points s
##    denominator: sum over all probabilities of being in a state s', starting from state s and summed acrross all previous states s
Bayesian_Thm_updater_element<-function(lambda_, P_st, morph, p_delta )
{
  numerator   <-  matrix(, TOT_STATES) 
  if (TOT_STATES ==2) 
  {
  numerator[s=1] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=1, morph] * p_delta[morph, s_prime=1] + 
                     P_st[s_prime=2] * lambda_[s_prime=2, s=1, morph] * p_delta[morph, s_prime=2]      
       
  numerator[s=2] <-  P_st[s_prime=1] * lambda_[s_prime=1, s=2, morph] * p_delta[morph, s_prime=1] + 
                     P_st[s_prime=2] * lambda_[s_prime=2, s=2, morph] * p_delta[morph, s_prime=2]      

  } else
  {
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
