#######################################################################################################################
# Author: 
# Purpose: Gets simulation results for webshop with one Morphing Opportunities
# 
# Note: 
#  - Make sure to have run 'Configuration.R' and 'Functions.R' before running this file
#  - After that, ensure the working directory is Replication_Morphing_HMM/Study2/3. Simulation Code
#
# Overview:
#     A) Set parameters and load data
#     B) Gather HMM estimates
#     C) Loop over all visitors
#
#
#
#######################################################################################################################


##########
# A) Set common parameters and load data
##########


# register parallal cores
fregisterDoParallel(cores=10)

# Context:  "PHONE STORE" or "HULB2009"
EMPIRICAL_SETTING     <- "PHONE STORE"  

# Model specs: "HULB","HMM_MAB_NO_DP","MAB_NO_HMM","HMM_BANDIT",  "HMM_BANDIT_STORE", "RANDOM" 
MODEL       <- "HMM_BANDIT_STORE" 

# HMM model to update states. This is also the folder used to store the HMM input files
HMM_MODEL   <- "geomtricHMM_estimates_April2021"

# Benchmarks
#    test1: the states are estimated by the HMM modek
#    test2: state 1 is the initial; state 2 is when a person add a product to the comparison or to the trolley
BENCHMARK   <-  "test1"   

# Terminal or non-terminal state
TERMINAL              <- TRUE 

# Number of states: 1 is HMM_MAB_NO_DP. 2 is  our chosen model
TOT_STATES            <- 2  

# When the morph decision is made
MORPHING_OPPORTUNITY  <- 7

# Do we have covariates in the transition matrix?
TRANSITION_COVARIATES <- TRUE
## the probability of clicking on "add-to-comparison" at click 7
COVARIATE_PROBABILITY <- .229 

# Arrival rates
QS_ARRIVAL        <- if(TOT_STATES==2) {c(1, 0)} else {c(1, 0, 0) } # Everyone starts in state 1 in the RCT
QS_ARRIVAL_NATURE <- if(TOT_STATES==2) {c(1, 0)} else {c(1, 0, 0)}  

# Global parameters, rarely change  
TOT_VISITORS <- 100000 
TOT_MORPHS   <- 2  
K_FULL       <- 14  #  4 in MBA/first two rounds: 4  periods so we allow to change morphs every 4 sets of clicks. Alina suggests (a) using 15 (clicks, RCT median=15, mean=20). (b) decouple clicks from periods. This should be thought thorough
TOT_PERIODS=K_FULL
TOT_CONSIDERED_PERIODS<-4


# path for hmm estimate data
PATH_HMM_EST = '../2. HMM Estimates'

  
# HMM Configuration files  
if (TRANSITION_COVARIATES) 
{
  FILENAME_MU  <- paste(PATH_HMM_EST,  "/", HMM_MODEL,"/mu_c14_2ST_RCT_April2021.csv" ,sep="")   
  FILENAME_RHO <- paste(PATH_HMM_EST,  "/", HMM_MODEL,"/rho_c14_2ST_RCT_April2021.csv",  sep="")  
}else
{
  FILENAME_MU  <- paste(PATH_HMM_EST,  "/", HMM_MODEL,"/mu_c14_nocov_2ST_RCT_April2021.csv" ,sep="") 
}


FILENAME_PSM_per_state  <- paste(PATH_HMM_EST, "/", HMM_MODEL,file= "/psm_condition_user_level_RCT_April2021_prePost.csv", sep="")
FILENAME_BOUNCE_PROBS  <- paste(PATH_HMM_EST, "/", HMM_MODEL,file= "/pmf_bounce_geometric_RCT_April2021.csv", sep="")  
FILENAME_GEOM_EM_PROBS  <- paste(PATH_HMM_EST, "/", HMM_MODEL,file= "/geometric_emission_probs_RCT_April2021.csv", sep="")


# Process the parameters already loaded
STATES      <- 1:TOT_STATES # needed for the draw_state funtion

HULB <- "HULB"; HMM_MAB_NO_DP <- "HMM_MAB_NO_DP"; MAB_NO_HMM <- "MAB_NO_HMM"; HMM_BANDIT_STORE <- "HMM_BANDIT_STORE"; HMM_BANDIT_4P <- "HMM_BANDIT_4P"; HMM_BANDIT <- "HMM_BANDIT" ; RANDOM <- "RANDOM" # do not change this
if (TRANSITION_COVARIATES) XA_dim  <- 1 else XA_dim=NULL
if (TOT_STATES==1) {all_s <- matrix(1,ncol=TOT_MORPHS)}; if (TOT_STATES==2) {all_s <- cbind(1,c(-1,1)) }; if (TOT_STATES==3) {all_s <- cbind(1, c(-1, 1, -1), c(-1, -1,1))}
INIT_STATES <- TOT_STATES

# Load G table
Gmatrix <- read.table(DROPBOX_GMATRIX_LINK)

# Load HMM estimated parameters - geometric emission probs (model-based probabilities to bounce)
geometric_emission_probs <- read.csv(FILENAME_GEOM_EM_PROBS, header=T)

# Load bounce probability based on raw number of clicks
temp_pbounce_vec  <- read.csv(FILENAME_BOUNCE_PROBS , header=T)


##########
# B) Gather HMM estimates
##########

# Load HMM estimated parameters ?
if (TOT_STATES > 1 ) {mu_vec <-read.csv(FILENAME_MU, header=T) }

# Load HMM estimated parameters - rho_vec is ...?
if (TRANSITION_COVARIATES) rho_vec <- read.csv( FILENAME_RHO , header=T)

# Load purchase probability 
temp_psm <- read.csv(FILENAME_PSM_per_state , header=T)  

# II.  Curate HMM Parameters: rho, mu, XA, Q0, Q0_K, psm   ####
if (TOT_STATES > 1 )
{ 
  # to start with I go with first morph matters: condition 1+ condition 2 and condition 3 + condition 4
  pbounce_vec <- matrix(, nc=2, nr=3);  rownames(pbounce_vec) <- c("State 1", "State 2", "State 3")
  pbounce_vec[state=1,] <- c(temp_pbounce_vec[1,2], temp_pbounce_vec[3,2]) # c(mean(st1-c1 & st1-c2),mean(st1-c3 & st1-c4))
  pbounce_vec[state=2,] <- c(temp_pbounce_vec[2,2], temp_pbounce_vec[4,2]) # c(mean(st2-c1 & st2-c2),mean(st2-c3 & st2-c4))
  pbounce_vec[state=3,] <- c(temp_pbounce_vec[5,2], temp_pbounce_vec[6,2]) # c(mean(st3-c1 & st3-c2),mean(st3-c3 & st3-c4))
  pbounce_vec <- cbind( c("S1M1", "S1M2", "S2M1", "S2M2","S3M1", "S3M2"), t(cbind(t(pbounce_vec[1,]),t(pbounce_vec[2,]),t(pbounce_vec[3,]))) )  
  pbounce_vec[,2]=as.numeric(pbounce_vec[,2])
  if (TOT_STATES==2) {pbounce_vec=pbounce_vec[-c(5, 6),]}
  
  # Not used - remove WTP page as exit will still be handled separately given the positioning of the paper
  # Ckjn <- Ckjn[ -NROW(Ckjn),]
  
  # Parameters 
  if (TERMINAL) {DIM_NON_TERMINAL_STATES=TOT_STATES-1} else {DIM_NON_TERMINAL_STATES=TOT_STATES}
  
  mu  <- array(mu_vec[1:2, 3], dim=c(DIM_NON_TERMINAL_STATES, (TOT_STATES-1), TOT_MORPHS)) 
  
  # Important - special case: if there are transition covariates, replace mu with the right mu
  if (TRANSITION_COVARIATES) {
    XA_dim=XA_dim
    rho <- array(rho_vec[  ,3],    dim=c(DIM_NON_TERMINAL_STATES,  XA_dim, TOT_MORPHS) )  
  }
  
  # Not used - XA: Covariates of the transition matrix - FIXED to XA=0 as start 
  if (TRANSITION_COVARIATES) {XA <-  c(0) }
  
  # Transition probabilities at the last state K: 0/1, because we re using terminal states.
  Q0_K      <- if(TOT_STATES==3) {array(rep(c(0,0,1), TOT_STATES), dim=c(1, TOT_STATES, TOT_MORPHS))} else
  {array(rep(c(0,1), TOT_STATES), dim=c(1, TOT_STATES, TOT_MORPHS))}
  
  # Not used - Q0 binds the transition matrix for previous state 1 to K-1, with the last row of the transition matrix, which is 0/1.
  # Q0=Transition_Matrix(XA, rho, mu, Q0_K, DIM_NON_TERMINAL_STATES, TERMINAL)
}  

# III. Compute key variables: Omega, Psm_exposure, p_bounce   ####
if (TOT_STATES ==1) 
{ # when TOT_STATES == 1, thus MODEL=MAB_NO_HMM 
  # Omega       <- matrix(c(0, 0, as.numeric(c(OOMEGA[5,2:3]))), ncol=2)
  psm_exposure <- matrix(c(temp_psm_exposure[2, 2], temp_psm_exposure[5, 2]), nrow=TOT_MORPHS)
  p_bounce_click <-matrix(as.numeric(c(temp_pbounce_vec[2,2], temp_pbounce_vec[5,2])), nc=TOT_MORPHS, nr=TOT_STATES )  
} 

if (TOT_STATES > 1)  
{  
  # true success probabilities at Morph / Timing (Pre/Post) / State level
  P_smt_true <- matrix(temp_psm[,3], ncol=TOT_STATES, byrow = T)
  ## true success probabilities at Morph / Timing (Pre/Post) level for state 2
  P_mt_post  <- P_smt_true[1:4 ,2]
  P_mt_pre  <- P_smt_true[5:6 ,2]
  
  # Compute the probability of bouncing. FORMAT: p_bounce_mt[period,init_stage, morph] 
  p_bounce_click <- matrix(as.numeric(pbounce_vec[,2]), nrow=TOT_STATES, byrow = T)
}

# Loads the psm - true purchase per click
if (BENCHMARK == "test2") {   }  

##########
# C) Loop over all visitors
##########

set.seed(9000)
trials=1000


# Loop over replicates -----
ptime <- system.time({
  sim_1k_reps_one_morph <- foreach(icount(trials), .packages=c('expm','nnet', 'maotai', 'tidyverse') ,.combine=rbind) %dopar% {
    # IV.   Initialize data structures used in next section for storage and loop control ####
    N        <- matrix(0, ncol=TOT_VISITORS)
    last_click  <- rep(0, TOT_VISITORS) #matrix(K_FULL, ncol=TOT_VISITORS)
    I        <- matrix(, ncol=TOT_MORPHS)
    delta_sm <- matrix(0, nrow = TOT_STATES, ncol=TOT_MORPHS); colnames (delta_sm ) <-c("abstract", "concrete"); rownames (delta_sm )<- if(TOT_STATES==2) {c("early", "late") } ; rownames (delta_sm ) <- if(TOT_STATES==3)    {c("early", "mid","late") } 
    I_evolution       <- matrix(0,nrow=TOT_VISITORS, ncol=TOT_MORPHS)
    morph_chosen      <- matrix(0,nrow=TOT_VISITORS, ncol=K_FULL+1)  
    posterior         <- matrix( , nrow=TOT_VISITORS, ncol=K_FULL)  
    true_state        <- rep(0, K_FULL)
    alphabeta         <- array(rep(1, TOT_STATES*TOT_MORPHS*2), dim=c(TOT_STATES, TOT_MORPHS, 2))
    
    ## We are not using PRE at all, we are just using POST and consider that only the morph seen from click 7 to 14 matters. 
    alpha_POST<- matrix(alphabeta[,,1], ncol=TOT_MORPHS) ; beta_POST <- matrix(alphabeta[,,2], ncol=TOT_MORPHS)
    G_current_sm_POST   <- matrix(0,  nrow=TOT_STATES, ncol=TOT_MORPHS) 
    
    # V.  Loop over all visitors: initializations, then Nature draws state, we update state estimates, get best morph from DP w/G matrix ####
    for (visitor  in 1:TOT_VISITORS)
    {
      # for testing
      # visitor=1
      ######### 1 Initialize data structures to store state estimates (qs_system is HULB, qs_HMM is HMM) ####
      if (TOT_STATES>1) {qs_HMM <- qs_HULB <- matrix(, nrow=K_FULL+1, ncol=TOT_STATES); qs_HMM[1,] <- QS_ARRIVAL; qs_HULB[1,] <- QS_ARRIVAL}
      qs_HMM_nature  <- matrix(, nrow=K_FULL+1, ncol=TOT_STATES)  
      qs_HMM_nature[1,] <- QS_ARRIVAL_NATURE
      
      ####  Draw clickstream - user-level click on the covariate ######
      XA <- sample(c(0, 1), 1, prob=c((1-COVARIATE_PROBABILITY), COVARIATE_PROBABILITY), replace=T)
      
      ######### 2 Load current G (use asymptotical value if alpha or beta are greater than 3000)  ####
      
      for (st_ in 1: TOT_STATES) {for (mor in 1:TOT_MORPHS)  
      {if (any(alpha_POST[st_,mor] >= 3000, beta_POST[st_,mor] >= 3000))
      {  G_current_sm_POST[st_, mor] <- alpha_POST[st_,mor]/(alpha_POST[st_,mor]+beta_POST[st_, mor]) } 
        else {G_current_sm_POST[st_, mor]<- G_interpolation_fast(alpha_POST[st_,mor],beta_POST[st_,mor])  } } }
      
      ######### 3 Assign the initial morph - random in all scenarios ####
      best_morph <- which(rmultinom(1,1,matrix(c(1/TOT_MORPHS),ncol=TOT_MORPHS) ) == 1)  ## do we need to update here with G-PRE?
      morph_chosen[visitor,1]  <- best_morph 
      
      ######### 4 Loop over all clicks: Nature draws state, update our state estimates, find best morph from DP w/ G matrix    ####
      for (click in 1:K_FULL)     
      { 
        # for testing
        # click=1 
        # Draw true state:  S(t)= f(qs_HMM(t) )  but this does not take into account effects of morph exposure ####
        if (BENCHMARK == "test1") {true_state[click] <- Draw_state(STATES, qs_HMM_nature[click, ]) }
        if (BENCHMARK == "test2") {true_state[click] <- Draw_state(STATES, qs_HMM_nature[click, ]) }
        
        # Draw bounce: last_click[visitor] will have the number of the last click ####
        if (click < K_FULL) 
        { if (TOT_STATES==1) {temp_prob <- p_bounce_click[morph_chosen[visitor, click]] } 
          if (TOT_STATES >1) {temp_prob <- p_bounce_click[true_state[click], morph_chosen[visitor, click]] }
          bounced_at_click <- sample(c(1,0), 1, c(temp_prob, 1-temp_prob), replace=TRUE)
        } 
        
        # Stop clickstream if bounced ####
        if (bounced_at_click == 1) {break}  
        
        # Update state estimates and find best morph  #### 
        if (TOT_STATES==1) {for (m in 1:TOT_MORPHS) {if(any(beta[ ,m]> 2999) | any(alpha[ ,m]> 2999)) {I[m]<-alpha[1, m]/(alpha[1, m]+beta[1, m])} else {I[m]<-G_interpolation_fast(alpha[1, m], beta[1, m]) } }; I_evolution[visitor,] <- I;  best_morph <- which.is.max(I[])}
        if (TOT_STATES> 1)
        { 
          # A. Computes transition matrix including all probs. Stack it ####
          lambda   <- Transition_Matrix(XA, rho, mu, Q0_K, DIM_NON_TERMINAL_STATES, TERMINAL) 
          lambda_m <- lambda[,,best_morph]  # Select the matrix corresponding to the morph served 
          
          # B. Update NATURE state for next click estimate $QS_NATURE$ = f($Q_{mt}$)  ####
          qs_HMM_nature[click+1,] <- lambda_m[true_state[click], ]  
          
          # C. Update algo state estimates: HMM, HULB, HMM of future states ###
          if (is.element(MODEL, c(HMM_BANDIT, HMM_MAB_NO_DP, HMM_BANDIT_4P, HMM_BANDIT_STORE, RANDOM) ) & BENCHMARK  == "test1"  ) 
          {qs_HMM[click+1,] <- Bayesian_Thm_updater_element_cond(lambda, qs_HMM[click,], best_morph, geometric_emission_probs) }   
          
          if (is.element(MODEL, c(HMM_BANDIT, HMM_MAB_NO_DP, HMM_BANDIT_4P, HMM_BANDIT_STORE) ) & BENCHMARK  == "test2"  ) 
          {qs_HMM[click+1,] <- 1  } # add here test 2 condition 
          
          # D. Update HULB state estimate: qs_HULB= f(click (t+1), m(t)  )  ####
          #    qs_HULB[click+1,] <- softmax(t(c(1, log(click+1))) %*% Omega) 
          
          # E  Update HMM estimate of future states: qs_HMM_sys_future = f(click (t+1), m(t) , lambda ) ####
          if (is.element(MODEL, c(HMM_BANDIT, HMM_MAB_NO_DP, HMM_BANDIT_4P, HMM_BANDIT_STORE)) & BENCHMARK  == "test1") { qs_HMM_sys_future <- estimate_future_period_qs_general(lambda)} #, qs_HMM[click+1,] ) } # Gui replaced lambda instead of q0
          if (MODEL == HMM_BANDIT & BENCHMARK  == "test2") { qs_HMM_sys_future <- estimate_future_period_qs_test2(lambda)} #, qs_HMM[click+1,] ) } # Gui replaced lambda instead of q0
          
          # F. Time to Morph? Then (1) compute EGI for second morph (_POST), (2) Solve DP, (3) Assign best_morph  ####
          if (click == MORPHING_OPPORTUNITY) 
          {
            # F.1. Compute EGI: m*(t=1) = Argmax EGI(qs_hulb(t+1)  ) ####
            if (MODEL == HMM_MAB_NO_DP)    {qs <- qs_HMM[click,  ] } 
            if (MODEL == HULB)             {qs <- qs_HULB[click, ] } 
            if ( is.element(MODEL, c(HULB, HMM_MAB_NO_DP) )    )     
            {
              for (m in 1:TOT_MORPHS) 
              {
                if (any(beta_POST[ ,m]> 2999) | any(alpha_POST[ ,m]> 2999)) 
                {
                  I[m] <- sum( (alpha_POST[ ,m]/ (alpha_POST[ ,m]+beta_POST[ ,m])) * qs )
                } else {
                  EGI  <- mapply(G_interpolation_fast, alpha_POST[,m], beta_POST[,m]) 
                  I[m] <- sum(EGI *qs)     
                } # close if
              }# close for
              I_evolution[visitor,] <- I 
            }        
            
            # F.2. Solve DP using HMM's m*(t=1) = Argmax DP(qs_HMM(t+1), qs_future(t+1) )  ####
            if (MODEL == HMM_BANDIT) {best_morph_uncond <- DP(click, qs_HMM_sys_future, qs_HMM[click,], G_current_sm_POST, p_bounce_click )  } # , delta_dp= delta_dp[period+1]) }
            if  (MODEL == HMM_BANDIT_4P) {best_morph_uncond <- DP_4periods(click, qs_HMM_sys_future, qs_HMM[click,], G_current_sm_POST, p_bounce_click ) } # , delta_dp= delta_dp[period+1]) }
            if  (MODEL == HMM_BANDIT_STORE) {best_morph_uncond <- DP_webstore_implementation(G_current_sm_POST) }  # , delta_dp= delta_dp[period+1]) }
            
            # F.3. Assign best morph for the next click according to the method of choice (breaking ties randomly)####
            if (MODEL == HMM_BANDIT)       { best_morph <- best_morph_uncond[click] }  # note that index should be click, not TOT_PERIODS
            if  (MODEL == HMM_BANDIT_4P)    {best_morph<- best_morph_uncond[1] } # , delta_dp= delta_dp[period+1]) }
            if  (MODEL == HMM_BANDIT_STORE) {best_morph <- best_morph_uncond }  # , delta_dp= delta_dp[period+1]) }
            if (MODEL == HMM_MAB_NO_DP)    { best_morph <- which.is.max(I[]) }     
            if (MODEL == HULB)             { best_morph <- which.is.max(I[]) }   
            if (MODEL == RANDOM)           { best_morph <- which(rmultinom(1,1,matrix(c(1/TOT_MORPHS),ncol=TOT_MORPHS) ) == 1) }   
          } # if time to morph 
          
        } #  close if TOT_STATES >1
        
        # Store next click's new best morph
        morph_chosen[visitor, click+1] <-best_morph
        
      }  # close click loop  
      
      # mark last click: bounced or purchase
      last_click[visitor] <-LAST_CLICK<- click  
      
      ######### 5 Compute Psm_i ####
      
      # 5.1 If only saw the random morph 0, take the probbaility from m11, and m22
      if (click < MORPHING_OPPORTUNITY ) 
      { 
        if (best_morph==1) psm_i <- P_mt_pre[1] else psm_i <- P_mt_pre[2]
      }
      
      # 5.2 If saw the random morph and one more then use the corresponding probability using the four conditions measured in the field
      if (click >=  MORPHING_OPPORTUNITY )  
      {  
        
        ## take P_mt_true: the true success probabilities at state 2 only
        if (morph_chosen[visitor, MORPHING_OPPORTUNITY]== 1 & morph_chosen[visitor, LAST_CLICK]== 1) { psm_i <- P_mt_post[1]  } # m1m1
        if (morph_chosen[visitor, MORPHING_OPPORTUNITY]== 1 & morph_chosen[visitor, LAST_CLICK]== 2) { psm_i <- P_mt_post[2]  } # m1m2
        if (morph_chosen[visitor, MORPHING_OPPORTUNITY]== 2 & morph_chosen[visitor, LAST_CLICK]== 1) { psm_i <- P_mt_post[3]  } # m2m1
        if (morph_chosen[visitor, MORPHING_OPPORTUNITY]== 2 & morph_chosen[visitor, LAST_CLICK]== 2) { psm_i <- P_mt_post[4]  } # m2m2
      }
      
      
      ######### 6 Draw purchase using Psm_i given probs from conditions and the person's specific condition ####
      if (TOT_STATES == 1) {delta <-  rbinom(1, 1, psm_i); delta_sm[state=1, morph_chosen[visitor,LAST_CLICK]]<- delta_sm[1, morph_chosen[visitor,LAST_CLICK] ] + delta; fraction <- rep(1,TOT_STATES) } 
      if (TOT_STATES > 1 )
      {  
        # 6.3 Draws purchase  #### 
        
        # Gets prm from morph seen
        delta <-  rbinom(1, 1, psm_i)
        delta_sm[true_state[LAST_CLICK], morph_chosen[visitor,LAST_CLICK] ]<- delta_sm[true_state[LAST_CLICK], morph_chosen[visitor,LAST_CLICK] ] + delta # changed to LAST_CLICK
        
        if (is.element(MODEL, c(RANDOM, HMM_BANDIT, HMM_MAB_NO_DP, HMM_BANDIT, HMM_BANDIT_4P, HMM_BANDIT_STORE) ) )  {qs_algo <-qs_HMM[LAST_CLICK,] } else {qs_algo <- qs_HULB[LAST_CLICK,] }   # set which state estimate to use
        ## use fractional updating with state probability estimates to update alpha and beta
        fraction <- qs_algo       
      }  
      
      ######### 7 Learn: updates alpha, beta  with qs_algo depending  on method being HULB or HMM ####
      for (s in 1:TOT_STATES) 
      {
        alpha_POST[s, best_morph] <- alpha_POST[s, best_morph]+      delta * fraction[s]
        beta_POST[s, best_morph]  <- beta_POST[s, best_morph] + (1 - delta)* fraction[s]   
      }        
      
    } # close visitor loop  
    
    Success_rate = sum(delta_sm)/TOT_VISITORS
    sumdelta=sum(delta_sm)
    c(Success_rate, sumdelta)
  }
  
})[3]


