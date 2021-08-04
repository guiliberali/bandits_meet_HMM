##################### Morphing Consumer Dynamics - DP code    #############                                       
#  Coded by Gui Liberali in December 2018/January 2019 using HULB & functions + HMM csv's from Alina Ferecatu
#  Updated by Gui Liberali and Alina Ferecatu several times. See file "Log of updates.txt" for details.

#  IMPORTANT: Run Config.R before running this code.

# I.   Loaded parameters (some were initialized in the Config files) ####

# Load functions and support code. IMPORTANT: Run Config.R before running this code.
source(paste (PATH, "/Functions.R" , sep="") )   

# Process the parameters already loaded
STATES      <- 1:TOT_STATES # needed for the draw_state funtion
 
HULB <- "HULB"; HMM_MAB_NO_DP <- "HMM_MAB_NO_DP"; MAB_NO_HMM <- "MAB_NO_HMM"; HMM_BANDIT_STORE <- "HMM_BANDIT_STORE"; HMM_BANDIT_4P <- "HMM_BANDIT_4P"; HMM_BANDIT <- "HMM_BANDIT" ; RANDOM <- "RANDOM" # do not change this
if (TRANSITION_COVARIATES) XA_dim  <- 1 else XA_dim=NULL
if (TOT_STATES==1) {all_s <- matrix(1,ncol=TOT_MORPHS)}; if (TOT_STATES==2) {all_s <- cbind(1,c(-1,1)) }; if (TOT_STATES==3) {all_s <- cbind(1, c(-1, 1, -1), c(-1, -1,1))}
INIT_STATES <- TOT_STATES

# Load G table
setwd(PATH_IN);setwd("../..") # Gmatrix.out should not be in the github folder due to size
#if (!GUI_IS_RUNNING) {setwd(PATH_G_MATRIX)}
Gmatrix <- read.table(paste0('Github_HMM/', FILENAME_G))

# Load Ckjn with website structure
# This is not needed here. I commented it out.
# Ckjn <- read.csv(FILENAME_CKJN, header=T)  

# Load HMM estimated parameters - Omega, to learn states from clicks
# OOMEGA  <- read.csv( FILENAME_OOMEGA,  header=T) 

# Load HMM estimated parameters - geometric emission probs (model-based probabilities to bounce)
geometric_emission_probs <- read.csv(FILENAME_GEOM_EM_PROBS, header=T)

# Load bounce probability based on raw number of clicks
temp_pbounce_vec  <- read.csv(FILENAME_BOUNCE_PROBS , header=T)

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
      # For now:1st morph matters only. Details on Github spreadsheet labeled "From Conditions to states.xslx"
      # This needs to be adjusted for exposure and ordering once it is clear how the updating will assign credits
      # Current procedure:
        # psm (m1 |state 1)=  mean(cond1|state1,cond2|state1)       psm (m2 |state 1)=  mean(cond3|state1,cond4|state1) 
        # psm (m1 |state 2)=  mean(cond1|state2,cond2|state2)       psm (m2 |state 2)=  mean(cond3|state2,cond4|state2) 
        # psm_true <- matrix(, nc=2,nr=2); colnames(psm_true) <- c("morph 1", "morph 2") ; rownames(psm_true) <- c("State 1", "State 2" )
        # psm_true[state=1,] <- c(mean(c(temp_psm[1,2],temp_psm[3,2]), na.rm=T), mean(temp_psm[5,2],temp_psm[7,2]))  
        # psm_true[state=2,] <- c(mean(c(temp_psm[2,2],temp_psm[4,2]), na.rm=T), mean(temp_psm[6,2],temp_psm[8,2])) 
        
        ## PROPOSAL 
        # true success probabilities at Morph / Timing (Pre/Post) / State level
        P_smt_true <- matrix(temp_psm[,3], ncol=TOT_STATES, byrow = T)
        ## true success probabilities at Morph / Timing (Pre/Post) level for state 2
        P_mt_post  <- P_smt_true[1:4 ,2]
        P_mt_pre  <- P_smt_true[5:6 ,2]
        
       # Compute the probability of bouncing. FORMAT: p_bounce_mt[period,init_stage, morph] 
       p_bounce_click <- matrix(as.numeric(pbounce_vec[,2]), nrow=TOT_STATES, byrow = T)
       
       # Computes click probabilities, inspect them 
       # zerosXA         <- rep(0, XA_dim*TOT_STATES)
       # Omega - intercept + slope showing the impact of the click number on the probability of being in state 1 vs. state 2
       # Baseline is state 1, so the coeficients are 0.
       # Omega           <- matrix(c(zerosXA, as.numeric(c(OOMEGA[5,2:3]))), ncol=TOT_STATES) # cbind(zerosXA, as.matrix(OOMEGA[,2:TOT_PAGES])) 
     }

# Loads the psm - true purchase per click
if (BENCHMARK == "test2") {   }  

# IV.   Initialize data structures used in next section for storage and loop control ####
N        <- matrix(0, ncol=TOT_VISITORS)
last_click  <- rep(0, TOT_VISITORS) #matrix(K_FULL, ncol=TOT_VISITORS)
I        <- matrix(, ncol=TOT_MORPHS)
delta_sm_pre<- delta_sm_post <- matrix(0, nrow = TOT_STATES, ncol=TOT_MORPHS); colnames (delta_sm_pre ) <- colnames(delta_sm_post) <- c("abstract", "concrete"); 
rownames (delta_sm_pre )<- rownames (delta_sm_post )<- if(TOT_STATES==2) {c("early", "late") } ; 
rownames (delta_sm_pre )<- rownames (delta_sm_post ) <- if(TOT_STATES==3)    {c("early", "mid","late") } 
I_evolution       <- matrix(0,nrow=TOT_VISITORS, ncol=TOT_MORPHS)
morph_chosen      <- matrix(0,nrow=TOT_VISITORS, ncol=K_FULL+1)  
posterior         <- matrix( , nrow=TOT_VISITORS, ncol=K_FULL)  
true_state        <- rep(0, K_FULL)
true_states_everyone     <- qs_HULB_everyone_error <- qs_HMM_everyone_error <- 
qs_HMM_nature_everyone_error   <- matrix(0, nrow=TOT_VISITORS, ncol=K_FULL) 
error_state       <- matrix(0, nrow=TOT_VISITORS, ncol= 9); colnames(error_state)<- c("visitor", "qr_t1","qr_t2","qr_t3","qr_t4","qs_t1","qs_t2","qs_t3","qs_t4")
alphabeta         <- array(rep(1, TOT_STATES*TOT_MORPHS*2), dim=c(TOT_STATES, TOT_MORPHS, 2))
success_per_visitor <- rep(0, TOT_VISITORS)
success_prob_runAverage_pre <-success_prob_runAverage_post <- rep(0, TOT_VISITORS)

## We are not using PRE at all, we are just using POST and consider that only the morph seen from click 7 to 14 matters. 
alpha_PRE         <- alpha_POST<- matrix(alphabeta[,,1], ncol=TOT_MORPHS) ; beta_PRE<- beta_POST <- matrix(alphabeta[,,2], ncol=TOT_MORPHS)
G_current_sm_PRE  <- G_current_sm_POST   <- matrix(0,  nrow=TOT_STATES, ncol=TOT_MORPHS) 
G_asym_m1_PRE     <- G_asym_m1_POST     <- G_asym_m2_PRE <- G_asym_m2_POST    <-  matrix(0, nrow=TOT_VISITORS, ncol=TOT_STATES)  # current value of the  index
monitor_G_m1_s_PRE<- monitor_G_m2_s_PRE <- monitor_G_m1_s_POST<- monitor_G_m2_s_POST <- matrix(0, nrow=TOT_VISITORS, ncol=TOT_STATES)
MODEL
set.seed(9000)
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
        {if (any(alpha_PRE[st_,mor] >= 3000,  beta_PRE[st_,mor] >= 3000)) 
      {  G_current_sm_PRE[st_, mor]  <- alpha_PRE[st_,mor]  /(alpha_PRE[st_,mor] +beta_PRE[st_,mor])  } 
        else {G_current_sm_PRE[st_, mor] <- G_interpolation_fast(alpha_PRE[st_,mor], beta_PRE[st_,mor])  } } }
      
      for (st_ in 1: TOT_STATES) {for (mor in 1:TOT_MORPHS)  
      {if (any(alpha_POST[st_,mor] >= 3000, beta_POST[st_,mor] >= 3000))
      {  G_current_sm_POST[st_, mor] <- alpha_POST[st_,mor]/(alpha_POST[st_,mor]+beta_POST[st_, mor]) } 
        else {G_current_sm_POST[st_, mor] <- G_interpolation_fast(alpha_POST[st_,mor],beta_POST[st_,mor])  } } }

      ######### 3 Assign the initial morph - random in all scenarios ####
      best_morph <- which(rmultinom(1,1,matrix(c(1/TOT_MORPHS),ncol=TOT_MORPHS) ) == 1)  ## do we need to update here with G-PRE?
      
      if (MODEL == HMM_BANDIT) {
        lambda   <- Transition_Matrix(XA, rho, mu, Q0_K, DIM_NON_TERMINAL_STATES, TERMINAL) 
        qs_HMM_sys_future <- estimate_future_period_qs_general(lambda)
        best_morph_uncond <- DP(1, qs_HMM_sys_future, qs_HMM[1,], G_current_sm_POST, G_current_sm_PRE, p_bounce_click )  # , delta_dp= delta_dp[period+1]) }

        best_morph <- best_morph_uncond[1]
        } # , delta_dp= delta_dp[period+1]) }
      
      
      morph_chosen[visitor, 1]  <- best_morph 
          
      ######### 4 Loop over all clicks: Nature draws state, update our state estimates, find best morph from DP w/ G matrix    ####
      for (click in 1:K_FULL)     
        { 
        # for testing
        # click=8
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
                       if (MODEL == HMM_BANDIT) {best_morph_uncond <- DP(click, qs_HMM_sys_future, qs_HMM[click,], G_current_sm_POST, G_current_sm_PRE, p_bounce_click ) } # , delta_dp= delta_dp[period+1]) }
                       if  (MODEL == HMM_BANDIT_4P) {best_morph_uncond<- DP_4periods(click, qs_HMM_sys_future, qs_HMM[click,], G_current_sm_POST, p_bounce_click ) } # , delta_dp= delta_dp[period+1]) }
                       if  (MODEL == HMM_BANDIT_STORE) {best_morph_uncond <- DP_webstore_implementation(G_current_sm_POST) }  # , delta_dp= delta_dp[period+1]) }
                   
                      # F.3. Assign best morph for the next click according to the method of choice (breaking ties randomly)####
                       if  (MODEL == HMM_BANDIT_4P)    {best_morph<- best_morph_uncond[1] } # , delta_dp= delta_dp[period+1]) }
                       if  (MODEL == HMM_BANDIT_STORE) {best_morph <- best_morph_uncond }  # , delta_dp= delta_dp[period+1]) }
                       if (MODEL == HMM_BANDIT)       { best_morph <- best_morph_uncond[click] }  # note that index should be click, not TOT_PERIODS
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

      # 5.2 If saw a random morph twice then use the corresponding probability using the four conditions measured in the field
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
             # 6.1 Stores true state ####
             true_states_everyone[visitor, ]     <- true_state
             
             # 6.2 Stores all states at 7th click (when we make a choice)
             # qs_HULB_everyone_error[visitor,  7]        <- (qs_HULB[LAST_CLICK,      true_states_everyone[visitor, 7]  ] - true_states_everyone[visitor, 7])^2
             ## In this equation, I replaced 7 with LAST_CLICK
             # qs_HMM_everyone_error[visitor,  7]        <- (qs_HMM[LAST_CLICK,      true_states_everyone[visitor, 7]  ] - true_states_everyone[visitor, 7])^2
             ## Second term should be 1, as we are taking the probability at the true state
             qs_HMM_everyone_error[visitor,  LAST_CLICK]        <- (qs_HMM[LAST_CLICK,       true_states_everyone[visitor, LAST_CLICK]  ] - 1)^2 
             qs_HMM_nature_everyone_error[visitor, LAST_CLICK]  <- (qs_HMM_nature[LAST_CLICK,true_states_everyone[visitor, LAST_CLICK]  ] - 1)^2 
                
             # 6.3 Draws purchase  #### 
             
             # Gets prm from morph seen
             delta <-  rbinom(1, 1, psm_i)
             
             if (is.element(MODEL, c(RANDOM, HMM_BANDIT, HMM_MAB_NO_DP, HMM_BANDIT, HMM_BANDIT_4P, HMM_BANDIT_STORE) ) )  
             {qs_algo <-qs_HMM[LAST_CLICK,] 
              qs_algo_before_morphing<-qs_HMM[MORPHING_OPPORTUNITY-1, ]
               } else {
                 qs_algo <- qs_HULB[LAST_CLICK,] 
                 qs_algo_before_morphing <- qs_HULB[MORPHING_OPPORTUNITY-1, ]
              }   # set which state estimate to use
             ## use fractional updating with state probability estimates to update alpha and beta
             fraction <- qs_algo      
             fraction_before_morphing <- qs_algo_before_morphing
          }  
       
      ######### 7 Learn: updates alpha, beta  with qs_algo depending  on method being HULB or HMM ####
        if (click < MORPHING_OPPORTUNITY ) {
          delta_sm_pre[true_state[LAST_CLICK], morph_chosen[visitor, LAST_CLICK] ]<- delta_sm_pre[true_state[LAST_CLICK], morph_chosen[visitor,LAST_CLICK] ] + delta # changed to LAST_CLICK
           eta_1=1
            for (s in 1:TOT_STATES) 
                {
                 alpha_PRE[s, best_morph] <- alpha_PRE[s, best_morph]+      delta * fraction[s]*eta_1
                 beta_PRE[s, best_morph]  <- beta_PRE[s, best_morph] + (1 - delta)* fraction[s]*eta_1   
            }    
                  
          }else{ 
            eta_1=(MORPHING_OPPORTUNITY-1)/LAST_CLICK
            eta_2=(LAST_CLICK-MORPHING_OPPORTUNITY+1)/LAST_CLICK
            
            delta_sm_pre[true_state[(MORPHING_OPPORTUNITY-1)], morph_chosen[visitor,(MORPHING_OPPORTUNITY-1)] ]<- delta_sm_pre[true_state[(MORPHING_OPPORTUNITY-1)], morph_chosen[visitor, (MORPHING_OPPORTUNITY-1)] ] + delta *eta_1
            delta_sm_post[true_state[LAST_CLICK], morph_chosen[visitor,LAST_CLICK] ]<- delta_sm_post[true_state[LAST_CLICK], morph_chosen[visitor,LAST_CLICK] ] + delta*eta_2 # changed to LAST_CLICK
            
              for (s in 1:TOT_STATES) 
                {
              alpha_PRE[s, morph_chosen[visitor, (MORPHING_OPPORTUNITY-1)]] <- alpha_PRE[s, morph_chosen[visitor, (MORPHING_OPPORTUNITY-1)]]+      delta * fraction_before_morphing[s]*eta_1
              beta_PRE[s, morph_chosen[visitor, (MORPHING_OPPORTUNITY-1)]]  <- beta_PRE[s, morph_chosen[visitor, (MORPHING_OPPORTUNITY-1)]] + (1 - delta)* fraction_before_morphing[s]*eta_1 
               
              alpha_POST[s, best_morph] <- alpha_POST[s, best_morph]+      delta * fraction[s]*eta_2
              beta_POST[s, best_morph]  <- beta_POST[s, best_morph] + (1 - delta)* fraction[s]*eta_2  
            }
         
        }
          
      ######## 8 Monitor G for pre and post ####
      for (s in 1:TOT_STATES) 
        {
         monitor_G_m1_s_PRE[visitor,s]  <- G_current_sm_PRE[s,1]; monitor_G_m2_s_PRE[visitor,s]  <- G_current_sm_PRE[s,2]
         G_asym_m1_PRE[visitor,s] <-  alpha_PRE[s,1]/(alpha_PRE[s,1]+beta_PRE[s,1]); G_asym_m2_PRE[visitor,s]       <-  alpha_PRE[s,2]/(alpha_PRE[s,2]+beta_PRE[s,2])
         
         monitor_G_m1_s_POST[visitor,s]  <- G_current_sm_POST[s,1]; monitor_G_m2_s_POST[visitor,s]  <- G_current_sm_POST[s,2]
         G_asym_m1_POST[visitor,s] <-  alpha_POST[s,1]/(alpha_POST[s,1]+beta_POST[s,1]);         G_asym_m2_POST[visitor,s] <-  alpha_POST[s,2]/(alpha_POST[s,2]+beta_POST[s,2])
      } 
      
      success_per_visitor[visitor] <- delta
      success_prob_runAverage_pre[visitor]=sum(delta_sm_pre)/visitor
      success_prob_runAverage_post[visitor]=sum(delta_sm_post)/visitor
      
              
} # close visitor loop  
      
## monitor G for m1 vs m2 POST -----
colnames(monitor_G_m1_s_POST)<- colnames(monitor_G_m2_s_POST)<-c("State 1", "State 2")
monitor_G=as_tibble(monitor_G_m1_s_POST) %>% 
  mutate(Visitor=1:TOT_VISITORS) %>% 
  pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
  mutate(Morph=c("Morph 1")) %>% 
  bind_rows(as_tibble(monitor_G_m2_s_POST) %>% 
              mutate(Visitor=1:TOT_VISITORS) %>% 
              pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
              mutate(Morph=c("Morph 2")) ) %>% 
  mutate(state_morph=case_when(State== 'State 1' & Morph== "Morph 1" ~ "S1M1",
                               State== 'State 2' & Morph== "Morph 1" ~ "S2M1",
                               State== 'State 1' & Morph== "Morph 2" ~ "S1M2",
                               TRUE ~ "S2M2"))

monitor_G %>% 
  filter(Visitor %in% 2000:100000) %>% 
  ggplot(aes(x=Visitor, y=Index, color=state_morph))+
  labs(title="Gittins index POST")+
  geom_line()


## monitor G for m1 vs m2 PRE -----
colnames(monitor_G_m1_s_PRE)<- colnames(monitor_G_m2_s_PRE)<-c("State 1", "State 2")
monitor_G=as_tibble(monitor_G_m1_s_PRE) %>% 
  mutate(Visitor=1:TOT_VISITORS) %>% 
  pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
  mutate(Morph=c("Morph 1")) %>% 
  bind_rows(as_tibble(monitor_G_m2_s_PRE) %>% 
              mutate(Visitor=1:TOT_VISITORS) %>% 
              pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
              mutate(Morph=c("Morph 2")) ) %>% 
  mutate(state_morph=case_when(State== 'State 1' & Morph== "Morph 1" ~ "S1M1",
                               State== 'State 2' & Morph== "Morph 1" ~ "S2M1",
                               State== 'State 1' & Morph== "Morph 2" ~ "S1M2",
                               TRUE ~ "S2M2"))

monitor_G %>% 
  filter(Visitor %in% 2000:50000) %>% 
  ggplot(aes(x=Visitor, y=Index, color=state_morph))+
  labs(title="Gittins index PRE")+
  geom_line()

## true success rate vs. model estimated one
mean(P_mt_post); mean(P_mt_pre); 
(success_rate_pre  <- sum(delta_sm_pre)/TOT_VISITORS)
(success_rate_post  <- sum(delta_sm_post)/TOT_VISITORS)
table(morph_chosen[,K_FULL])  

(overall_success_rate = (sum(delta_sm_pre) + sum(delta_sm_post)) / TOT_VISITORS)

## monitor G for m1 vs m2, in state 2
plot(monitor_G_m1_s_POST[10000:40000, 1], type="l")  
lines(monitor_G_m2_s_POST[10000:40000, 1], col="red")  

plot(G_asym_m1_POST[5000:40000, 1], type="l")  
lines(G_asym_m2_POST[5000:40000, 1], col="red")  

plot(G_asym_m1_POST[10000:40000, 1], type="l")  
lines(G_asym_m2_POST[10000:40000, 2], col="red")  

plot( success_prob_runAverage[3000:40000], type="l")


# VI. Results - successes per state and morph #############
if (!GUI_IS_RUNNING) { setwd(PATH_OUT) } 
save.image("HMM_BANDIT_single_run.RData")

load("Random_single_run.RData")

### PLOTS RANDOM vs. various treatments
colnames(monitor_G_m1_s_POST)<- colnames(monitor_G_m2_s_POST)<-c("State 1", "State 2")
monitor_Success_Prob=tibble(success_prob=success_prob_runAverage, 
                            Treatment="HMM_Bandit_Webstore")

monitor_Success_Prob=monitor_Success_Prob %>% 
  bind_rows(tibble(success_prob=success_prob_runAverage, 
         Treatment="Random") )%>% 
  mutate(Visitor=rep(1:TOT_VISITORS, 2) )

monitor_Success_Prob %>% 
  filter(Visitor %in% 5000:20000) %>% 
  ggplot(aes(x=Visitor, y=success_prob, color=Treatment))+
  geom_line()

monitor_Success_Prob %>% filter(Visitor %in% 2000:4000)
  # Check last-state membership
  last_state <-last_morph <- last_qs_HULB <- last_qs_HMM <- last_qs_nature <-  true_state_binary <-  matrix (0, ncol = TOT_VISITORS)
  if (TOT_STATES==1)
    { 
     last_state <- rep(1,ncol = TOT_VISITORS)
    } else
    {
      for (iii in 1: TOT_VISITORS ) 
          {  last_state[iii]     <- true_states_everyone[iii,bounced[iii]]
             last_morph[iii]     <- morph_chosen[iii,bounced[iii]]
            # last_qs_HULB[iii]   <- qs_HULB_everyone[iii,bounced[iii]]
             last_qs_HMM[iii]    <- qs_HMM_everyone[iii,bounced[iii]]
             last_qs_nature[iii] <- qs_HMM_nature_everyone[iii,bounced[iii]]
             if(last_state[iii]==1) { true_state_binary[iii] <- 1} else { true_state_binary[iii] <- 0 } 
           }
    }
 
   
   # Get RMSE
  # RMSE_HULB_7 <-sqrt(mean(qs_HULB_everyone_error[ ,  7],na.rm=T ))
   RMSE_HMM_7  <-sqrt(mean(qs_HMM_everyone_error[ ,  7],na.rm=T ))
     
   
   # Store replicate results
    last_state_and_morph           <- t(rbind(last_state, last_morph, last_qs_HULB, last_qs_HMM, last_qs_nature, true_state_binary ))
    colnames(last_state_and_morph) <- c("state", "morph", "qs_HULB_sl", "qs_HMM_sl", "qs_hMM_nature_sl", "last_true_state_binary")
    test                            <- last_state_and_morph
    HULB_abs_error_replicates       <- sum(abs(test[,3]-test[,6]))
    HMM_abs_error_replicates        <- sum(abs(test[,4]-test[,6]))
    HMM_NATURE_abs_error_replicates <- sum(abs(test[,5]-test[,6])) 
    
    delta_replicates_post                <- c(delta_sm_post)
    sum_delta_replicates_post            <- sum(delta_sm_post)
    delta_replicates_pre                <- c(delta_sm_pre)
    sum_delta_replicates_pre            <- sum(delta_sm_pre)
    success_rate_pre                  <- sum(delta_sm_pre)/TOT_VISITORS
    success_rate_post                  <- sum(delta_sm_post)/TOT_VISITORS
    names_states                    <-  paste("s_click", 1:K_FULL, sep="") 
    names_morphs                    <-  paste("m_click", 0:K_FULL, sep="")    
    
    # Bounce rates
    bounce1=length(bounced[bounced==1])/TOT_VISITORS
    bounce2=length(bounced[bounced==2])/TOT_VISITORS
    bounce3=length(bounced[bounced==3])/TOT_VISITORS
    bounce4=length(bounced[bounced==4])/TOT_VISITORS
    bounce7=length(bounced[bounced==7])/TOT_VISITORS 
    bounce15=length(bounced[bounced==15])/TOT_VISITORS 
   
# VII.Who recovers best? Who performs best? #############
   
  EMPIRICAL_SETTING
  mean(HULB_abs_error_replicates)
  mean(HMM_abs_error_replicates)
  mean(HMM_NATURE_abs_error_replicates)
 
  #  results for this run (see parameters above)
  QS_ARRIVAL
  QS_ARRIVAL_NATURE
  MODEL      
  TERMINAL  
  TOT_STATES
  sum_delta_replicates
  success_rate
  c(bounce1, bounce2,bounce3, bounce4,bounce7,bounce15)   
  alpha
  beta
  #RMSE
  RMSE_HULB_7
  RMSE_HMM_7
    
  # for table 6
  delta_sm
     
  # for table 7
    summarys
    summary_sm
      
   
   
     
# VIX. Write results to file   #############
   setwd(paste(PATH,"/_results",sep="") )
   write.csv(summary_sm3, "summary_sm3.csv")
   write.csv(summary_sm2, "summary_sm2.csv")
   write.csv(summary_sm1, "summary_sm1.csv")
   write.csv(delta_sm,"delta_sm")
   write.csv(summarys,"summarys.csv")
   write.csv(summary_sm, "summary_sm.csv")
   write.csv(t(sum_delta_replicates), "sum_delta_replicates.csv")
   write.csv(t(HULB_abs_error_replicates), "HULB_abs_error_replicates.csv")
   write.csv(t(HMM_abs_error_replicates), "HMM_abs_error_replicates.csv")
   write.csv(t(HMM_NATURE_abs_error_replicates), "HMM_NATURE_abs_error_replicates.csv")
   write.csv(delta_replicates, "delta_replicates.csv")
   write.csv(monitor_G_m1_s,   "monitor_last_G_m1_s.csv")
   write.csv(monitor_G_m2_s,   "monitor_last_G_m2_s.csv")
 #  write.csv(t(last_state_and_morph), "last_state_and_morph.csv", col.names = TRUE)
   save.image("Round2_Nov1_2019.RData")


 
   