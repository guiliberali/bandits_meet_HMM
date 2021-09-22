# Replication code for "Morphing for Consumer Dynamics: Bandits Meet HMM"
# Author: Gui Liberali and Alina Ferecatu 
# Date: October 2021
# Simulation code for Application 1 - MBA study

## Load required packages ----

getwd()



PATH_HMM_EST <- '../2. HMM Estimates'
PATH_HMM_EST_2S <- '../2. HMM Estimates/Estimates for 2 states'

registerDoParallel(cores=5) ## change depending on server capacity
  
  
  HMM_BANDIT <- "HMM_BANDIT" ;
  MAB_FIXED_STATE <- "MAB_FIXED_STATE"; 
  RANDOM <- "RANDOM"; ## baseline
  
## 1. Import data, parameters, functions ----
  # Key Frequently-Used Parameters  
  MODEL    <- MAB_FIXED_STATE # options are: HMM_BANDIT, MAB_FIXED_STATE, RANDOM
  TERMINAL <- TRUE # our chosen model
  
  TOT_STATES        <- 2  
  QS_ARRIVAL        <- if(TOT_STATES==2) {c(0.55, 0.45)}   else  {c(0.34, 0.33, 0.33) }
  QS_ARRIVAL_NATURE <- if(TOT_STATES==2) {c(0.98, 0.02)} else {c(.98,.015, 0.005)}  
  K_FULL            <- 4 # number of clicks in total - serve optimal morph after K_FULL clicks
      
  # Internal parameters 
  TOT_VISITORS    <- 40000 
  TOT_MORPHS      <- 2  
  TOT_PERIODS     <- 4
  STATES          <- 1:TOT_STATES # needed for the draw_state funtion
  XA_dim          <- 4
  if(TOT_STATES==2) {all_s <- cbind(1,c(-1,1)) }else{all_s=cbind(1, c(-1, 1, -1), c(-1, -1,1))} 
  INIT_STATES     <- TOT_STATES

  # Loading dock: data and functions
  #source(paste (PATH_CODE,                     "/Application1_Functions_for_Sims.R" , sep="") )   
  Gmatrix   <- read.table( DROPBOX_GMATRIX_LINK) 
  Ckjn      <- read.csv(paste (PATH_HMM_EST, file="/c_ij10.csv", sep="") , header=T) # import true ckjn: 10 pages for now

  mu_vec       <- read.csv(paste (PATH_HMM_EST_2S, file= "/mu_2Sterminal.csv",   sep=""),    header=T)
  rho_vec      <- read.csv(paste (PATH_HMM_EST_2S, file= "/rho_2Sterminal.csv",  sep=""),    header=T)
  pbounce_vec  <- read.csv( paste (PATH_HMM_EST_2S,file= "/pmf_bounce_user_2Sterminal.csv", sep="") , header=T)
  OOMEGA       <- read.csv(paste (PATH_HMM_EST_2S, file= "/omega_2Sterminal.csv", sep=""),    header=T)
  NEW_OMEGA       <- read.csv(paste (PATH_HMM_EST_2S, file= "/ArrivalState_omega_2Sterminal.csv",       sep="") , header=T)
  psm_exposure <- read.csv(paste (PATH_HMM_EST_2S, file= "/psm_2Sterminal.csv",   sep=""),    header=T)
 
  
  # Parameters 
  DIM_STATES=TOT_STATES-1
  rho       <- array(rho_vec[,2], dim=c(DIM_STATES,  XA_dim, TOT_MORPHS))  
  mu        <- array(mu_vec[,2], dim=c(DIM_STATES, (TOT_STATES-1), TOT_MORPHS)) 
  omega     <- matrix(OOMEGA[,2], ncol=TOT_STATES)
  new_omega <- matrix(NEW_OMEGA[,2], ncol=TOT_STATES)
  Ckjn      <- as.matrix(Ckjn%>% filter(page_id %in% 1:K_FULL)) 

  ## Transition probabilities at the last state K: 0/1 (terminal model)
  Q0_K      <- if(TOT_STATES==3) {array(rep(c(0,0,1), TOT_STATES), dim=c(1, TOT_STATES, TOT_MORPHS))} else
  {array(rep(c(0,1), TOT_STATES), dim=c(1, TOT_STATES, TOT_MORPHS))}
  
  # DP parameters
  TOT_LINKS <- nrow(Ckjn) ## number of links differs per page, thus we define the overall number of links and vectorize them per page
  psm_true  <- matrix(psm_exposure[,2],nrow=TOT_MORPHS, ncol=TOT_STATES, byrow=T) # emission probs 
  
  # Set the probability of bouncing
  p_bounce_click <- matrix(pbounce_vec[,2], nrow=TOT_STATES, byrow = T)
  
  # 2. Computes click probabilities, inspect them ----
  # This uses omega, because omega is based on actual states at every click, and nature is aware of that.
  p      <-  array(,c(TOT_STATES+1, TOT_LINKS, TOT_MORPHS)); p <-  Compute_click_prob(K_FULL)
  for (m in 1:TOT_MORPHS) {for (s in 1: TOT_STATES) { for (k in 1:K_FULL) { if(!round(sum(p[s+1, , m][p[1,,m]==k]),10) ==1) ERROR } } }  

 set.seed(9000)  
 trials=10

  # 3. Simulation: Loop over replicates and visitors -----
  ptime <- system.time({
      sim_1k_reps_2S_fixed <- foreach(icount(trials),.packages=c('expm','nnet', 'maotai', 'tidyverse')  ,.combine=rbind) %dopar% {
   # 3.1* Initialize data structures used for storage ----
   bounced           <- rep(K_FULL, TOT_VISITORS)
   I                 <- matrix(, ncol=TOT_MORPHS)
   delta_sm          <- matrix(0, nrow = TOT_STATES, ncol=TOT_MORPHS)
   colnames (delta_sm ) <-  c("abstract", "concrete")
   rownames (delta_sm ) <- if(TOT_STATES==2) {c("early", "late") } else  {c("early", "mid","late") } 
   I_evolution       <- matrix(0,nrow=TOT_VISITORS, ncol=TOT_MORPHS)
   morph_chosen      <- matrix(0,nrow=TOT_VISITORS, ncol=K_FULL+1)  
   true_state        <- rep(0, K_FULL)
   true_states_everyone <- matrix(0, nrow=TOT_VISITORS, ncol=K_FULL) 

   alphabeta <- array(rep(1, TOT_STATES*TOT_MORPHS*2), dim=c(TOT_STATES, TOT_MORPHS, 2))
   alpha     <- alphabeta[,,1] ; beta   <- alphabeta[,,2]
   G_current_sm   <- matrix(0,  nrow=TOT_STATES, ncol=TOT_MORPHS) 
   
   G_asym_m1      <- G_asym_m2    <-  matrix(0, nrow=TOT_VISITORS, ncol=TOT_STATES)  # current value of the  index
   monitor_G_m1_s <- monitor_G_m2_s <- matrix(0, nrow=TOT_VISITORS, ncol=TOT_STATES)

   # 3.2* Loop over visitors ----
   # ptime <- system.time({ ## use it to check run times
   # set.seed(9000)
   for (visitor  in 1:TOT_VISITORS)
      {
        # 1 Initialize clickstream, data structures to store state estimates
        Y                 <- rep(0, K_FULL);  Y_vec <-  NULL 
        qs_HMM_nature     <- qs_HMM    <- qs_FIXED_STATE   <- matrix(, nrow=K_FULL+1, ncol=TOT_STATES) # click-by-click state estimates 
        qs_HMM_nature[1,]           <- QS_ARRIVAL_NATURE
        qs_HMM[1,]  <- QS_ARRIVAL  
        qs_FIXED_STATE[1,] <- QS_ARRIVAL    

        # 2  Assign the initial morph - random in all scenarios
        ## Under RANDOM the morph does not change
        best_morph               <- which(rmultinom(1,1,matrix(c(1/TOT_MORPHS),ncol=TOT_MORPHS) ) == 1)  
        morph_chosen[visitor,1]  <- best_morph 
        
        # 3 Load current G (use asymptotic value if alpha or beta are greater than 3000) 
        for (st_ in 1: TOT_STATES) {  for (mor in 1:TOT_MORPHS)  { if(any(alpha[st_,mor] >= 3000, beta[st_,mor] >= 3000)) {
          G_current_sm[st_, mor] <- alpha[st_,mor]/(alpha[st_,mor]+beta[st_,mor])  } else {   
            G_current_sm[st_, mor] <- G_interpolation_fast(alpha[st_,mor],beta[st_,mor])  } } }
        
        # 3 Loop over all clicks
        BOUNCED <- FALSE # set flag as bot bouncing 
        for (click in 1:K_FULL)     
        {
          if (!BOUNCED) 
          {
           # A. DRAW TRUE STATE:  S(t)= f(qs_HMM(t) )
           true_state[click] <- Draw_state(STATES, qs_HMM_nature[click,]) 
          
           # B. DRAW CLICK: Click(t+1) = f(S(t)) # REMOVEd INITIALIZATION TEMP <- NULL
           temp     <- rmultinom (1,1,p[true_state[click]+1, ,best_morph][p[1,,best_morph]==click])
           Y_vec    <- c(Y_vec, temp)
           Y[click] <- which(temp==1)  
           Y_mat    <- cbind(Ckjn[,1][Ckjn[,1] %in% 1:click], Y_vec)
          
           #  HMM PRODUCTION: UPDATE HMM STATE  # variables used here:
           #  XA: Covariates of the transition matrix
           #  Lambda: transition matrices (one per morph)
          
           ## C. Get the covariates for the transition matrix: XAs are specific to the links a consumer clicked on  
           XA      <- Ckjn[Ckjn[,1]==click, 2:ncol(Ckjn)] ;  XA <- XA[Y[click],];   XA[1] <- log(XA[1])
          
           ## D. computes transition probabilities for previous state 1 to state K-1
           lambda  <- Transition_Matrix_study1(XA, rho, mu, DIM_STATES, TERMINAL) 
          
           ## E. stacks the transition matrix from state 1 to K-1, with the last row of the transition matrix, which is 0/1. Lambda is computed after every click
           if (TERMINAL) {lambda   <- array(c(rbind(lambda[,,1], Q0_K[,,1]), rbind(lambda[,,2], Q0_K[,,2])), dim=c(TOT_STATES, TOT_STATES, TOT_MORPHS))} 
         
           # F. Selects the transition matrix corresponding to the morph served and stores as Lambda_m
           lambda_m                <- lambda[,,best_morph]     
          
           # G. Nature knows the true state CHANGE!!!
           qs_HMM_nature[click+1,] <- lambda_m[true_state[click], ]  
          
           # H. Computes HMM estimate              
           if (MODEL==HMM_BANDIT) {   qs_HMM[click+1,]  <- qs_HMM[click+1,] <- Bayesian_Thm_updater_element(lambda, qs_HMM[click,],best_morph, psm_true)  }  # system does not know the true state so it uses lambda_m
           
           # I. UPDATE HULB STATE: qs_HULB= f(click (t+1), m(t)  ) 
           if (MODEL==MAB_FIXED_STATE) {qs_FIXED_STATE[click+1, ]     <- Bayesian_updater(click, best_morph, QS_ARRIVAL, new_omega, Ckjn, all_s, Y_mat )[click,] }
          
           # J. UPDATE FUTURE STATE: qs_future = f(click (t+1), m(t) , lambda )
           if (MODEL == HMM_BANDIT) { qs_HMM_sys_future <- estimate_future_period_qs_general(lambda)}
          
           # K. Computes EGI: m*(t=1) = Argmax EGI(qs_hulb(t+1)  )
           if (MODEL == MAB_FIXED_STATE)             {qs <- qs_FIXED_STATE[click, ] } 

           if (MODEL == MAB_FIXED_STATE)
           {
             for (m in 1:TOT_MORPHS) 
              {
                if (any(beta[ ,m]> 2999)) 
                 {
                   I[m] <- sum( (alpha[ ,m]/ (alpha[ ,m]+beta[ ,m])) * qs )
                  } else {
                    EGI  <- mapply(G_interpolation_fast, alpha[,m], beta[,m]) 
                    I[m] <- sum(EGI *qs)     
                  } # close if
               }# close for
               I_evolution[visitor,] <- I  
            } #close if loop over models        
          
           # L. Computes DP using HMM's m*(t=1) = Argmax DP(qs_HMM(t+1), qs_future(t+1) ) 
           if (MODEL == HMM_BANDIT) { best_morph_uncond <- DP(click, qs_HMM_sys_future, qs_HMM[click+1,], G_current_sm, p_bounce_click) }  # , delta_dp= delta_dp[period+1]) }
          
           # M. Assigns best morph according to the method of choice ( breaking ties randomly).
           if (MODEL == HMM_BANDIT)       { best_morph <- best_morph_uncond[click] }
           if (MODEL == MAB_FIXED_STATE)             { best_morph <- which.is.max(I[]) }   

           morph_chosen[visitor, click+1] <-best_morph
          
           # N. Bounce away or stay: bounced[visitor] will have the number of the last click 
           if (click <K_FULL)
           { temp_prob   <- p_bounce_click[true_state[click], morph_chosen[visitor, click]] 
             draw_bounce <- sample(c(1,0), 1, c(temp_prob, 1-temp_prob), replace=TRUE)
            if (draw_bounce == 1)
              {
               bounced[visitor] <- click
               BOUNCED <- TRUE # she is a bouncer, so get out of the inner loop of click, but continues processing this visitor to store data
              }
            }  # close link <4 test
          
           } # close if bounced check
          
        } # close clickstream loop 
        
        # 3.5 random exit (bouncing) so let's figure out the last click
        LAST_CLICK  <- bounced[visitor]
        
        # 3.6 Stores true state
        true_states_everyone[visitor, ]     <- true_state

        # 3.7 Draws purchase: delta_sm(t+1) = f( Psm(t), m(t), S(t) ) 
        psm_i <- psm_true[morph_chosen[visitor, LAST_CLICK], true_state[LAST_CLICK]] #true probability of purchase at every click using true state which changed at every click via HMM (unobserved by system)
        delta <-  rbinom(1, 1, psm_i)
        delta_sm[true_state[LAST_CLICK], morph_chosen[visitor, LAST_CLICK] ]<- delta_sm[true_state[LAST_CLICK], morph_chosen[visitor,LAST_CLICK] ] + delta 
     
        if (MODEL==HMM_BANDIT)  {qs_algo <- qs_HMM[LAST_CLICK,] }
        if (MODEL==MAB_FIXED_STATE) {qs_algo <- qs_FIXED_STATE[LAST_CLICK,] }  
        if (MODEL==RANDOM) {qs_algo <- rep(1, TOT_STATES) } # set just so that the loop can run
        
        # 3.8 Updates alpha, beta  = f( qs_FIXED_STATE(t+1) or f(qs_HMM(t+1) ) depending  on method being MAB_FIXED_STATE or HMM
        for (s in 1:TOT_STATES) {
          alpha[s,best_morph] <- alpha[s,best_morph] +      delta   * qs_algo[s] 
          beta[s, best_morph]  <- beta[s,best_morph] + (1 - delta) * qs_algo[s]       
        }

        # 3.9 monitor G  
        for (s in 1:TOT_STATES) {
          monitor_G_m1_s[visitor,s]  <- G_current_sm[s,1]
          monitor_G_m2_s[visitor,s]  <- G_current_sm[s,2]
          G_asym_m1[visitor,s]       <-  alpha[s,1]/(alpha[s,1]+beta[s,1])
          G_asym_m2[visitor,s]       <-  alpha[s,2]/(alpha[s,2]+beta[s,2])
        }
        # Reset bouncing flag 
        BOUNCED <- FALSE
        
      } # close visitor loop 
   
  # })[3]
   # ptime
   
   # 3.3* Gather results  ----
   # successes per state and morph.
   delta_replicates                <- c(delta_sm)
   sum_delta_replicates            <- sum(delta_sm)
   success_rate                    <- sum(delta_sm)/TOT_VISITORS
   # state evolution by click
   state_morph_combo=true_states_everyone %>%     cbind(morph_chosen)
   names_states                    <-  paste("s_click", 1:K_FULL, sep="")
   names_morphs                    <-  paste("m_click", 0:K_FULL, sep="")
   colnames(state_morph_combo)<-c(names_states,names_morphs)
   state_morph_combo=data.frame(state_morph_combo) %>% mutate(visitor=1:TOT_VISITORS, s_click2=ifelse(bounced<2, NA, s_click2), s_click3=ifelse(bounced<3, NA, s_click3), s_click4=ifelse(bounced<4, NA, s_click4),
             m_click2=ifelse(bounced<2, NA, m_click2), m_click3=ifelse(bounced<3, NA, m_click3), m_click4=ifelse(bounced<4, NA, m_click4)                         )

   summary_s1 = state_morph_combo %>% filter(!is.na(s_click1)) %>%   group_by(s_click1) %>%   summarise(count=n())
   summary_s2 = state_morph_combo %>% filter(!is.na(s_click2)) %>%   group_by(s_click2) %>%   summarise(count=n())
   summary_s3 = state_morph_combo %>% filter(!is.na(s_click3)) %>%   group_by(s_click3) %>%   summarise(count=n())
   summary_s4 = state_morph_combo %>% filter(!is.na(s_click4)) %>%   group_by(s_click4) %>%   summarise(count=n())
   summarys   = cbind(Click_1=summary_s1$count, Click_2=summary_s2$count,   Click_3=summary_s3$count, Click_4=summary_s4$count)
   summary_sm = state_morph_combo %>% filter(!is.na(s_click4) &!is.na(m_click4)) %>%  group_by(m_click4, s_click4) %>%  summarise(count=n())
   summary_statemorph=c(summary_sm$count)
   summary_state = c(summarys)
   bounce1=length(bounced[bounced==1])/TOT_VISITORS
   bounce2=length(bounced[bounced==2])/TOT_VISITORS
   bounce3=length(bounced[bounced==3])/TOT_VISITORS
   bounce4=length(bounced[bounced==4])/TOT_VISITORS

   Overall_results=c(delta_replicates, success_rate,  bounce1, bounce2, bounce3, bounce4, summary_state, summary_statemorph)

  } # close 1k reps loop
  
})[3]
