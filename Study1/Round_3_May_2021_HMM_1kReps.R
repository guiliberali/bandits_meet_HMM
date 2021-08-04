#################################################################################################################
#  Morphing Consumer Dynamics - Simulation Code                                                        
#  March 13, 2019
#  This file was coded from scratch by Gui in December 2018/January 2019 using HULB & functions + csv's from Alina
#
#  Steps taken to create this simulation code
#
#   1. Started with the original HULB code
# 
#   2. Adapted the original HULB code to reflect the core DP framework and solution developed by Gui in the 
#      Excel file "DP Feb 7 2018 revised.xlsx"  (emailed on February 8, 2018 and improved later by all). For 
#      this, Gui used the DP function that Alina wrote that implements the "DP Feb 7 2018 revised.xlsx" solution.
#
#   3. Created the flag for which run to do (HULB OR WM_HTM). One single code, two alternative methods to run.
#
#   4. Adapted the nature-simulation elements in this code to reflect dynamic states. Recall this was based 
#      on cognitive styles.
#    
#   5. Injected core HMM functions written by Alina:
#        - individualAij_terminal_march2019
#        - estimate_future_period_qs
#        - the more elegant Bayesian_updater 
#      
#   6. Coded HMM production, including using t-specific p_bouncing as the DP discount rate
#
#   7. Coded the handler that takes as input the following  .csv files remarkably estimated by Alina:
#        - mu_terminal.csv  from Alina's email of March 8, 2019
#        - rho_terminal.csv  from Alina's email of March 8, 2019
#        - q0_terminal.csv from Alina's email of March 9,2019
#        - mega_terminal.csv from Alina's email of March 9, 2019
#        - c_ij10.csv fom Alina's email of March 9, 2019
#        - p_bounce_long_terminal_user_levelMarch2019.csv from Alina's email of March 13, 2019
#
#  This code
#     - calls functions hosted in  'CJ_functions_2019.r' and 'Gmatrix.out'
#     - This is based on 4 clicks, terminal state, random exit
#
# NOTE BY GUI ON OCTOBER 9 2019:
#    THIS VERSION OF THE CODE HAS ADAPTATIONS BY ALINA BOTH HERE AND IN THE FUNCTION FILES. ADAPTATIONS WERE NOT
#    DOCUMENTED, SO I AM KEEPING THEM AND INSPECTING ONE BY ONE. JUST CHANGED THE FOLDERS & FILES SECTION AND COMMENTED OUT THE LOOP OVER TRIALS
#
# NOTE BY GUI ON OCT 17
#    Now uses Bayes thm from MONTOYA ET AL. (2010) eq. 12 IN PRODUCTION. It's the FUNCTION Bayesian_Thm_updater_element  
# About p
#   The first line is the page_id, which was necessary for subsetting. 
#   Turns out, subsetting arrays in R is a pain, and if you just do something like p[p[1,,]==1], it breaks down the dimensionality and it does not work. 
#   So you would have to first break the array into matrices, and then subset the matrix. Here is my solution, in the simulation code at line 231. 
#   p[true_state[click]+1, ,best_morph_next_t][p[1,,best_morph_next_t]==click]
#   So, to check all links of a page you would need to look morph by morph, and state by state. If you do :
#   p[true_state+1, , morph_id][p[1,,morph_id]==page_id], you will get the probabilities of clicking for state true_state, for morph_id, and for page_id (which is the same as the click number).
#   Since there are different number of links in each page, this was the only way I found to work with an array. Otherwise, I would have to make p a list, and that slower to work with.
# NOTE BY GUI ON OCT 24
#   Changed the code to correctly handle bouncing (the break was not working)
#   Changed the code to use LAST_CLICK instead of click in delta_sm.
# Key state variables
#    qs_system: click-by-click system estimate of state - depends on the WM_HMM flag  
#    qs_HMM: click-by-click HMM estimate of state 
#    qs_HMM_sys_future: HMM estimate of future qs's 
#
#  Note by Gui on Oct 28
#     - added flag and changed step M to run Oded's suggestion (an HMM only, without our DP)
#     - added the flag MODEL to make it easier to benchmark different models 
#
################################################################################################################

#  Morphing Consumer Dynamics 
#  May 13, 2021
## scp -r ~/Desktop/HMM_Bandit_codes/RSM_MBA_app/sims/codes/source ubuntu@145.100.58.221:/data/morphing/RSM_MBA_app 
## scp -r ubuntu@145.100.58.221:/data/morphing/RSM_MBA_app/source/results ~/Desktop/HMM_Bandit_codes/RSM_MBA_app/sims

  rm(list=ls()) 
# clean up R envirnoment
  set.seed(9000)
  
  library(doParallel)
  registerDoParallel(cores=6)
  library(expm) ## for matrix power
  library(nnet) # to break ties randomly
  library(tidyverse)
  library(reshape2)
  library(magrittr)
  

  server=F
  # Paths
  if (server)
  {
    PATH<- "/data/morphing/RSM_MBA_app"
    PATH_DATA <- paste0(PATH, "/source") 
    PATH_CODE <- paste0(PATH, "/source/hmmBandit_functions" )
    PATH_PARAMETERS   <- paste0(PATH, "/source/rsm_sim" )
    PATH_RESULTS   <- paste0(PATH_DATA, "/results" )
  }  else { 
    PATH          <-  "~/Desktop/HMM_Bandit_codes/RSM_MBA_app/sims/codes" 
    PATH_DATA <- paste0(PATH, "/source") 
    PATH_CODE <- paste0(PATH, "/source/hmmBandit_functions" )
    PATH_PARAMETERS   <- paste0(PATH, "/source/rsm_sim" )
    PATH_RESULTS   <- "~/Desktop/HMM_Bandit_codes/RSM_MBA_app/sims/results"
  } 
  
  HMM_BANDIT <- "HMM_BANDIT" ;   HULB <- "HULB" ;   HMM_MAB_NO_DP <- "HMM_MAB_NO_DP"   # for simplicity & clarity
  FIXED_STATE <- "FIXED_STATE"; RANDOM <- "RANDOM";  MAB_ONLY <- "MAB_ONLY"
  
  ##############################################################################
  # Key Frequently-Used Parameters  
  
  MODEL    <- RANDOM # options are: HMM_BANDIT, HULB, HMM_MAB_NO_DP  BTW, no need for ""
  TERMINAL <- TRUE # our chosen model: false
  
  # less-used parameters
  TOT_STATES        <- 2  # our chosen model: 3
  QS_ARRIVAL        <- if(TOT_STATES==2) {c(0.55, 0.45)}   else  {c(0.34, 0.33, 0.33) } # c{1/3,1/3,1/3) }  or {c(6/11, 6/22, 6/33)}  
  QS_ARRIVAL_NATURE <- if(TOT_STATES==2) {c(0.98, 0.02)} else {c(.98,.015, 0.005)}  
  K_FULL            <- 4 # number of clicks in total - serve optimal morph after K_FULL clicks
      
  # Internal parameters, rarely change  
  TOT_VISITORS    <- 40000 
  TOT_MORPHS      <- 2  
  TOT_PERIODS     <- 4
  STATES          <- 1:TOT_STATES # needed for the draw_state funtion
  XA_dim          <- 4
  if(TOT_STATES==2) {all_s <- cbind(1,c(-1,1)) }else{all_s=cbind(1, c(-1, 1, -1), c(-1, -1,1))} 
  INIT_STATES     <- TOT_STATES

  # Loading dock: data and functions
  source(paste (PATH_CODE,                     "/Round_2_functions_Oct_22_2019.R" , sep="") )   
  Gmatrix   <- read.table(paste (PATH_DATA,    "/Gmatrix.out", sep="") )   
  Ckjn      <- read.csv(paste (PATH_DATA, file="/c_ij10.csv", sep="") , header=T) # import true ckjn: 10 pages for now
  
  if (TOT_STATES==2){
    if (TERMINAL)
    {
      mu_vec       <- read.csv(paste (PATH_PARAMETERS, file= "/2S_terminal/mu_2Sterminal.csv",   sep=""),    header=T)
      rho_vec      <- read.csv(paste (PATH_PARAMETERS, file= "/2S_terminal/rho_2Sterminal.csv",  sep=""),    header=T)
      pbounce_vec  <- read.csv( paste (PATH_PARAMETERS,file= "/2S_terminal/pmf_bouce_user_2Sterminal.csv", sep="") , header=T)
      OOMEGA       <- read.csv(paste (PATH_PARAMETERS, file= "/2S_terminal/omega_2Sterminal.csv", sep=""),    header=T)
      NEW_OMEGA       <- read.csv(paste (PATH_PARAMETERS, file= "/2S_terminal/ArrivalState_omega_2Sterminal.csv",       sep="") , header=T)
      psm_exposure <- read.csv(paste (PATH_PARAMETERS, file= "/2S_terminal/psm_2Sterminal.csv",   sep=""),    header=T)
    } else{
      mu_vec       <- read.csv(paste (PATH_PARAMETERS, file= "/2S_nonterminal/mu_2Snonterminal.csv",   sep=""),   header=T)
      rho_vec      <- read.csv(paste (PATH_PARAMETERS, file= "/2S_nonterminal/rho_2Snonterminal.csv",  sep=""),   header=T)
      pbounce_vec  <- read.csv(paste (PATH_PARAMETERS, file= "/2S_nonterminal/pmf_bounce_user_2Snonterminal.csv", sep="") , header=T)
      OOMEGA       <- read.csv(paste (PATH_PARAMETERS, file= "/2S_nonterminal/omega_2Snonterminal.csv", sep=""),   header=T)
      psm_exposure <- read.csv(paste (PATH_PARAMETERS, file= "/2S_nonterminal/psm_2Snonterminal.csv",   sep=""), header=T)
    } 
  }else{
    if (TERMINAL)
    {
      mu_vec       <- read.csv(paste (PATH_PARAMETERS, file= "/3S_terminal/mu_3Sterminal.csv",          sep="") , header=T)
      rho_vec      <- read.csv(paste (PATH_PARAMETERS, file= "/3S_terminal/rho_3Sterminal.csv",         sep="") , header=T)
      pbounce_vec  <- read.csv( paste (PATH_PARAMETERS,file= "/3S_terminal/pmf_bouce_user_3Sterminal.csv", sep="") , header=T)
      OOMEGA       <- read.csv(paste (PATH_PARAMETERS, file= "/3S_terminal/omega_3Sterminal.csv",       sep="") , header=T)
      NEW_OMEGA       <- read.csv(paste (PATH_PARAMETERS, file= "/3S_terminal/ArrivalState_omega_3Sterminal.csv",       sep="") , header=T)
      psm_exposure <- read.csv(paste (PATH_PARAMETERS, file= "/3S_terminal/psm_3Sterminal.csv",         sep="") , header=T)
    } else{
      mu_vec       <- read.csv(paste (PATH_PARAMETERS, file= "/3S_nonterminal/mu_3Snonterminal.csv",     sep="") , header=T)
      rho_vec      <- read.csv(paste (PATH_PARAMETERS, file= "/3S_nonterminal/rho_3Snonterminal.csv",    sep="") , header=T)
      pbounce_vec  <- read.csv(paste (PATH_PARAMETERS, file= "/3S_nonterminal/pmf_bouce_user_3Snonterminal.csv", sep="") , header=T)
      OOMEGA       <- read.csv(paste (PATH_PARAMETERS, file= "/3S_nonterminal/omega_3Snonterminal.csv",  sep="") , header=T)
      psm_exposure <- read.csv(paste (PATH_PARAMETERS, file= "/3S_nonterminal/psm_3Snonterminal.csv",    sep="") , header=T)
    } 
  }
  
  # Parameters used in this terminal state version
  if (TERMINAL) {DIM_STATES=TOT_STATES-1}else{DIM_STATES=TOT_STATES}
  rho       <- array(rho_vec[,2], dim=c(DIM_STATES,  XA_dim, TOT_MORPHS))  
  mu        <- array(mu_vec[,2], dim=c(DIM_STATES, (TOT_STATES-1), TOT_MORPHS)) 
  omega     <- matrix(OOMEGA[,2], ncol=TOT_STATES)
  new_omega <- matrix(NEW_OMEGA[,2], ncol=TOT_STATES)
  Ckjn      <- as.matrix(Ckjn%>% filter(page_id %in% 1:K_FULL)) 
 
  #  XA: Covariates of the transition matrix
  XA        <- c(log(6), 0, 0, 0)
  ## Transition probabilities at the last state K: 0/1, because we re using terminal states.
  Q0_K      <- if(TOT_STATES==3) {array(rep(c(0,0,1), TOT_STATES), dim=c(1, TOT_STATES, TOT_MORPHS))} else
  {array(rep(c(0,1), TOT_STATES), dim=c(1, TOT_STATES, TOT_MORPHS))}
  
  
  ## Q0 binds the transition matrix for previous state 1 to K-1, with the last row of the transition matrix, which is 0/1.
  if (TERMINAL) {Q0 <- array(c(rbind(Transition_Matrix(XA, rho, mu, DIM_STATES, TERMINAL)[,,1], Q0_K[,,1]),rbind(Transition_Matrix(XA, rho, mu, DIM_STATES, TERMINAL)[,,2], Q0_K[,,2])), dim=c(TOT_STATES, TOT_STATES, TOT_MORPHS))} else
  {Q0 <- array(Transition_Matrix(XA, rho, mu, DIM_STATES, TERMINAL), dim=c(TOT_STATES, TOT_STATES, TOT_MORPHS))}# Q0 <- array(Q0_vec[,2], dim=c(INIT_STATES, TOT_STATES, TOT_MORPHS))
  
  # DP parameters
  TOT_LINKS <- nrow(Ckjn) ## number of links differs per page, thus we define the overall number of links and vectorize them per page
  psm_true  <- matrix(psm_exposure[,2],nrow=TOT_MORPHS, ncol=TOT_STATES, byrow=T) # emission probs psm was psm_exposure_newhack in non-terminal world
  
  # 1. Compute the probability of bouncing, i.e the switching cost FORMAT: p_bounce_mt[period,init_stage, morph]
  # pbounce_vec4 <- pbounce_vec[--which(pbounce_vec$click<= (TOT_PERIODS+1)),]
  # p_bounce_mt  <- array(rep(0, (TOT_PERIODS+1)*TOT_STATES*TOT_MORPHS),dim=c(TOT_PERIODS+1, TOT_STATES, TOT_MORPHS)) 
  # for (tt in 1:(TOT_PERIODS+1) ) {p_bounce_mt[tt,,]<- pbounce_vec4$value[which(pbounce_vec4$click==tt)]} # rows are states, columns are morphs
  # p_unbounce_mt <- 1-p_bounce_mt
  p_bounce_click <- matrix(pbounce_vec[,2], nrow=TOT_STATES, byrow = T)
  
  # 2. Computes click probabilities, inspect them
  p      <-  array(,c(TOT_STATES+1, TOT_LINKS, TOT_MORPHS)); p <-  Compute_click_prob(K_FULL)
  for (m in 1:TOT_MORPHS) {for (s in 1: TOT_STATES) { for (k in 1:K_FULL) { if(!round(sum(p[s+1, , m][p[1,,m]==k]),10) ==1) ERROR } } }  

set.seed(9000)  
 trials=12
 
 # Loop over replicates -----
  ptime <- system.time({
      sim_1k_reps <- foreach(icount(trials), .combine=rbind) %dopar% {
     #Initialize data structures used for storage
   N                 <- matrix(0, ncol=TOT_VISITORS)
   bounced           <- rep(K_FULL, TOT_VISITORS) #matrix(K_FULL, ncol=TOT_VISITORS)
   I                 <- matrix(, ncol=TOT_MORPHS)
   delta_sm          <- matrix(0, nrow = TOT_STATES, ncol=TOT_MORPHS)
   colnames (delta_sm ) <-  c("abstract", "concrete")
   rownames (delta_sm ) <- if(TOT_STATES==2) {c("early", "late") } else  {c("early", "mid","late") } 
   I_evolution       <- matrix(0,nrow=TOT_VISITORS, ncol=TOT_MORPHS)
   morph_chosen      <- matrix(0,nrow=TOT_VISITORS, ncol=K_FULL+1)  
   posterior         <- matrix( , nrow=TOT_VISITORS, ncol=K_FULL)  
   inspect_psm       <- inspect_delta      <- NULL
   true_state        <- rep(0, K_FULL)
   true_states_everyone <-qs_HULB_everyone <- qs_HMM_everyone <- 
   qs_HMM_nature_everyone   <- matrix(0, nrow=TOT_VISITORS, ncol=K_FULL) 
   error_state       <- matrix(0, nrow=TOT_VISITORS, ncol= 9); colnames(error_state)<- c("visitor", "qr_t1","qr_t2","qr_t3","qr_t4","qs_t1","qs_t2","qs_t3","qs_t4")
      
   alphabeta <- array(rep(1, TOT_STATES*TOT_MORPHS*2), dim=c(TOT_STATES, TOT_MORPHS, 2))
   alpha     <- alphabeta[,,1] ; beta   <- alphabeta[,,2]
   G_current_sm   <- matrix(0,  nrow=TOT_STATES, ncol=TOT_MORPHS) 
   
   #if(MODEL==MAB_ONLY) G_current_sm   <- matrix(0,    ncol=TOT_MORPHS)
   G_asym_m1      <- G_asym_m2    <-  matrix(0, nrow=TOT_VISITORS, ncol=TOT_STATES)  # current value of the  index
   monitor_G_m1_s <- monitor_G_m2_s <- matrix(0, nrow=TOT_VISITORS, ncol=TOT_STATES)

   # 2. Loop over all visitors
   # visitor      <-  1  # for testing purposes only 
   # ptime <- system.time({
   for (visitor  in 1:TOT_VISITORS)
      {
        # 3.1 Initialize clickstream, data structures to store state estimates (qs_system is HULB, qs_HMM is HMM) 
        Y                 <- rep(0, K_FULL);  Y_vec <-  NULL 
        qs_HMM_nature     <- qs_HMM    <- qs_HULB   <- matrix(, nrow=K_FULL+1, ncol=TOT_STATES) # click-by-click state estimates 
        qs_HMM_nature[1,]           <- QS_ARRIVAL_NATURE
        qs_HMM[1,]  <- QS_ARRIVAL  
        qs_HULB[1,] <- QS_ARRIVAL    

        # 3.2  Assign the initial morph - random in all scenarios
          ## Under RANDOM the morph does not change
        if (MODEL == FIXED_STATE)      {qs <- QS_ARRIVAL_NATURE } # state doesn't change throughout the visit
        if (MODEL == MAB_ONLY)         {qs <- rep(1, TOT_STATES)} # no fractional updating, just MAB
        
        if ( MODEL==FIXED_STATE )  
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
          best_morph <- which.is.max(I[])
          
        }else  if ( MODEL == MAB_ONLY)   {
          
          for (m in 1:TOT_MORPHS) 
          {
            if (beta[ 1,m]> 2999)
            {
              I[m] <- (alpha[1 ,m]/ (alpha[1 ,m]+beta[1 ,m]))
            } else {
              EGI  <- G_interpolation_fast(alpha[1,m], beta[1,m] )
              I[m] <-EGI     
            } # close if
          }# close for
          I_evolution[visitor,] <- I  
          best_morph <- which.is.max(I[])
          
        }else{ #close if loop over models       
        best_morph               <- which(rmultinom(1,1,matrix(c(1/TOT_MORPHS),ncol=TOT_MORPHS) ) == 1)  
        }
        
        morph_chosen[visitor,1]  <- best_morph 
        
        # 3.3 Load current G (use asymptotical value if alpha or beta are greater than 3000) 
        for (st_ in 1: TOT_STATES) {  for (mor in 1:TOT_MORPHS)  { if(any(alpha[st_,mor] >= 3000, beta[st_,mor] >= 3000)) {
          G_current_sm[st_, mor] <- alpha[st_,mor]/(alpha[st_,mor]+beta[st_,mor])  } else {   
            G_current_sm[st_, mor] <- G_interpolation_fast(alpha[st_,mor],beta[st_,mor])  } } }
        
        ## if (MODEL==FIXED_STATE) do I even need to draw a state for system????
        
        # 3.4 Loop over all clicsks
        # click =1 # for testing purposes only  
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
           #  Lambda: transition matrices (one per morph so this is current state x  state x morph)
          
           ## C. Get the covariates for the transition matrix: XAs are specific to the links a consumer clicked on  
           XA      <- Ckjn[Ckjn[,1]==click, 2:ncol(Ckjn)] ;  XA <- XA[Y[click],];   XA[1] <- log(XA[1])
          
           ## D. computes transition probabilities for previous state 1 to state K-1/ in our case we only have p11 and p12. Because if one starts in previous state 2, she stays there. So p21=0 and p22=1, this makes Q0_K=c(0,1)
           lambda  <- Transition_Matrix(XA, rho, mu, DIM_STATES, TERMINAL) 
          
           ## E. stacks the transition matrix from state 1 to K-1, with the last row of the transition matrix, which is 0/1. Lambda is computed after every click
           if (TERMINAL) {lambda   <- array(c(rbind(lambda[,,1], Q0_K[,,1]), rbind(lambda[,,2], Q0_K[,,2])), dim=c(TOT_STATES, TOT_STATES, TOT_MORPHS))} 
         
           # F. Selects the transition matrix corresponding to the morph served and stores as Lambda_m
           lambda_m                <- lambda[,,best_morph]     
          
           # G. Nature knows the true state CHANGE!!!
           qs_HMM_nature[click+1,] <- lambda_m[true_state[click], ]  
          
           # H. Computes HMM estimate              
           #if (HMM_BANDIT) {   qs_HMM[click+1,]  <- qs_HMM[click,] %*% lambda_m }  # system does not know the true state so it uses lambda_m
           if (is.element(MODEL, c(HMM_BANDIT, HMM_MAB_NO_DP) )) {qs_HMM[click+1,] <- Bayesian_Thm_updater_element(lambda, qs_HMM[click,],best_morph, psm_true) }   
          
           # I. UPDATE HULB STATE: qs_HULB= f(click (t+1), m(t)  ) 
           if (MODEL==HULB) {qs_HULB[click+1, ]     <- Bayesian_updater(click, best_morph, QS_ARRIVAL, new_omega, Ckjn, all_s, Y_mat )[click,] }
          
           # J. UPDATE FUTURE STATE: qs_future = f(click (t+1), m(t) , lambda )
           if (MODEL == HMM_BANDIT) { qs_HMM_sys_future <- estimate_future_period_qs_general(lambda)} #, qs_HMM[click+1,] ) } # Gui replaced lambda instead of q0
          
           # K. Computes EGI: m*(t=1) = Argmax EGI(qs_hulb(t+1)  )
           if (MODEL == HMM_MAB_NO_DP)    {qs <- qs_HMM[click,  ] } 
           if (MODEL == HULB)             {qs <- qs_HULB[click, ] } 
           # if (MODEL == FIXED_STATE)      {qs <- QS_ARRIVAL_NATURE } # state doesn't change throughout the visit
           # if (MODEL == MAB_ONLY)         {qs <- rep(1, TOT_STATES)} # no fractional updating, just MAB
           
           if ( is.element(MODEL, c(HULB, HMM_MAB_NO_DP) )    )     #FIXED_STATE, 
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
          
           ## for MAB only, there is no state, use only the first alpha/beta 
           #if ( MODEL == MAB_ONLY)    
           #{
           #  for (m in 1:TOT_MORPHS) 
           #  {
           #    if (beta[ ,m]> 2999)
           #    {
           #      I[m] <- (alpha[1 ,m]/ (alpha[1 ,m]+beta[1 ,m]))
           #    } else {
           #      EGI  <- G_interpolation_fast(alpha[1,m], beta[1,m] )
           #      I[m] <-EGI     
           #    } # close if
           #  }# close for
           #  I_evolution[visitor,] <- I  
           #} #close if loop over models        
           
           # L. Computes DP using HMM's m*(t=1) = Argmax DP(qs_HMM(t+1), qs_future(t+1) ) 
           ## TESTING  - changed qs_HMM[click+1,] to qs_HMM[click,]
           ## In the DP, changed exponent in (1-pbounce)^(next_period-K) into (period-K)
           if (MODEL == HMM_BANDIT) { best_morph_uncond <- DP(click, qs_HMM_sys_future, qs_HMM[click+1,], G_current_sm, p_bounce_click) }  # , delta_dp= delta_dp[period+1]) }
          
           # M. Assigns best morph according to the method of choice ( breaking ties randomly).
           if (MODEL == HMM_BANDIT)       { best_morph <- best_morph_uncond[click] }  # note that index should be click, not TOT_PERIODS
           if (MODEL == HMM_MAB_NO_DP)    { best_morph <- which.is.max(I[]) }     
           if (MODEL == HULB)             { best_morph <- which.is.max(I[]) }   
           #if (MODEL == FIXED_STATE)      { best_morph <- which.is.max(I[]) }
           #if (MODEL == MAB_ONLY)         { best_morph <- which.is.max(I[]) }
           
           morph_chosen[visitor, click+1] <-best_morph
          
           # N. Bounce away or stay: bounced[visitor] will have the number of the last click 
           if (click <K_FULL)
           { temp_prob   <- p_bounce_click[true_state[click], morph_chosen[visitor, click]] # p_bounce_mt[click, true_state[click], morph_chosen[visitor, click]]
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
        qs_HULB_everyone[visitor,  ]       <- qs_HULB[LAST_CLICK,1]
        qs_HMM_everyone[visitor,  ]        <- qs_HMM[LAST_CLICK,1]
        qs_HMM_nature_everyone[visitor,  ] <- qs_HMM_nature[LAST_CLICK,1]
        
        # 3.7 Draws purchase: delta_sm(t+1) = f( Psm(t), m(t), S(t) ) 
        psm_i <- psm_true[morph_chosen[visitor, (LAST_CLICK)], true_state[LAST_CLICK]] #true probability of purchase at every click using true state which changed at every click via HMM (unobserved by system)
        delta <-  rbinom(1, 1, psm_i)
        delta_sm[true_state[LAST_CLICK], morph_chosen[visitor, LAST_CLICK] ]<- delta_sm[true_state[LAST_CLICK], morph_chosen[visitor,LAST_CLICK] ] + delta # changed to LAST_CLICK
     
        if (is.element(MODEL, c(HMM_BANDIT, HMM_MAB_NO_DP)))  {qs_algo <- qs_HMM[LAST_CLICK,] }
        if (MODEL==HULB) {qs_algo <- qs_HULB[LAST_CLICK,] }    # set which state estimate to use
        if (MODEL == FIXED_STATE) {qs_algo <- QS_ARRIVAL_NATURE } 
        if (is.element(MODEL, c(RANDOM, MAB_ONLY)) ) {qs_algo <- rep(1, TOT_STATES) } 
        ## for MAB only, there is no state, so there is no fractional updating
        # for random it does not matter, qs_algo is there so that the loop can run
        
        # 3.8 Updates alpha, beta  = f( qs_hulb(t+1) or f(qs_HMM(t+1) ) depending  on method being HULB or HMM
        for (s in 1:TOT_STATES) {
          alpha[s,best_morph] <- alpha[s,best_morph] +      delta   * qs_algo[s] 
          beta[s, best_morph]  <- beta[s,best_morph] + (1 - delta) * qs_algo[s]       
        }
        #inspect_psm <- c(inspect_psm, psm_i);  inspect_delta <- c(inspect_delta, delta)  # stores qr at last click  
        
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
   
   ## FULL RESULTS
   # 4. Results - successes per state and morph.
   # 4.1 check last-state membership
   last_state <-last_morph <- last_qs_HULB <- last_qs_HMM <- last_qs_nature <-  true_state_binary <-  matrix (0, ncol = TOT_VISITORS)
   for (iii in 1: TOT_VISITORS )
      {  last_state[iii]     <- true_states_everyone[iii,bounced[iii]]
         last_morph[iii]     <- morph_chosen[iii,bounced[iii]]
         last_qs_HULB[iii]   <- qs_HULB_everyone[iii,bounced[iii]]
         last_qs_HMM[iii]    <- qs_HMM_everyone[iii,bounced[iii]]
         last_qs_nature[iii] <- qs_HMM_nature_everyone[iii,bounced[iii]]
         if(last_state[iii]==1) { true_state_binary[iii] <- 1} else { true_state_binary[iii] <- 0}
      }
   last_state_and_morph           <- t(rbind(last_state, last_morph, last_qs_HULB, last_qs_HMM, last_qs_nature ,true_state_binary ))
   colnames(last_state_and_morph) <- c("state", "morph", "qs_HULB_sl", "qs_HMM_sl", "qs_hMM_nature_sl", "last_true_state_binary")

   # 4.2 Store replicate results
   test                            <- last_state_and_morph
   HULB_abs_error_replicates       <- sum(abs(test[,3]-test[,6]))
   HMM_abs_error_replicates        <- sum(abs(test[,4]-test[,6]))
   HMM_NATURE_abs_error_replicates <- sum(abs(test[,5]-test[,6]))
   delta_replicates                <- c(delta_sm)
   sum_delta_replicates            <- sum(delta_sm)
   success_rate                    <- sum(delta_sm)/TOT_VISITORS
   state_morph_combo=true_states_everyone %>%     cbind(morph_chosen)
   names_states                    <-  paste("s_click", 1:K_FULL, sep="")
   names_morphs                    <-  paste("m_click", 0:K_FULL, sep="")
   colnames(state_morph_combo)<-c(names_states,names_morphs)
   state_morph_combo=data.frame(state_morph_combo) %>% mutate(visitor=1:TOT_VISITORS, s_click2=ifelse(bounced<2, NA, s_click2), s_click3=ifelse(bounced<3, NA, s_click3), s_click4=ifelse(bounced<4, NA, s_click4),
             m_click2=ifelse(bounced<2, NA, m_click2), m_click3=ifelse(bounced<3, NA, m_click3), m_click4=ifelse(bounced<4, NA, m_click4)                         )

   # 4.3 state evolution by click
   summary_s1 = state_morph_combo %>% filter(!is.na(s_click1)) %>%   group_by(s_click1) %>%   summarise(count=n())
   summary_s2 = state_morph_combo %>% filter(!is.na(s_click2)) %>%   group_by(s_click2) %>%   summarise(count=n())
   summary_s3 = state_morph_combo %>% filter(!is.na(s_click3)) %>%   group_by(s_click3) %>%   summarise(count=n())
   summary_s4 = state_morph_combo %>% filter(!is.na(s_click4)) %>%   group_by(s_click4) %>%   summarise(count=n())
   summarys   = cbind(Click_1=summary_s1$count, Click_2=summary_s2$count,   Click_3=summary_s3$count, Click_4=summary_s4$count)
   summary_sm = state_morph_combo %>% filter(!is.na(s_click4) &!is.na(m_click4)) %>%  group_by(m_click4, s_click4) %>%  summarise(count=n())
   summary_sm3 = state_morph_combo %>% filter(!is.na(s_click3) &!is.na(m_click3)) %>%  group_by(m_click3, s_click3) %>%  summarise(count=n())
   summary_sm2 = state_morph_combo %>% filter(!is.na(s_click2) &!is.na(m_click2)) %>%  group_by(m_click2, s_click2) %>%  summarise(count=n())
   summary_sm1 = state_morph_combo %>% filter(!is.na(s_click1) &!is.na(m_click1)) %>%  group_by(m_click1, s_click1) %>%  summarise(count=n())
   summary_statemorph=c(summary_sm$count)
   summary_state = c(summarys)
   bounce1=length(bounced[bounced==1])/TOT_VISITORS
   bounce2=length(bounced[bounced==2])/TOT_VISITORS
   bounce3=length(bounced[bounced==3])/TOT_VISITORS
   bounce4=length(bounced[bounced==4])/TOT_VISITORS

   Overall_results=c(delta_replicates, success_rate,  bounce1, bounce2, bounce3, bounce4, summary_state, summary_statemorph  )

   ### SUCCESS RATES ONLY
   #success_rate = sum(delta_sm)/TOT_VISITORS
   # c(success_rate, c(delta_sm))
  } # close 1k reps loop
  
})[3]
  
  ptime
  
  ## Results ----
  
# load("sim_fixed_state_4K_3T_1kRun.RData")
  ## how many people are moved to the late state at clicks 1, 2, 3, 4
  ## for two state nature
  summarys=colMeans(sim_1k_reps[,c(11, 13, 15, 17)])
  summarys/400

  ## for three state nature
  ## state 2
  summarys=colMeans(sim_1k_reps[,c(13, 16, 19, 22)])
  summarys/400
  ## state 3
  summarys=colMeans(sim_1k_reps[,c(14, 17, 20, 23)])
  summarys/400
  
  delta_sm_average = colMeans(sim_1k_reps[,1:4])
  delta_replicates=rowSums(sim_1k_reps[,1:4])
  
  ## CHANGE to [,1:6] for the 3T
  sim_1k_reps[1:2,1:17]
  hist(delta_replicates)
  median(delta_replicates)
  c(mean(delta_replicates), 
    t.test(delta_replicates)$conf.int)
  
  ## how many people are served m2 at last click 
  ## for two state nature
  summary_m2=mean(rowSums(sim_1k_reps[,18:21]))
  summary_m2 
  # HMM bandit: 22703 out of 24013
  # HULB: 19120.8 out of 24025
  
  # setwd("~/Desktop/HMM_Bandit_codes/RSM_MBA_app/sims/results")
  setwd(PATH_RESULTS)
  # save.image("sim_random_1kRun.RData")
   
  save.image("sim_random_4K_2T_1kRun.RData")
   
  ## monitor G for m1 vs m2, in state 1
  morph_chosen[1:5,]
 table(morph_chosen[,1])
 table(morph_chosen[,5])
 
 if(TOT_STATES==2) {
    colnames(monitor_G_m1_s)<- colnames(monitor_G_m2_s)<-c("State 1", "State 2")
  monitor_G=as_tibble(monitor_G_m1_s) %>% 
    mutate(Visitor=1:TOT_VISITORS) %>% 
    pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
    mutate(Morph=c("Morph 1")) %>% 
    bind_rows(as_tibble(monitor_G_m2_s) %>% 
                mutate(Visitor=1:TOT_VISITORS) %>% 
                pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
                mutate(Morph=c("Morph 2")) ) %>% 
  mutate(state_morph=case_when(State== 'State 1' & Morph== "Morph 1" ~ "S1M1",
                               State== 'State 2' & Morph== "Morph 1" ~ "S2M1",
                               State== 'State 1' & Morph== "Morph 2" ~ "S1M2",
                               TRUE ~ "S2M2"))
  }else{
    colnames(monitor_G_m1_s)<- colnames(monitor_G_m2_s)<-c("State 1", "State 2", "State 3")
  monitor_G=as_tibble(monitor_G_m1_s) %>% 
    mutate(Visitor=1:TOT_VISITORS) %>% 
    pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
    mutate(Morph=c("Morph 1")) %>% 
    bind_rows(as_tibble(monitor_G_m2_s) %>% 
                mutate(Visitor=1:TOT_VISITORS) %>% 
                pivot_longer(-Visitor, names_to="State", values_to="Index") %>% 
                mutate(Morph=c("Morph 2")) ) %>% 
    mutate(state_morph=case_when(State== 'State 1' & Morph== "Morph 1" ~ "S1M1",
                                 State== 'State 2' & Morph== "Morph 1" ~ "S2M1",
                                 State== 'State 3' & Morph== "Morph 1" ~ "S3M1",
                                 State== 'State 1' & Morph== "Morph 2" ~ "S1M2",
                                 State== 'State 2' & Morph== "Morph 2" ~ "S2M2",
                                 TRUE ~ "S3M2"))
}  
  monitor_G %>% 
    #filter(State == 'State 3') %>% 
    #filter(state_morph %in% c("S1M1", "S1M2")) %>% 
    filter(Visitor %in% 400:40000) %>% 
    ggplot(aes(x=Visitor, y=Index, color=state_morph))+
    geom_line()
  

   # 5. who recovers best? Who performs best?
   mean(HULB_abs_error_replicates)
   mean(HMM_abs_error_replicates)
   mean(HMM_NATURE_abs_error_replicates)
   # sum(delta_replicates )/TOT_REPLICATES# got to summarize this


   #     results for this run (see paramters above)
   QS_ARRIVAL
   QS_ARRIVAL_NATURE
   MODEL
   TERMINAL
   TOT_STATES
   sum_delta_replicates
   success_rate

   # for table 6
   delta_sm

   # for table 7
   summarys
   summary_sm









   # write results to file
   setwd(PATH)
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
   write.csv(t(last_state_and_morph), "last_state_and_morph.csv", col.names = TRUE)
   save.image("Round2_Nov1_2019.RData")






   ### check computation times ----

   # 2. Loop over all visitors
   # visitor      <-  1  # for testing purposes only 
   
   MODEL=RANDOM
   TOT_VISITORS=1
   
   #Initialize data structures used for storage
   N                 <- matrix(0, ncol=TOT_VISITORS)
   bounced           <- rep(K_FULL, TOT_VISITORS) #matrix(K_FULL, ncol=TOT_VISITORS)
   I                 <- matrix(, ncol=TOT_MORPHS)
   delta_sm          <- matrix(0, nrow = TOT_STATES, ncol=TOT_MORPHS)
   colnames (delta_sm ) <-  c("abstract", "concrete")
   rownames (delta_sm ) <- if(TOT_STATES==2) {c("early", "late") } else  {c("early", "mid","late") } 
   I_evolution       <- matrix(0,nrow=TOT_VISITORS, ncol=TOT_MORPHS)
   morph_chosen      <- matrix(0,nrow=TOT_VISITORS, ncol=K_FULL+1)  
   posterior         <- matrix( , nrow=TOT_VISITORS, ncol=K_FULL)  
   inspect_psm       <- inspect_delta      <- NULL
   true_state        <- rep(0, K_FULL)
   true_states_everyone <-qs_HULB_everyone <- qs_HMM_everyone <- 
     qs_HMM_nature_everyone   <- matrix(0, nrow=TOT_VISITORS, ncol=K_FULL) 
   error_state       <- matrix(0, nrow=TOT_VISITORS, ncol= 9); colnames(error_state)<- c("visitor", "qr_t1","qr_t2","qr_t3","qr_t4","qs_t1","qs_t2","qs_t3","qs_t4")
   
   alphabeta <- array(rep(1, TOT_STATES*TOT_MORPHS*2), dim=c(TOT_STATES, TOT_MORPHS, 2))
   alpha     <- alphabeta[,,1] ; beta   <- alphabeta[,,2]
   G_current_sm   <- matrix(0,  nrow=TOT_STATES, ncol=TOT_MORPHS) 
   
   #if(MODEL==MAB_ONLY) G_current_sm   <- matrix(0,    ncol=TOT_MORPHS)
   G_asym_m1      <- G_asym_m2    <-  matrix(0, nrow=TOT_VISITORS, ncol=TOT_STATES)  # current value of the  index
   monitor_G_m1_s <- monitor_G_m2_s <- matrix(0, nrow=TOT_VISITORS, ncol=TOT_STATES)
   
  ptime <- system.time({
   for (visitor  in 1:TOT_VISITORS)
   {
     # 3.1 Initialize clickstream, data structures to store state estimates (qs_system is HULB, qs_HMM is HMM) 
     Y                 <- rep(0, K_FULL);  Y_vec <-  NULL 
     qs_HMM_nature     <- qs_HMM    <- qs_HULB   <- matrix(, nrow=K_FULL+1, ncol=TOT_STATES) # click-by-click state estimates 
     qs_HMM_nature[1,]           <- QS_ARRIVAL_NATURE
     qs_HMM[1,]  <- QS_ARRIVAL  
     qs_HULB[1,] <- QS_ARRIVAL    
     
     # 3.2  Assign the initial morph - random in all scenarios
     ## Under RANDOM the morph does not change
     best_morph               <- which(rmultinom(1,1,matrix(c(1/TOT_MORPHS),ncol=TOT_MORPHS) ) == 1)  
     morph_chosen[visitor,1]  <- best_morph 
     
     # 3.3 Load current G (use asymptotical value if alpha or beta are greater than 3000) 
     for (st_ in 1: TOT_STATES) {  for (mor in 1:TOT_MORPHS)  { if(any(alpha[st_,mor] >= 3000, beta[st_,mor] >= 3000)) {
       G_current_sm[st_, mor] <- alpha[st_,mor]/(alpha[st_,mor]+beta[st_,mor])  } else {   
         G_current_sm[st_, mor] <- G_interpolation_fast(alpha[st_,mor],beta[st_,mor])  } } }
     
     # 3.4 Loop over all clicsks
     # click =1 # for testing purposes only  
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
         #  Lambda: transition matrices (one per morph so this is current state x  state x morph)
         
         ## C. Get the covariates for the transition matrix: XAs are specific to the links a consumer clicked on  
         XA      <- Ckjn[Ckjn[,1]==click, 2:ncol(Ckjn)] ;  XA <- XA[Y[click],];   XA[1] <- log(XA[1])
         
         ## D. computes transition probabilities for previous state 1 to state K-1/ in our case we only have p11 and p12. Because if one starts in previous state 2, she stays there. So p21=0 and p22=1, this makes Q0_K=c(0,1)
         lambda  <- Transition_Matrix(XA, rho, mu, DIM_STATES, TERMINAL) 
         
         ## E. stacks the transition matrix from state 1 to K-1, with the last row of the transition matrix, which is 0/1. Lambda is computed after every click
         if (TERMINAL) {lambda   <- array(c(rbind(lambda[,,1], Q0_K[,,1]), rbind(lambda[,,2], Q0_K[,,2])), dim=c(TOT_STATES, TOT_STATES, TOT_MORPHS))} 
         
         # F. Selects the transition matrix corresponding to the morph served and stores as Lambda_m
         lambda_m                <- lambda[,,best_morph]     
         
         # G. Nature knows the true state CHANGE!!!
         qs_HMM_nature[click+1,] <- lambda_m[true_state[click], ]  
         
         # H. Computes HMM estimate              
         #if (HMM_BANDIT) {   qs_HMM[click+1,]  <- qs_HMM[click,] %*% lambda_m }  # system does not know the true state so it uses lambda_m
         if (is.element(MODEL, c(HMM_BANDIT, HMM_MAB_NO_DP) )) {qs_HMM[click+1,] <- Bayesian_Thm_updater_element(lambda, qs_HMM[click,],best_morph, psm_true) }   
         
         # I. UPDATE HULB STATE: qs_HULB= f(click (t+1), m(t)  ) 
         if (MODEL==HULB) {qs_HULB[click+1, ]     <- Bayesian_updater(click, best_morph, QS_ARRIVAL, new_omega, Ckjn, all_s, Y_mat )[click,] }
         
         # J. UPDATE FUTURE STATE: qs_future = f(click (t+1), m(t) , lambda )
         if (MODEL == HMM_BANDIT) { qs_HMM_sys_future <- estimate_future_period_qs_general(lambda)} #, qs_HMM[click+1,] ) } # Gui replaced lambda instead of q0
         
         # K. Computes EGI: m*(t=1) = Argmax EGI(qs_hulb(t+1)  )
         if (MODEL == HMM_MAB_NO_DP)    {qs <- qs_HMM[click,  ] } 
         if (MODEL == HULB)             {qs <- qs_HULB[click, ] } 

         if ( is.element(MODEL, c(HULB, HMM_MAB_NO_DP) )    ) 
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
         ## TESTING  - changed qs_HMM[click+1,] to qs_HMM[click,]
         ## In the DP, changed exponent in (1-pbounce)^(next_period-K) into (period-K)
         if (MODEL == HMM_BANDIT) { best_morph_uncond <- DP(click, qs_HMM_sys_future, qs_HMM[click+1,], G_current_sm, p_bounce_click) }  # , delta_dp= delta_dp[period+1]) }
         
         # M. Assigns best morph according to the method of choice ( breaking ties randomly).
         if (MODEL == HMM_BANDIT)       { best_morph <- best_morph_uncond[click] }  # note that index should be click, not TOT_PERIODS
         if (MODEL == HMM_MAB_NO_DP)    { best_morph <- which.is.max(I[]) }     
         if (MODEL == HULB)             { best_morph <- which.is.max(I[]) }   

         morph_chosen[visitor, click+1] <-best_morph
         
         # N. Bounce away or stay: bounced[visitor] will have the number of the last click 
         if (click <K_FULL)
         { temp_prob   <- p_bounce_click[true_state[click], morph_chosen[visitor, click]] # p_bounce_mt[click, true_state[click], morph_chosen[visitor, click]]
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
     qs_HULB_everyone[visitor,  ]       <- qs_HULB[LAST_CLICK,1]
     qs_HMM_everyone[visitor,  ]        <- qs_HMM[LAST_CLICK,1]
     qs_HMM_nature_everyone[visitor,  ] <- qs_HMM_nature[LAST_CLICK,1]
     
     # 3.7 Draws purchase: delta_sm(t+1) = f( Psm(t), m(t), S(t) ) 
     psm_i <- psm_true[morph_chosen[visitor, (LAST_CLICK)], true_state[LAST_CLICK]] #true probability of purchase at every click using true state which changed at every click via HMM (unobserved by system)
     delta <-  rbinom(1, 1, psm_i)
     delta_sm[true_state[LAST_CLICK], morph_chosen[visitor, LAST_CLICK] ]<- delta_sm[true_state[LAST_CLICK], morph_chosen[visitor,LAST_CLICK] ] + delta # changed to LAST_CLICK
     
     if (is.element(MODEL, c(HMM_BANDIT, HMM_MAB_NO_DP)))  {qs_algo <- qs_HMM[LAST_CLICK,] }
     if (is.element(MODEL, c(RANDOM, HULB))) {qs_algo <- qs_HULB[LAST_CLICK,] } 
     # for random it does not matter, qs_algo is there so that the loop can run
     
     # 3.8 Updates alpha, beta  = f( qs_hulb(t+1) or f(qs_HMM(t+1) ) depending  on method being HULB or HMM
     for (s in 1:TOT_STATES) {
       alpha[s,best_morph] <- alpha[s,best_morph] +      delta   * qs_algo[s] 
       beta[s, best_morph]  <- beta[s,best_morph] + (1 - delta) * qs_algo[s]       
     }
     # Reset bouncing flag 
     BOUNCED <- FALSE
     
   } # close visitor loop 
   
 })
 ptime
   
ptime*1000














 #  ## scp ubuntu@145.100.59.55:/datastore/simhmm4_3NT.RData ~/Documents/server_files
 #  ## for the replicates
 #  ## scp -r ferecatu@lisa.surfsara.nl:rsm_sim/bin ~/Documents/server_files/out_files
 #  # setwd("/Users/alinaferecatu/Dropbox/LODE/morphing@RSM/analysis/simulation using synthetic data/Alina_sims_updateNov2019/results 1000 replicates/2S full with sep. delta bounce states")
 #  #  load("simhmm4_2Tfull_transitions.RData")
 #  #  load("simhulb4_2Tfull_transitions.RData")
 #   
 #   #  trials=1000
 #   #  simHMM=data.frame(resfull=simhmm4_2Tfull[,8], Cond=rep("HMM_4K", trials))
 #   #  simHULB=data.frame(resfull=simhulb4_2Tfull[,8], Cond=rep("HULB_4K", trials))
 #   #  simF=rbind(simHMM, simHULB)
 #   #  simF %>% group_by(Cond) %>% summarise(meanC=mean(resfull))
 #   #  summary(aov(resfull~Cond, data=simF))
 #   #  library(dplyr)
 #   #  simF=simF %>% mutate(Model=c(rep("Forward looking", trials), rep("Myopic", trials), rep("Random", trials)))
 #   #  dim(simF)
 #   #  head(simF)
 #   #  library("Hmisc")
 #   #  violinplot=ggplot(simF, aes(x = Cond, y =resfull, group=Cond))+
 #   #    geom_boxplot(adjust = 0.5, draw_quantiles = c(0.5)) +
 #   #    stat_summary(aes(group=Cond), fun.y=mean, geom="point",
 #   #                 fill="red", shape=21, size=3, position = position_dodge(width = .9))+
 #   #       stat_summary(aes(group=Cond), fun.data=mean_cl_normal, geom = "errorbar", mult = 1, fill="red")+
 #   #  facet_wrap(~ Model, scale="free")+ theme_grey()+
 #   #  scale_x_continuous(limits = c(0.03,0.06)) +
 #   #  scale_fill_manual(values = c("#00BFC4","#F8766D")) +
 #   #  labs(x = "Scenario", y = "Expected Success")
 #   #  violinplot
 #  
 #  ### bouncing hist
 #  bouncedHMM=simhmm4_2Tfull[,17:20]
 #  colMeans(bouncedHMM)
 #  bouncedHULB=simhulb4_2Tfull[,9:12]
 #  colMeans(bouncedHULB)
 # 
 #  
 #  simHMM=data.frame(resfull=bouncedHMM, Cond=rep("HMM_4K", trials))
 #  simHULB=data.frame(resfull=bouncedHULB, Cond=rep("HULB_4K", trials))
 #  
 #  simF=rbind(simHMM, simHULB)
 #  simF %>% group_by(Cond) %>% summarise(meanC=mean(resfull))
 #  summary(aov(resfull~Cond, data=simF))
 #  library(dplyr)
 #  # simF=simF %>% mutate(Model=c(rep("Forward looking", trials), rep("Myopic", trials), rep("Random", trials)))
 #  dim(simF)
 #  head(simF)
 #  
 #  library("Hmisc")
 #  bounce_plot=ggplot(simF, aes(x=resfull,fill=Cond))+
 #    geom_histogram(position="dodge2") 
 #  bounce_plot
 #  
 #  stat_summary(aes(group=Cond), fun.y=mean, geom="point",
 #               fill="red", shape=21, size=3, position = position_dodge(width = .9))+
 #  stat_summary(aes(group=Cond), fun.data=mean_cl_normal, geom = "errorbar", mult = 1, fill="red")+
 #    #facet_wrap(~ Model, scale="free")+ theme_grey()+
 #    #scale_x_continuous(limits = c(0.03,0.06)) +
 #    #scale_fill_manual(values = c("#00BFC4","#F8766D")) +
 #   labs(x = "Scenario", y = "Expected Success")
 #   violinplot
 #  
 #  
 #  
 #  #  successes per state and morph?
 #  simF %>% group_by(Cond) %>% summarize(a=mean(resfull))
 #  delta_sm
 #  (overall_success_rate=sum(delta_sm)/TOT_VISITORS)
 #   
 #  
 #  ## state transitions
 #  stateHMM=array(c(simhmm4_2Tfull[,9:16]), dim=c(1000,2, 4))
 #  stateHMM[1,,]
 #  stateClick1=stateHMM[,,1]  
 #  (propSC1=mean(stateClick1[,2]/rowSums(stateClick1)))
 #  stateClick2=stateHMM[,,2]  
 #  (propSC2=mean(stateClick2[,2]/rowSums(stateClick2)))
 #  stateClick3=stateHMM[,,3]  
 #  (propSC3=mean(stateClick3[,2]/rowSums(stateClick3)))
 #  stateClick4=stateHMM[,,4]  
 #  (propSC4=mean(stateClick4[,2]/rowSums(stateClick4)))
 #  c(mean(stateClick1[,1]), mean(stateClick2[,1]), mean(stateClick3[,1]), mean(stateClick4[,1]))
 #  c(mean(stateClick1[,2]), mean(stateClick2[,2]), mean(stateClick3[,2]), mean(stateClick4[,2]))
 #  
 #  ## state transitions
 #  stateHULB=array(c(simhulb4_2Tfull[,13:20]), dim=c(1000,2, 4))
 #  stateHULB[1,,]
 #  stateClick1=stateHULB[,,1]  
 #  (propSC1=mean(stateClick1[,2]/rowSums(stateClick1)))
 #  stateClick2=stateHULB[,,2]  
 #  (propSC2=mean(stateClick2[,2]/rowSums(stateClick2)))
 #  stateClick3=stateHULB[,,3]  
 #  (propSC3=mean(stateClick3[,2]/rowSums(stateClick3)))
 #  stateClick4=stateHULB[,,4]  
 #  (propSC4=mean(stateClick4[,2]/rowSums(stateClick4)))
 #  c(mean(stateClick1[,1]), mean(stateClick2[,1]), mean(stateClick3[,1]), mean(stateClick4[,1]))
 # c(mean(stateClick1[,2]), mean(stateClick2[,2]), mean(stateClick3[,2]), mean(stateClick4[,2]))
 # 
 # ## state-morph combo at the last click
 # ## average number of m2 served at the last click
 # statemHULB=rowSums(simhulb4_2Tfull[,23:24])
 # hist(statemHULB)
 # mean(statemHULB)
 # 
 # 
 # statemHMM=rowSums(simhmm4_2Tfull[,31:32])
 # hist(statemHMM)
 # mean(statemHMM)
 # 

   
   