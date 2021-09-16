# Configuration file
# Run this file before running code from Study1 or Study2 folders

# clean up R environment
rm(list=ls()) 
set.seed(9000)

# Context:  "PHONE STORE" or "HULB2009"
EMPIRICAL_SETTING     <- "PHONE STORE"  

# Model specs: "HULB","HMM_MAB_NO_DP","MAB_NO_HMM","HMM_BANDIT",  "HMM_BANDIT_STORE", "RANDOM" 
MODEL       <- "HMM_BANDIT_STORE" #"HMM_MAB_NO_DP" #"HMM_BANDIT" # "RANDOM" #    

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
QS_ARRIVAL        <- if(TOT_STATES==2) {c(0.55, 0.45)} else {c(0.34, 0.33, 0.33) } # c{1/3,1/3,1/3) }  or {c(6/11, 6/22, 6/33)}  
QS_ARRIVAL_NATURE <- if(TOT_STATES==2) {c(0.98, 0.02)} else {c(.98,.015, 0.005)}     # Alina suggestion: the HMM estimation is throwing everyone at state 1 at t0. Gui will think more about it as it can have huge pushback
QS_ARRIVAL        <- if(TOT_STATES==2) {c(1, 0)} else {c(1, 0, 0) } # Alina suggestion: everyone starts in state 1 in the RCT
QS_ARRIVAL_NATURE <- if(TOT_STATES==2) {c(1, 0)} else {c(1, 0, 0)}  

# Global parameters, 
TOT_VISITORS <- 100000 
TOT_MORPHS   <- 2  
K_FULL       <- 14  #  4 in MBA/first two rounds: 4  periods so we allow to change morphs every 4 sets of clicks. Alina suggests (a) using 15 (clicks, RCT median=15, mean=20). (b) decouple clicks from periods. This should be thought thorough
TOT_PERIODS=K_FULL
TOT_CONSIDERED_PERIODS<-4


# paths to use
PATH          <-  "~/Github/Replication_Morphing_HMM/" 
PATH_RAW       <-  paste("1. Raw Data",sep="")
PATH_HMM_EST   <-  paste("2. HMM Estimates", sep="")
PATH_OUT      <-  paste(PATH,"_results",sep="")
DROPBOX_GMATRIX_LINK= 'https://dl.dropboxusercontent.com/s/5gv7jquou6y3tlw/Gmatrix.out'

# set the working directory to path
setwd(PATH)



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


# set the memory limit
memory.limit(size=50000)

# load functions
source('Functions.R') 

# Load libraries 
# Packes required for subsequent functions. P_load ensures these will be installed and loaded. 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               stringr,
               stargazer,
               qdapRegex, 
               anytime,
               stringr,
               writexl,
               reshape2,
               maotai,
               nnet,
               expm,
               tidyverse,
               doParallel,
               rstan,
               loo,
               devtools) 

# special case: install from source
# if already installed: just load with library()
install.packages("RcppParallel", type = "source")
library(RcppParallel)
