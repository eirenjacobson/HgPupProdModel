# Script to run TMB models using simulated data
# Last updated 2024-02-14 by EKJ

library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(lubridate)

##### USER INPUTS
gendata = TRUE
seed <- 20240214
model <- "HgPupProd_TMBV1b"
datasets <- "SimulatedData"
version <- "V1b"
simpopproc <- "detcont"
simobsproc <- "normal"
ndays <- 125

# Function to run TMB model 
source("./Scripts/runTMB_SimData_V1b.R")

# Function to simulate data
source("./Scripts/HgSim_V1b.R")

# IF gendata = TRUE, these are the parameters that will be used for simulation

# Parameters that will be estimated by the models
N.vals <- seq(100, 1000, by = 100)
#N.vals <- seq(100, 1000, by = 100)
#N.vals <- 100000
skewbday.vals <- 5
mubday.vals <- 50
sdbday.vals <- 10

# obs pars
nobsday.vals <- c(1:5)
start.vals <- c(25, 30, 35, 40, 45, 50)
#interval.vals <- 10 # for 5 obsdays
#interval.vals <- 5 # for 10 obsdays
interval.vals <- c(1, 5, 10, 15, 20) # 1 for 90 obsdays
set.seed(seed)

sim.df <- expand.grid(N.vals, 
                        mubday.vals, sdbday.vals, skewbday.vals, 
                        nobsday.vals, start.vals, interval.vals)
  
names(sim.df) <- c("N.set", "mubday.set", "sdbday.set", "skewbday.set", "nobsday.set", "obs.start", "obs.int")
  
sim.list <- list()

sim.df$simpopproc <- simpopproc
sim.df$simobsproc <- simobsproc
sim.df$TMB.conv <- NA
sim.df$TMB.NEst <- NA
sim.df$TMB.NSE <- NA
sim.df$TMB.mubday <- NA
sim.df$TMB.sdbday <- NA
sim.df$TMB.skewbday <- NA
sim.df$TMB.estNLL <- NA
sim.df$TMB.trueNLL <- NA
sim.df$N.real <- NA

sim.df$TMBModel <- model

# Data simulation script
# Used only to generate vectors of Pm and Pl

  sim.out <- simSeals()
  problist <- list("Pm" = sim.out$Pm, "Pl" = sim.out$Pl)
  
  for (s in 1:nrow(sim.df)){
    
    #set.seed(s)
    print(paste("Beginning iteration s =", s, "of", nrow(sim.df), sep = " "))
    
    if(gendata == FALSE){sim.out <- sim.list[[s]]}
    
    if(gendata == TRUE){
      
      obsdays <- seq(sim.df$obs.start[s], by = sim.df$obs.int[s], length.out = 100)
      obsdays <- as.numeric(na.omit(obsdays[which(obsdays<ndays)][1:sim.df$nobsday.set[s]]))
      
      #Run simulation
      
      sim.out <- simSeals(N = sim.df$N.set[s],
                          skew.bday = sim.df$skewbday.set[s],
                          mu.bday = sim.df$mubday.set[s],
                          sd.bday =  sim.df$sdbday.set[s],
                          popproc = simpopproc, obsproc = simobsproc, 
                          obsdays = obsdays, ndays = ndays)
      
      # store simulated datasets for future use
      sim.list[[s]] <- sim.out
      
    }
    
    true <- data.frame("N" = sim.df$N.set[s],
                       "mubday" = sim.df$mubday.set[s],
                       "sdbday" = sim.df$sdbday.set[s], 
                       "alphabday" = sim.df$skewbday.set[s])
                       
    
    tmb.out <- runTMB(data = sim.out$greydata,
                      model = model,
                      problist = problist, 
                      true = true)
    
    # Compile results
    sim.df$N.real[s] <- sum(sim.out$realBorn)
    
    sim.df$TMB.conv[s] <- tmb.out$conv
    sim.df$TMB.NEst[s] <- tmb.out$rep$par[1]
    sim.df$TMB.NSE[s] <-   summary(tmb.out$rep)[1,2]
    sim.df$TMB.mubday[s] <- tmb.out$rep$par[2]
    sim.df$TMB.sdbday[s] <- tmb.out$rep$par[3]
    sim.df$TMB.skewbday[s] <- tmb.out$rep$par[4]
    sim.df$TMB.estNLL[s] <- tmb.out$estNLL
    sim.df$TMB.trueNLL[s] <- tmb.out$trueNLL
    
    print(paste("Completed iteration s =", s, "of", nrow(sim.df), sep = " "))
    
  } # end for s
  
  sim.df$nobsday <- NA
  
  for (i in 1:nrow(sim.df)){
    
    sim.df$nobsday[i] <- nrow(sim.list[[i]]$greydata)
  }
  
  sim.df$N.CV <- as.numeric(sim.df$TMB.NSE)/sim.df$TMB.NEst
  sim.df$C <- exp(1.96 * sqrt(log(1+sim.df$N.CV^2)))
  sim.df$N.LCI <- sim.df$TMB.NEst/sim.df$C
  sim.df$N.UCI <- sim.df$TMB.NEst*sim.df$C
  
  sim.df$N.PErr <- (sim.df$TMB.NEst-sim.df$N.set)/sim.df$N.set
  sim.df$mubday.PErr <- (sim.df$TMB.mubday-sim.df$mubday.set)/sim.df$mubday.set
  sim.df$sdbday.PErr <- (sim.df$TMB.sdbday-sim.df$sdbday.set)/sim.df$sdbday.set
  sim.df$skewbday.PErr <- (sim.df$TMB.skewbday-sim.df$skewbday.set)/sim.df$skewbday.set
  
  save(sim.df, file = paste0("./Data/", "SimResults_", simpopproc, simobsproc,"_", model, ".RData"))

  
  

  
  