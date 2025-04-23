# Function to run TMB Hg Pup Production model V1b on simulated data
# Last updated 2022-06-14 by EKJ
# Note that many things (e.g., upper and lower bounds for some parameters and the
# number of days in the season) are hardcoded in this version. 

runTMB <- function(data, problist, model, true){
  
  # PURPOSE
  # Run V1b of the TMB Hg Pup Production model given simulated data
  # INPUTS
  # - data:   A dataframe containing all observations for a given colony and year
  # - problist: list of Pm and Pl
  # - model:  a string identifying the model filename (without the .cpp extension)
  #           the model file is assumed to be in the ./Models/ directory 
  # - version: a string, either "V1" or "V2" "V1b"
  # OUTPUTS
  # - rep:      the output of sdreport(fun)
 
  
  library(TMB)
  library(optimx)

  # Compile and load TMB model
  compile(paste0("./Models/", model, ".cpp"), "-O1 -g",DLLFLAGS="")
  dyn.load(dynlib(paste0("./Models/", model)))
  
  # Create a list of data and fixed parameters to be passed to TMB
  TMBdata <- list()

  # Number of individuals counted
  TMBdata$nwhite <- data$WhiteObs
  TMBdata$nmoult <- data$MoultObs

  # Days on which counts happened (relative to Oct 1st)
  TMBdata$obsdays <- data$Day - 1

  # Probability of moulting and leaving at each age
  # have to add zeroes because season = 123 days, problist = 100 days
  TMBdata$pmoult <- c(problist$Pm, rep(0, 2))
  TMBdata$pleave <- c(problist$Pl, rep(0, 2))
  TMBdata$psmoult <- c(problist$PS4, rep(0, 2))
  
  # Probability of detecting white and moulted pups 
  TMBdata$powhite <- 0.95
  TMBdata$pomoult <- 0.95

  # Probability of correctly classifying white and moulted pups
  TMBdata$pcwhite <- 1
  TMBdata$pcmoult <- 0.91
    
  # Number of days in the season
  TMBdata$ndays <- 125

  # Vector of 1-indexed birthdays 
  TMBdata$bday <- 1:125

  
  # these are the intervals at which the birth curve will be evaluated
  # to create the cumulative pborn
  TMBdata$slices <- seq(0, 125, length.out = 12500)
  TMBdata$start <- seq(0, 100*125, by = 100)

  # Set lower and upper bounds for parameters

  lower <- list("N" = max(data$WhiteObs+data$MoultObs), 
              "mubday" = 30, # mubday in days after oct 1st
              "sdbday" = 1, # 
              "alphabday" = 1) # skew parameter

  upper <- list("N" = max(data$WhiteObs+data$MoultObs)*5,
              "mubday" = 70, # mubday in days after oct 1st
              "sdbday" = 20,
              "alphabday" = 10)# skew parameter 

  default <- list("N" = mean(c(lower$N, upper$N)),
                "mubday" = mean(c(lower$mubday, upper$mubday)),
                "sdbday" = mean(c(lower$sdbday, upper$sdbday)),
                "alphabday" = mean(c(lower$alphabday, upper$alphabday)))

  # Initialize model at default parameter values

  fun <- TMB::MakeADFun(data=TMBdata, 
                      parameters = default,
                      DLL = model, 
                      checkParameterOrder = FALSE, silent=TRUE)

  # OPTIMIZATION

  # use grid search to choose 5 starting values
  
    start.df <- expand.grid("N" = seq(lower[[1]], upper[[1]], length.out = 12)[2:10],
                            "mubday" = seq(lower[[2]], upper[[2]], length.out = 12)[2:10],
                            "sdbday" = seq(lower[[3]], upper[[3]], length.out = 12)[2:10],
                            "alphabday" = seq(lower[[4]], upper[[4]], length.out = 12)[2:10])
    
    
    start.df$value <- NA
    for (i in 1:nrow(start.df)){
      
      start.df$value[i] <- fun$fn(start.df[i,1:4])
      
    }
    
    sort.start <- start.df[order(start.df$value),]
    sort.start <- sort.start[!is.na(sort.start$value),]
    best.gstart <- sort.start[1:5,]
    
    # nse random draws to choose 5 values
    
    nrstart <- 100
    rstart.df <- data.frame("N" = runif(nrstart, lower[[1]], upper[[1]]),
                            "mubday" = runif(nrstart, lower[[2]], upper[[2]]),
                            "sdbday" = runif(nrstart, lower[[3]], upper[[3]]),
                            "alphabday" = runif(nrstart, lower[[4]], upper[[4]]))
    
    rstart.df$value <- NA
    for (i in 1:nrow(rstart.df)){
      
      rstart.df$value[i] <- fun$fn(rstart.df[i,1:4])
      
    }
    
    sort.rstart <- rstart.df[order(rstart.df$value),]
    sort.rstart <- sort.rstart[!is.na(sort.rstart$value),]
    
    # stick together best gridded and random starting values
    best.start <- rbind.data.frame(best.gstart, sort.rstart[1:5,])

  best.start <- na.omit(na.omit(best.start))
    
  # Now optimize with BFGS 

  optim.df <- data.frame("Conv" = rep(NA, nrow(best.start)), 
                       "Value" = rep(NA, nrow(best.start)))
  
  optim.list <- list()

  for (i in 1:nrow(best.start)){
    toofun <- function(...){
      val <- fun$fn(...)
      #cat("lnl=", val, "\n")
      #cat("par=", paste(..., collapse=", "), "\n")
      val
    }
    optim.out <- optimx(par = as.numeric(best.start[i,1:4]), fn = toofun, #gr = fun$gr, 
                        method = "L-BFGS-B", lower = as.numeric(lower), upper = as.numeric(upper), 
                        control = list(maxit = 100000))
    
    optim.df[i,] <- c(optim.out$convcode, optim.out$value)
    optim.list[[i]] <- optim.out
  } # end for i in best start

  
  # Choose the best set of parameters and update the model
  # check which of best.start converged (and wasn't NULL)
  optim.nonnans <- which(!is.na(optim.df$Conv) & optim.df$Conv == 0)
  # filter the optim.df to include only valid parameter sets
  optim.df.nonnans <- optim.df[optim.nonnans, ]
  # filter the optim.list to include only valid parameter sets
  optim.list.nonans <- optim.list[optim.nonnans]
  # find the best of the parameter sets
  optim.best <- which.min(optim.df.nonnans[optim.df.nonnans$Conv==0, 2])
  
  # if none of the par sets meet the criteria, choose the first set
  #ifelse(length(optim.best) == 0, optim.out <- optim.list[[1]], 
  #   optim.out <- optim.list.nonans[[optim.best]])
  
  if(length(optim.best)>0){
    optim.out <- optim.list.nonans[[optim.best]]
  }else{
    return(list(rep = list(par=rep(NA, 4)), conv = 9999, 
                estNLL = NA, trueNLL = NA))
  }
  
  # fiddle with output so MakeADFun is happy
  optimpar <- as.list(optim.out[,1:4])
  names(optimpar) <- names(upper)
  
  # there is probably a more elegant way to do this
  fun <- TMB::MakeADFun(data=TMBdata, 
                        parameters = optimpar,
                        DLL = model, 
                        checkParameterOrder = FALSE)

  # Compile results

  rep <- sdreport(fun)
  
  return(list(rep = rep, conv = optim.out$convcode, 
              estNLL = fun$fn(optimpar), trueNLL = fun$fn(true)))
  
} # end runTMB  



