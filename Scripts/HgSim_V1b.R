############ Simulation of grey seal pup production V1b ###########################
# Last modified 2023-05-04 by EKJ

simSeals <- function(N = 500, 
   ndays = 123, amax = 123, 
   birthcurve = "skewnormal",
   skew.bday = 5, 
   start = 32, mu.bday = 45, sd.bday = 10, 
   mu.mday = 23, sd.mday = 5,
   mu.lday = 31.5, sd.lday = 7,
   PObs = 0.95,
   PCWhite = 1, PCMoult = 0.91, PCMoultSD = 0.000000000001, 
   obsdays = c(35, 40, 50, 65, 70, 80, 95),
   popproc = "detcont", # option for "stocdisc"
   obsproc = "exact", # options for "normal" and "binom"
   seed = sample(1:1E6, 1)){

   library(lubridate)

  #########
  
  if(birthcurve == "lognormal"){
  
  # convert birth parameters for lognormal dist
  lCV <- sd.bday/(mu.bday-start)
  lsd <- sqrt(log(lCV^2+1))
  lmu <- log(mu.bday-start)-(0.5*(lsd^2))
  
  # calculate probability of being born on each day
  Pb <- rep(0, ndays)
  
  for(t in 2:ndays){
    Pb[t] <- plnorm((t-start), lmu, lsd) - plnorm((t-start-1), lmu, lsd)
  }
  
  Pb <- Pb/sum(Pb)}
  
  if(birthcurve == "skewnormal"){
    
    library(brms)
  
    Pb <- rep(0, ndays)
    
    for(t in 2:ndays){
      Pb[t] <- pskew_normal(q = t, mu = mu.bday, sigma = sd.bday, alpha = skew.bday) - 
        pskew_normal(q = t-1, mu = mu.bday, sigma = sd.bday, alpha = skew.bday)
    }
    
    Pb <- Pb/sum(Pb)
    
  }
  
  # calculate the probability of moulting at each age
  
  Pm <- rep(0, amax)
  
  for (t in 1:amax){
    Pm[t] <- pnorm(t, mu.mday, sd.mday) - pnorm(t-1, mu.mday, sd.mday)
  }

  Pm <- Pm/sum(Pm)
  
  # calculate the probability of leaving at each age
  
  Pl <- rep(0, amax)
  
  for (t in 2:amax){
    Pl[t] <- pnorm(t, mu.lday, sd.lday) - pnorm(t-1, mu.lday, sd.lday)
  }
  
  Pl <- Pl/sum(Pl)
  
  # process model if deterministic and continuous
  if (popproc == "detcont"){
      
      realBorn <- Pb*N 
      
      whitePups <- rep(0, ndays)
      moultPups <- rep(0, ndays)
      
      for (d in 1:(ndays)){
        for (i in 1:d){
          
          age = d-i
          
          whitePups[d] <- sum(whitePups[d], realBorn[i] * (1-sum(Pm[1:age])) * (1-sum(Pl[1:age])))
          moultPups[d] <- sum(moultPups[d], realBorn[i] * sum(Pm[1:age]) * (1-sum(Pl[1:age])))
          
        } # end for d
      } # end for i
      
  } # end if detcont
  
  # process model if stochastic and discrete
  if (popproc == "stocdisc"){
    
      realBorn <- 0
      while(sum(realBorn) != N){
        realBorn <- rbinom(n = ndays, size = N, prob = Pb)
      }
    
      leavers <- matrix(rep (0, amax*amax), nrow = amax)
      realMoult <- matrix(rep(0, ndays*amax), nrow = ndays, ncol = amax)
      realLeave <- matrix(nrow = ndays, ncol = amax)
      
      for (i in 1:ndays){
        
        while(sum(realMoult[i,]) != sum(realBorn[i])){
          realMoult[i, ] <- rbinom(n = amax, size = realBorn[i], prob = Pm)
        }
        
        for (j in 1:amax){
          while(sum(leavers[j,]) != realMoult[i, j]){
            leavers[j, ] <- c(rep(0, j), rbinom(size = realMoult[i, j], prob = Pl, n = amax)[j:(amax-1)])
          } # end while
        } # end for j
        
        realLeave[i,] <- colSums(leavers)
        
      } # end for i
      
      # now count number of individuals per day
      
      # initialize vectors
      
      tmin <- min(which(realBorn!=0))
      
      moultvec <- rep(0, ndays) # number that have moulted by each day
      leavevec <- rep(0, ndays) # number that have left by each day
      
      for (i in tmin:ndays){

        mv <- 0
        lv <- 0
        
        rowmax <- i - tmin + 1
        
        for (j in tmin:i){
          
          for (k in 1:rowmax){
            
            mv <- sum(mv, realMoult[j, k])
            lv <- sum(lv, realLeave[j, k])
            
          } # end for k
          
          rowmax <- rowmax - 1
        }
        
        moultvec[i] <- mv
        leavevec[i] <- lv
        
      } # end for i
      
      whitePups <- cumsum(realBorn) - moultvec
      moultPups <- moultvec - leavevec
      allPups <- cumsum(realBorn) - leavevec
  
  } # end if stocdisc 
  
  if (obsproc == "exact"){
  
    whiteObs <- whitePups*PObs*PCWhite + moultPups*PObs*(1-PCMoult)
    moultObs <- moultPups*PObs*PCMoult + whitePups*PObs*(1-PCWhite) }
  
  if (obsproc == "binom"){
    
    whiteDetected <- rbinom(length(whitePups), whitePups, PObs)
    moultDetected <- rbinom(length(moultPups), moultPups, PObs)
    
    WW <- rbinom(length(whiteDetected), whiteDetected, PCWhite)
    WM <- whiteDetected - WW
    MM <- rbinom(length(moultDetected), moultDetected, PCMoult)
    MW <- moultDetected - MM
    whiteObs <- WW + MW
    moultObs <- MM + WM
    
  }
    
  if (obsproc == "normal"){

    whiteDetectSD <- sqrt(whitePups*PObs * (1-PObs))
    moultDetectSD <- sqrt(moultPups*PObs * (1-PObs))

    whiteDetect <- rnorm(length(whitePups),
                         mean = whitePups*PObs,
                         sd = whiteDetectSD)

    moultDetect <- rnorm(length(moultPups),
                         mean = moultPups*PObs,
                         sd = moultDetectSD)

    whiteClassSD <- sqrt(abs(whiteDetect)*PCWhite * (1-PCWhite))
    moultClassSD <- sqrt(abs(moultDetect)*PCMoult * (1-PCMoult))
    WW <- rnorm(length(whiteClassSD),
                mean = whiteDetect*PCWhite,
                sd = whiteClassSD)
    WM <- whiteDetect - WW
    MM <- rnorm(length(moultDetect),
                mean = moultDetect*PCMoult,
                sd = moultClassSD)
    MW <- moultDetect - MM

    whiteObs <- abs(WW + MW)
    moultObs <- abs(MM + WM)

  }
  
    # then subsample for observations on certain days
    
    wObs <- whiteObs[obsdays]
    mObs <- moultObs[obsdays]
    
    # create a data frame
    greydata <- data.frame("Day" = obsdays, "WhiteObs" = wObs, "MoultObs" = mObs)
    
    # create output needed for Hiby model
    
    sim.out <- list(realBorn = realBorn, Pb = Pb, Pm = Pm, Pl = Pl, 
                    whitePups = whitePups, moultPups = moultPups,
                    greydata = greydata, 
                    whiteObs = whiteObs, moultObs = moultObs)
    
    return(sim.out)

}



