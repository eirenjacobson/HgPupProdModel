# Script to run TMB model using real data
# Last updated 2024-08-02 by EKJ

##### USER INPUTS
datafile <- "./Data/Hg_pup_counts_2018_temp.csv"
seed <- 513807
model <- "HgPupProd_TMBV1b"

# Function to run TMB model 
source("./Scripts/runTMB_RealData_V1b.R")

# Data simulation script
# Used only to generate vectors of Pm and Pl

  source("./Scripts/HgSim_V1b.R")
  sim.out <- simSeals()
  problist <- list("Pm" = sim.out$Pm, "Pl" = sim.out$Pl)


#####

library(tidyr)
library(dplyr)
library(ggplot2)

set.seed(seed)

data <- read.csv(datafile)

colyears <- data %>% select(Colony, Year) %>% distinct()

results <- rbind.data.frame(cbind(colyears, "Model" = "TMB_V1b"))

results$N.Est <- NA; results$N.CV <- NA; results$N.SE <- NA; 
results$N.LCI <- NA; results$N.UCI <- NA
results$mubday.Est <- NA; results$sdbday.Est <- NA; results$alphabday.Est <- NA

realdataout <- list()

for (i in 1:nrow(results)){
    
    realdata <- subset(data, Colony == results$Colony[i] & Year ==results$Year[i])
    out <- runTMB(data = realdata, problist = problist, model = model)
    realdataout[[i]] <- out
    results$N.Est[i] <- out$rep$par.fixed[1]
    results$N.SE[i] <- summary(out$rep)[1,2]
    results$N.sample[i] <- nrow(realdata)
    results$N.LCI[i] <- results$N.Est[i] - (1.96*results$N.SE[i])
    results$N.UCI[i] <- results$N.Est[i] + (1.96*results$N.SE[i])
    results$mubday.Est[i] <- out$rep$par.fixed[2]
    results$sdbday.Est[i] <- out$rep$par.fixed[3]
    results$alphabday.Est[i] <- out$rep$par.fixed[4]
  }
  
for (i in 1:nrow(results)){
  
  realdata <- subset(data, Colony == results$Colony[i] & Year ==results$Year[i])
  results$samplesize[i] <- nrow(realdata)
}

results$N.CV <- results$N.SE/results$N.Est
results$C <- exp(1.96 * sqrt(log(1+results$N.CV^2)))
results$N.LCI <- results$N.Est/results$C
results$N.UCI <- results$N.Est*results$C

save(results, file = paste0("./Data/Summary_RealData_", model, ".RData"))
save(realdataout, file = paste0("./Data/Summary_RealData_", model, "_rep", ".RData"))

hibyresults <- read.csv("./Data/Hg_pup_prod_output_SET_DEFAULT_DIG2020-11-11.csv")

ggplot(results) +
  geom_point(aes(x=Colony, y = N.Est)) +
  geom_errorbar(aes(x=Colony, ymin = N.LCI, ymax = N.UCI), width = .25) +
  geom_point(data = hibyresults, aes(x = Colony, y = PROD), col = "red") +
  theme_bw()

#ggsave(plot = last_plot(), filename = paste0("./Figures/RealData_", model, "wHiby.png"), 
#       width = 6, height = 4, units = "in")
  