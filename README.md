# HgPupProdModel

This repository contains code to accompany Jacobson et al. (In Prep): A state-space model for estimating pinniped pup production from serial counts at breeding colonies.  

# Required packages

brms  
dplyr  
ggplot2  
lubridate  
magrittr  
optimx  
tidyr  
TMB  

# Running the model

- The file Scripts/runColonies_SimData_V1b.R runs the model on data from multiple simulated colonies
- The file Scripts/runColonies_RealData_V1b.R runs the model on data from multiple real colonies

# Repository contents

- Data/Hg_pup_counts_2018_temp.csv contains real pup counts from six grey seal colonies
- Data/Hg_pup_prod_output_SET_DEFAULT_DIG2020-11-11.csv contains output from the Hiby model as described in Russell et al. 2019
- Model/HgPupProd_TMBV1b.cpp contains the TMB model (V1b indicates that this is version 1b of the model; this is for internal reference)
- Scripts/HgSim_V1b.R contains a function to simulate serial counts of grey seal pups
- Scripts/runColonies_RealData_V1b.R contains a script to run the model on data from multiple real colonies
- Scripts/runColonies_SimData_V1b.R contains a script to run the model on data from multiple simulated colonies
- Scripts/runTMB_RealData.R contains a function to run the TMB model on a single real colony
- Scripts/runTMB_SimData_V1b.R contains a function to run the TMB model on a single simulated colony
