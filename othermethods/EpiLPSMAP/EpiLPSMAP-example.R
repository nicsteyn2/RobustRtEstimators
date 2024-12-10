
# Make sure your working directory is set to the main directory of the repo
# install.packages("EpiLPS", version="1.3.0")
# install.packages("hypergeo")
# install.packages("scoringRules)

rm(list=ls())

library(tidyverse)
library(foreach)
library(EpiLPS)
library(Rcpp)

sourceCpp("othermethods/EpiLPSMAP/KercubicBspline.cpp") # We need this C++ function to calculate the cubic B-spline basis functions
source("othermethods/EpiLPSMAP/LPSMAP-Incidence.R")
source("othermethods/EpiLPSMAP/LPSMAP-CRPS.R")

# Load simulated data
df = read.csv("othermethods/exampledata.csv")
simulations = c("Random walk", "Sinusoidal", "Step change")

# Create discretised serial interval (matching our simulated data)
w <- dgamma(1:100, shape = 2.3669, scale = 2.7463)
w = w/sum(w)

# Explicitly define the prior distributions at their default hyperparameter values
priordists = Rmodelpriors(list(
  a_delta = 10,
  b_delta = 10,
  phi = 2,
  a_rho = 1e-04,
  b_rho = 1e-04
))

# Pre-allocate dataframe to store parameter estimates
df_params = data.frame(simulation=simulations, lambda=NA, rho=NA)

# Fit model to each simulation
df_results = foreach (ii = seq(1, length(simulations)), .combine=rbind) %do% {
  
  print(paste0("Fitting simulation: ", simulations[ii]))
  
  # Extract incidence data for current simulation
  Ct = df %>% filter(simulation==simulations[ii]) %>% arrange(t) %>% pull(Cases)
  
  # Fit EpiLPS
  LPSfit = estimR(incidence = Ct, si = w, priors=priordists, K = 100)
  df_R = LPSfit$RLPS %>% rename(t=Time, Mean=R, Lower=Rq0.025, Upper=Rq0.975) %>% dplyr::select(t, Mean, Lower, Upper) %>% mutate(variable="R")
  
  # Fetch cases
  LPScases = observedIncidenceCurve(LPSfit)
  df_I = LPScases %>% rename(t=Time, Mean=I) %>% mutate(variable="I")
  
  # Fetch CRPS
  df_R$CRPS = NA
  print(paste0("Computing CRPS..."))
  df_I$CRPS = computeCRPS(LPSfit)

  # Fetch parameter values
  if (!(LPSfit$optimconverged)) {
    print("Optimisation did not converge")
  } else {
    print(paste0("lambda (smoothness parameter): ", LPSfit$penparam))
    print(paste0("rho (reporting overdispersion): ", LPSfit$NegBinoverdisp))
    df_params$simulation[ii] = simulations[ii]
    df_params$lambda[ii] = LPSfit$penparam
    df_params$rho[ii] = LPSfit$NegBinoverdisp
  }
  
  # Return results
  return(rbind(df_R, df_I) %>% mutate(simulation=simulations[ii]))
  
}

# Append true values
df_truth = df %>% rename(R = TrueRt, I=Cases) %>% pivot_longer(cols=c("R", "I"), names_to="variable", values_to="TrueValue")
df_results = left_join(df_results, df_truth, by=c("t", "variable", "simulation"))

# Store parameter estimates in standard format
df_params_out = df_params %>% pivot_longer(cols=c("lambda", "rho"), names_to="parameter", values_to="mode")

# Save results
write.csv(df_results, "othermethods/EpiLPSMAP/results.csv", row.names=FALSE)
write.csv(df_params_out, "othermethods/EpiLPSMAP/params.csv", row.names=FALSE)
