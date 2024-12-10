
# Make sure your working directory is set to the main directory of the repo
# install.packages("EpiLPS", version="1.3.0")
# install.packages("hypergeo")
# install.packages("scoringRules)

rm(list=ls())

library(tidyverse)
library(foreach)
library(EpiLPS)
library(Rcpp)

sourceCpp("othermethods/EpiLPSMALA/packagefiles/KercubicBspline.cpp")
sourceCpp("othermethods/EpiLPSMALA/packagefiles/KerLaplace.cpp")
sourceCpp("othermethods/EpiLPSMALA/packagefiles/KerMVN.cpp")
sourceCpp("othermethods/EpiLPSMALA/packagefiles/KerRpostmcmc.cpp")

source("othermethods/EpiLPSMALA/packagefiles/estimRmcmc-custom.R")
source("othermethods/EpiLPSMALA/packagefiles/KerIncidCheck.R")
source("othermethods/EpiLPSMALA/packagefiles/KerPtheta.R")
source("othermethods/EpiLPSMALA/packagefiles/KerLikelihood.R")
source("othermethods/EpiLPSMALA/packagefiles/KerMCMC.R")

source("othermethods/EpiLPSMALA/LPSMALA-incidence.R")
source("othermethods/EpiLPSMALA/LPSMALA-CRPS.R")

# Load simulated data
df = read.csv("othermethods/exampledata.csv")
simulations = c("Random walk", "Sinusoidal", "Step change")

# Create discretised serial interval (matching our simulated data)
w <- dgamma(1:100, shape = 2.3669, scale = 2.7463)
w = w/sum(w)

# Explicitly define the prior distributions at their default hyper parameter values
priordists = Rmodelpriors(list(
  a_delta = 10,
  b_delta = 10,
  phi = 2,
  a_rho = 1e-04,
  b_rho = 1e-04
))

# Pre-allocate dataframe to store parameter estimates
df_params = data.frame()

# Fit model to each simulation
df_results = foreach (ii = seq(1, length(simulations)), .combine=rbind) %do% {

  print(paste0("Fitting simulation: ", simulations[ii]))
  
  # Extract incidence data for current simulation
  Ct = df %>% filter(simulation==simulations[ii]) %>% arrange(t) %>% pull(Cases)
  
  # Fit EpiLPS
  LPSfit = estimRmcmc_custom(incidence = Ct, si = w, priors=priordists, K = 100)
  df_R = LPSfit$RLPS %>% rename(t=Time, Mean=R, Lower=Rq0.025, Upper=Rq0.975) %>% dplyr::select(t, Mean, Lower, Upper) %>% mutate(variable="R")
  
  # Fetch cases
  LPScases = observedIncidenceCurve(LPSfit)
  df_I = LPScases %>% rename(t=Time, Mean=I) %>% mutate(variable="I")
  
  # Fetch CRPS
  df_R$CRPS = NA
  print(paste0("Computing CRPS from MCMC samples"))
  df_I$CRPS = computeCRPS(LPSfit)
  # df_I$CRPS = 0
  
  # Print parameter values
  kde = density(LPSfit$MCMC$rho_mcmc)
  moderho = kde$x[which.max(kde$y)]
  meanrho = mean(LPSfit$MCMC$rho_mcmc)
  lowerrho = quantile(LPSfit$MCMC$rho_mcmc, 0.025)
  upperrho = quantile(LPSfit$MCMC$rho_mcmc, 0.975)
  print(paste0("MAP rho: ", moderho, " (95% CI: ", lowerrho, ", ", upperrho, ")"))
  
  kdelambda = density(LPSfit$MCMC$lambda_mcmc)
  modelambda = kdelambda$x[which.max(kdelambda$y)]
  meanlambda = mean(LPSfit$MCMC$lambda_mcmc)
  lowerlambda = quantile(LPSfit$MCMC$lambda_mcmc, 0.025)
  upperlambda = quantile(LPSfit$MCMC$lambda_mcmc, 0.975)
  print(paste0("MAP lambda: ", modelambda, " (95% CI: ", lowerlambda, ", ", upperlambda, ")"))
  
  # Store parameter values
  df_params = rbind(df_params, data.frame(
    simulation = c(simulations[ii], simulations[ii]),
    parameter = c("rho", "lambda"),
    mean = c(meanrho, meanlambda),
    mode = c(moderho, modelambda),
    lower = c(lowerrho, lowerlambda),
    upper = c(upperrho, upperlambda)
  ))
  
  # Return results
  return(rbind(df_R, df_I) %>% mutate(simulation=simulations[ii]))
  
}

# Append true values
df_truth = df %>% rename(R = TrueRt, I=Cases) %>% pivot_longer(cols=c("R", "I"), names_to="variable", values_to="TrueValue")
df_results = left_join(df_results, df_truth, by=c("t", "variable", "simulation"))

# Save results
write.csv(df_results, "othermethods/EpiLPSMALA/results.csv", row.names=FALSE)
write.csv(df_params, "othermethods/EpiLPSMALA/params.csv", row.names=FALSE)





