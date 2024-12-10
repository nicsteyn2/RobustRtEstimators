
rm(list=ls())

library(tidyverse)
library(foreach)
library(EpiNow2)
library(rstan)

# Define CRPS estimator using MCMC samples from Stan
CRPS = function(samples, obs) {
  
  mean_absolute_deviation = mean(abs(samples - obs))
  mean_pairwise_absolute_diffs = 0.5 * (1 / length(samples)^2) * sum(outer(samples, samples, FUN=function(x,y) abs(x-y)))
  return(mean_absolute_deviation - mean_pairwise_absolute_diffs)
  
}

# Function to transform phi onto natural scale
inv_square = function(x){ return(1/(x^2)) }

# Load simulated data
df = read.csv("othermethods/exampledata.csv") %>% arrange(t) %>% mutate(date = as.Date("2024-01-01") + t-1)
simulations = c("Random walk", "Sinusoidal", "Step change")

# Pre-allocate dataframe to store parameter estimates
df_params = data.frame()

df_results_invgam = foreach(ii = seq(1, length(simulations)), .combine=rbind) %do% {
  
  # Extract incidence data for current simulation and setup for EpiNow2 (including making a fake date vector)
  Ct = df %>% filter(simulation==simulations[ii]) %>% rename(confirm=Cases) %>% dplyr::select(date, confirm)
  
  # Fit the model
  EpiNow2Fit = estimate_infections(
    data = Ct,
    stan = stan_opts(return_fit=TRUE, cores=4),
    gp = gp_opts(ls_sd=0), # This forces an inverse-gamma prior on the length-scale
    gt_opts(dist=Gamma(shape=2.3669, scale=2.7463, max=21)),
    delay_opts(dist=Fixed(0)), # No delay
    horizon=0, # No forecasting
    CrIs = c(0.975)
  )
  
  # Extract posterior samples
  df_R = EpiNow2Fit$summarised %>% filter(variable=="R") %>%
    mutate(t = as.numeric(date - as.Date("2024-01-01") + 1)) %>%
    rename(Mean=mean, Lower=lower_97.5, Upper=upper_97.5) %>%
    dplyr::select(t, Mean, Lower, Upper, variable)
  
  df_I = EpiNow2Fit$summarised %>% filter(variable=="reported_cases") %>%
    mutate(t = as.numeric(date - as.Date("2024-01-01") + 1)) %>%
    rename(Mean=mean, Lower=lower_97.5, Upper=upper_97.5) %>%
    dplyr::select(t, Mean, Lower, Upper) %>%
    mutate(variable = "I")
  
  # Extract samples for CRPS calculation
  samples = EpiNow2Fit$samples %>% filter(variable == "reported_cases") %>% dplyr::select(time, value)
  
  # Calculate CRPS
  CRPSvalues = numeric(nrow(Ct))
  for (t in seq(1, nrow(Ct))) {
    CRPSvalues[t] = CRPS(samples$value[samples$time == t], Ct$confirm[t])
  }
  
  # Append CRPS to dataframes
  df_R$CRPS = NA
  df_I$CRPS = CRPSvalues
  
  # Extract MCMC samples
  samples = extract(EpiNow2Fit$fit)
  
  # Summarise lengthscale posterior
  kderho = density(samples$rho, from=0)
  mode_rho = kderho$x[which.max(kderho$y)]
  mean_rho = mean(samples$rho)
  lower_rho = quantile(samples$rho, 0.025)
  upper_rho = quantile(samples$rho, 0.975)
  print(paste("MAP lengthscale:", mode_rho, "95% CI:", lower_rho, upper_rho))
  
  # Summarise reporting overdispersion posterior 
  phi = inv_square(samples$rep_phi)
  kde_phi = density(phi, from=0)
  mode_phi = kde_phi$x[which.max(kde_phi$y)]
  mean_phi = mean(phi)
  lower_phi = quantile(phi, 0.025)
  upper_phi = quantile(phi, 0.975)
  print(paste("MAP reporting overdisp:", mode_phi, "95% CI:", lower_phi, upper_phi))
  
  # Store parameter estimates
  df_params = rbind(df_params, data.frame(
    simulation = c(simulations[ii], simulations[ii]),
    parameter = c("ell", "phi"),
    mean = c(mean_rho, mean_phi),
    mode = c(mode_rho, mode_phi),
    lower = c(lower_rho, lower_phi),
    upper = c(upper_rho, upper_phi)
  ))
  
  # Return results
  return(rbind(df_R, df_I) %>% mutate(simulation=simulations[ii]))
  
}

# Append the true values
df_truth = df %>% rename(R = TrueRt, I = Cases) %>% pivot_longer(cols=c("R", "I"), names_to="variable", values_to="TrueValue")
df_results_invgam = left_join(df_results_invgam, df_truth, by=c("t", "simulation", "variable"))

# Save results
write.csv(df_results_invgam, "othermethods/EpiNow2/results_invgamma.csv", row.names=FALSE)
write.csv(df_params, "othermethods/EpiNow2/params_invgamma.csv", row.names=FALSE)
