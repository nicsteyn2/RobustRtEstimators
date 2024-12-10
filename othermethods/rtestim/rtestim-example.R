
# remotes::install_github("dajmcdon/rtestim") # installed on 7 Dec using commit ID cf47415

rm(list=ls())

library(tidyverse)
library(rtestim)

# Load simulated data
df = read.csv("othermethods/exampledata.csv")
simulations = c("Random walk", "Sinusoidal", "Step change")

# Define trapezoidal numerical CRPS function
CRPStrapezoidal = function(Fx, x, obs) {
  
  indic = ifelse(x >= obs, 1, 0)
  CRPS = 0
  for (ii in seq(1, length(x)-1)) {
    CRPS = CRPS + (0.5 * (x[ii+1] - x[ii]) * ( (Fx[ii+1] - indic[ii+1])^2 + (Fx[ii] - indic[ii])^2 ) )
  }
  return(CRPS)
  
}

# Specify choice of lambda
lambda = "lambda.min"

# Pre-allocate dataframe to store parameter estimates
df_params = data.frame()


df_results = foreach (ii = seq(1, length(simulations)), .combine=rbind) %do% {
  
  # Extract incidence data for current simulation
  Ct = df %>% filter(simulation==simulations[ii]) %>% arrange(t) %>% pull(Cases)

  # Use cross-validation to find the optimal value of lambda
  rtestim_cv = cv_estimate_rt(Ct, dist_gamma=c(2.3669, 2.7463))
  
  # Calculate Rt estimates
  rtestim_R = confband(rtestim_cv, lambda=lambda, type="Rt", level=0.95)
  colnames(rtestim_R) = c("Mean", "Lower", "Upper")
  
  # Calculate It estimates
  rtestim_I = confband(rtestim_cv, lambda=lambda, type="Yt", level=0.95)
  colnames(rtestim_I) = c("Mean", "Lower", "Upper")
  
  # Calculate CRPS
  probs_in = seq(0.02,0.98,0.02)
  rtestim_I_levels = confband(rtestim_cv, lambda=lambda, type="Yt", level=probs_in) %>% dplyr::select(-fit)
  probs_out = setdiff(seq(0.01, 0.99, 0.01), 0.5)
  
  CRPSall = numeric(length(Ct))
  CRPSall[1] = NA
  for (tt in seq(2,length(Ct))) {
    CRPSall[tt] = CRPStrapezoidal(probs_out, as.numeric(rtestim_I_levels[tt,]), Ct[tt])
  }
  
  rtestim_R$CRPS = NA
  rtestim_I$CRPS = CRPSall
  
  rtestim_R$t = 1:100
  rtestim_I$t = 1:100
  
  rtestim_R$variable = "R"
  rtestim_I$variable = "I"
  
  lambda = rtestim_cv$lambda.min
  print(paste0("Simulation: ", simulations[ii], ", lambda = ", lambda))
  df_params = rbind(df_params, data.frame(
    simulation=simulations[ii],
    parameter="lambda",
    mode=lambda
  ))
  
  return(rbind(rtestim_R, rtestim_I) %>% mutate(simulation=simulations[ii]))
  
}

# Append true values
df_truth = df %>% rename(R = TrueRt, I=Cases) %>% pivot_longer(cols=c("R", "I"), names_to="variable", values_to="TrueValue")
df_results = left_join(df_results, df_truth, by=c("t", "variable", "simulation"))
write.csv(df_results, "othermethods/rtestim/results.csv", row.names=FALSE)
write.csv(df_params, "othermethods/rtestim/params.csv", row.names=FALSE)


