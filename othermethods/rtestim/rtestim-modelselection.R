
# An example of using CRPS for rtestim model selection instead of cross-validation
# very early work in progress

rm(list=ls())

library(tidyverse)
library(foreach)
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


# Extract incidence data for current simulation
ii = 1
Ct = df %>% filter(simulation==simulations[ii]) %>% arrange(t) %>% pull(Cases)

lambda_vals = seq(0.1, 10, 0.2)
CRPS_vals = numeric(length(lambda_vals))
rtestim_fit = estimate_rt(Ct, lambda=lambda_vals, dist_gamma=c(2.3669, 2.7463))

for (jj in seq(1, length(lambda_vals))) {
  
  # Calculate CRPS
  probs_in = seq(0.02,0.98,0.02)
  rtestim_I_levels = confband(rtestim_fit, lambda=lambda_vals[jj], type="Yt", level=probs_in) %>% dplyr::select(-fit)
  probs_out = setdiff(seq(0.01, 0.99, 0.01), 0.5)
  
  CRPSall = numeric(length(Ct))
  CRPSall[1] = NA
  for (tt in seq(2,length(Ct))) {
    CRPSall[tt] = CRPStrapezoidal(probs_out, as.numeric(rtestim_I_levels[tt,]), Ct[tt])
  }
  CRPS_vals[jj] = mean(CRPSall[10:length(Ct)])
  
}

plot(lambda_vals, CRPS_vals)
