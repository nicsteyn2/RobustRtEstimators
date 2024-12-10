
library(scoringRules)
library(MASS)

options(hypergeo_tol = 1e-10, hypergeo_maxiter = 5000)

# Function to compute CRPS with uncertainty in mu(t)
computeCRPS <- function(LPSfit, n_samples = 1000) {
  
  # Extract parameters
  numdays <- length(LPSfit$incidence)
  thetahat <- LPSfit$thetahat
  Sighat <- LPSfit$Sighat
  K <- length(thetahat)
  B <- Rcpp_KercubicBspline(x = seq_len(numdays), lower = 1, upper = numdays, K = K)
  rho <- LPSfit$NegBinoverdisp
  
  # Sample theta from posterior
  theta_samples <- MASS::mvrnorm(n_samples, mu = thetahat, Sigma = Sighat)
  
  # Compute mu(t) samples
  mu_samples <- apply(theta_samples, 1, function(theta) { exp(B %*% theta) })
  
  # Initialize vector to store CRPS values for each time point
  crps_values <- numeric(numdays)
  crps_values[1] = NaN
  
  # Loop over each time point
  pb <- txtProgressBar(min = 2, max = numdays, style = 3)
  for (t in seq(2, numdays)) {
    # Extract sampled mu(t)
    mu_t_samples <- mu_samples[t, ]
    
    # Compute predictive distribution
    crps_t <- mean(sapply(mu_t_samples, function(mu_t) {
      size <- rho
      prob <- rho / (rho + mu_t)
      result <- crps_nbinom(y = LPSfit$incidence[t], size = size, prob = prob)
    }))
    
    crps_values[t] <- crps_t
    setTxtProgressBar(pb, t)
  }
  close(pb)
  
  return(crps_values)
}
