
library(MASS)

observedIncidenceCurve = function(LPSfit, alpha = 0.05, n_samples = 1000) {
  
  # Extract key params
  numdays <- length(LPSfit$incidence)
  Time  <- seq_len(numdays)
  thetahat <- LPSfit$thetahat
  Sighat <- LPSfit$Sighat
  K <- length(thetahat)
  B <- Rcpp_KercubicBspline(x = seq_len(numdays), lower = 1, upper = numdays, K = K)
  
  # Extract the overdispersion parameter
  rho <- LPSfit$NegBinoverdisp
  
  # Sample theta from the (approximated) posterior
  theta_samples <- MASS::mvrnorm(1000, mu = thetahat, Sigma = Sighat)
  
  # Compute mu(t) for each sampled theta
  mu_samples <- apply(theta_samples, 1, function(theta) { exp(B %*% theta) })
  
  # Generate observables y_t for each sampled mu(t)
  y_samples <- apply(mu_samples, 2, function(mu) { rnbinom(numdays, mu = mu, size = rho) })
  
  # Compute central estimates and credible intervals
  musmooth <- rowMeans(mu_samples)
  lower <- apply(y_samples, 1, quantile, probs = alpha / 2)
  upper <- apply(y_samples, 1, quantile, probs = 1 - alpha / 2)
  
  # Return the smoothed incidence curve and credible intervals
  df_out = data.frame(Time, musmooth, lower, upper)
  colnames(df_out) = c("Time", "I", "Lower", "Upper")
  return(df_out)
}
