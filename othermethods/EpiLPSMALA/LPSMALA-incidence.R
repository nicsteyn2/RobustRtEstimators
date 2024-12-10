
# A custom function based on https://github.com/oswaldogressani/EpiLPS/blob/main/R/epicurve.R
# that returns the smoothed incidence curve and credible intervals rather than
# a ggplot object.

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
  
  # Sample theta and rho from the (sampled) posterior
  theta_samples <- LPSfit$MCMC$theta_mcmc
  rho_samples <- LPSfit$MCMC$rho_mcmc
  
  # Compute mu(t) for each sampled theta
  mu_samples <- apply(theta_samples, 1, function(theta) { exp(B %*% theta) })
  
  # Generate observables y_t for each sampled mu(t)
  y_samples = matrix(nrow=numdays, ncol=dim(mu_samples)[2])
  for (jj in seq(1, dim(mu_samples)[2])) {
    y_samples[,jj] = rnbinom(numdays, mu=mu_samples[,jj], size=rho_samples[jj])
  }
  
  # Compute central estimates and credible intervals
  musmooth <- rowMeans(mu_samples)
  lower <- apply(y_samples, 1, quantile, probs = alpha / 2)
  upper <- apply(y_samples, 1, quantile, probs = 1 - alpha / 2)
  
  # Return the smoothed incidence curve and credible intervals
  df_out = data.frame(Time, musmooth, lower, upper)
  colnames(df_out) = c("Time", "I", "Lower", "Upper")
  return(df_out)
}
