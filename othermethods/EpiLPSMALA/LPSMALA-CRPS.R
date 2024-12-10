

computeCRPSt = function(samples, obs) {
  
  mean_absolute_deviation = mean(abs(samples - obs))
  mean_pairwise_absolute_diffs = 0.5 * (1 / length(samples)^2) * sum(outer(samples, samples, FUN=function(x,y) abs(x-y)))
  return(mean_absolute_deviation - mean_pairwise_absolute_diffs)
  
}


computeCRPS = function(LPSfit) {
  
  # Extract key params
  numdays <- length(LPSfit$incidence)
  thetahat <- LPSfit$thetahat
  K <- length(thetahat)
  B <- Rcpp_KercubicBspline(x = seq_len(numdays), lower = 1, upper = numdays, K = K)
  
  # Extract samples of incidence
  theta_samples <- LPSfit$MCMC$theta_mcmc
  rho_samples <- LPSfit$MCMC$rho_mcmc
  mu_samples <- apply(theta_samples, 1, function(theta) { exp(B %*% theta) })
  y_samples = matrix(nrow=numdays, ncol=dim(mu_samples)[2])
  for (jj in seq(1, dim(mu_samples)[2])) {
    y_samples[,jj] = rnbinom(numdays, mu=mu_samples[,jj], size=rho_samples[jj])
  }
  
  # Compute CRPS for each day
  CRPS = numeric(numdays)
  pb <- txtProgressBar(min = 1, max = numdays, style = 3)
  for (tt in seq(1, numdays)) {
    CRPS[tt] = computeCRPSt(y_samples[tt,], LPSfit$incidence[tt])
    setTxtProgressBar(pb, tt)
  }
  close(pb)
  
  return(CRPS)
  
}