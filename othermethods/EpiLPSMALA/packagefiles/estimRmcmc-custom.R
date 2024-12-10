

estimRmcmc_custom <- function(incidence, si, K = 30, dates = NULL, niter = 5000,
                       burnin = 2000, CoriR = FALSE, WTR = FALSE,
                       priors = Rmodelpriors(), progressbar = TRUE){
  tic <- proc.time()             # Clock starts ticking
  y <- KerIncidCheck(incidence)  # Run checks on case incidence vector
  n <- length(y)                 # Total number of days of the epidemic
  simax <- length(si)            # Length of serial interval distribution
  B <- Rcpp_KercubicBspline(seq_len(n), lower = 1, upper = n, K = K) # C++ call
  D <- diag(K)
  penorder <- 2
  for(k in 1:penorder){D <- diff(D)}
  P <- t(D) %*% D           # Penalty matrix of dimension c(K,K)
  P <- P + diag(1e-06, K)   # Perturbation to ensure P is full rank
  
  # Prior specification
  priorspec <- priors
  a_delta <- priorspec$a_delta
  b_delta <- priorspec$b_delta
  phi <- priorspec$phi
  a_rho <- priorspec$a_rho
  b_rho <- priorspec$b_rho
  
  # Minus log-posterior of hyperparameter vector
  logphyper <- function(x) {
    w <- x[1] # w = log(rho)
    v <- x[2] # v = log(lambda)
    
    # Laplace approximation
    LL <- Rcpp_KerLaplace(theta0 = rep(1.5,K), exp(w), exp(v), K,
                          KerPtheta(Dobs = y, BB = B, Pen = P)$Dlogptheta,
                          KerPtheta(Dobs = y, BB = B, Pen = P)$D2logptheta)
    thetastar <- as.numeric(LL$Lapmode)
    logdetSigstar <- Re(LL$logdetSigma)
    
    equal <- (-1) * (0.5 * logdetSigstar + 0.5 * (K + phi) * v + a_rho * w -
                       (0.5 * phi + a_delta) * log(0.5 * phi * exp(v) + b_delta) +
                       KerLikelihood(Dobs = y, BB = B)$loglik(thetastar, exp(w)) -
                       0.5 * exp(v) * sum((thetastar * P) %*% thetastar) -
                       b_rho * exp(w))
    return(equal)
  }
  
  # Maximization of the posterior hyperparameters
  optim_method <- "NelderMead"
  optimhyper <- stats::optim(c(1, 5), fn = logphyper)
  hypermap <- optimhyper$par
  optimconverged <- (optimhyper$convergence == 0)
  disphat <- exp(hypermap[1])
  lambhat <- exp(hypermap[2])
  Lap_approx <- Rcpp_KerLaplace(theta0 = rep(1.5,K), disphat, lambhat, K,
                                KerPtheta(Dobs = y, BB = B, Pen = P)$Dlogptheta,
                                KerPtheta(Dobs = y, BB = B, Pen = P)$D2logptheta)
  thetahat <- as.numeric(Lap_approx$Lapmode)
  muhat <- as.numeric(exp(B %*% thetahat))
  Sighat <- Lap_approx$Lapvar
  
  # Call Metropolis-adjusted Langevin algorithm
  if(isTRUE(progressbar)){
    cat(paste0("Metropolis-adjusted Langevin algorithm running for ",niter,
               " iterations \n"))
    progbar <- utils::txtProgressBar(min = 1, max = niter, initial = 1,
                                     style = 3, char =">")
  } else{
    progbar <- NULL
  }
  MCMCout <- KerMCMC(Dobs = incidence, BB = B, Pen = P, Covar = Sighat,
                     thetaoptim = thetahat, penoptim = lambhat,
                     overdispoptim = disphat, progress = progbar,
                     priors = priorspec)$MALA(M=niter)
  # Chain extraction
  lambdaMALA <- coda::as.mcmc(MCMCout$lambda_mcmc[(burnin+1):niter])
  deltaMALA <- coda::as.mcmc(MCMCout$delta_mcmc[(burnin+1):niter])
  rhoMALA <- coda::as.mcmc(MCMCout$rho_mcmc[(burnin+1):niter])
  thetaMALA <- coda::as.mcmc(MCMCout$theta_mcmc[(burnin+1):niter,])
  accept_mcmc <- MCMCout$accept_rate
  
  # Point estimation
  thetahat_mcmc  <- colMeans(thetaMALA)
  lambdahat_mcmc <-  mean(lambdaMALA)
  rhohat_mcmc <- mean(rhoMALA)
  muMALA_mcmc <- matrix(0, nrow = (niter - burnin), ncol = n)
  for(j in 1:(niter-burnin)){
    muMALA_mcmc[j,] <- as.numeric(exp(B %*% thetaMALA[j, ]))
  }
  mu_estim <- colMeans(muMALA_mcmc)
  
  if(is.null(dates)) {
    Time <- seq_len(n)
  } else{
    Time <- dates
  }
  
  RLPS <- data.frame(matrix(0, nrow = n, ncol = 10))
  colnames(RLPS) <- c("Time", "R", "Rsd", "Rq0.025", "Rq0.05","Rq0.25",
                      "Rq0.50","Rq0.75", "Rq0.95", "Rq0.975")
  RLPS$Time <- Time
  HPD90_Rt <- data.frame(matrix(0, nrow = n, ncol = 2))
  HPD95_Rt <- data.frame(matrix(0, nrow = n, ncol = 2))
  colnames(HPD90_Rt) <- c("HPD90_low","HPD90_up")
  colnames(HPD95_Rt) <- c("HPD95_low","HPD95_up")
  rownames(HPD90_Rt) <- Time
  rownames(HPD95_Rt) <- Time
  
  for(j in 1:n){
    CppCall <- as.numeric(Rcpp_KerRpostmcmc(t = j, BB = B, sinter = si,
                                            thetasample = thetaMALA))
    RLPS$R[j] <- mean(CppCall)
    RLPS$Rsd[j] <- stats::sd(CppCall)
    RLPS[j, (4:10)] <- stats::quantile(CppCall,
                                       probs = c(0.025, 0.05, 0.25, 0.50, 0.75,
                                                 0.95, 0.975))
    HPD90_Rt[j, ] <- coda::HPDinterval(coda::as.mcmc(CppCall), prob = 0.90)
    HPD95_Rt[j, ] <- coda::HPDinterval(coda::as.mcmc(CppCall), prob = 0.95)
  }
  
  if (CoriR == TRUE) {# Use Cori method with weekly sliding windows
    RCori <- KerCori(Dobs = incidence, sinter = si)
  } else{
    RCori <- "Not called"
  }
  
  if(WTR == TRUE){# Use Wallinga-Teunis method to estimate R
    RWT <- KerWT(Dobs = incidence, sinter = si)
  } else{
    RWT <- "Not called"
  }
  
  toc <- proc.time() - tic
  toc <- round(toc[3], 3)
  
  #-- Output results in a list
  outputlist <- list(incidence = y, si = si, RLPS = RLPS,
                     thetahat = thetahat,
                     Sighat = Sighat,
                     RCori = RCori, RWT = RWT,
                     LPS_elapsed = toc, penparam = lambdahat_mcmc, K = K,
                     NegBinoverdisp = rhohat_mcmc,
                     optimconverged = optimconverged,
                     method = "LPSMALA",
                     optim_method = optim_method,
                     HPD90_Rt = HPD90_Rt,
                     HPD95_Rt = HPD95_Rt,
                     MCMC = MCMCout)
  
  attr(outputlist, "class") <- "Rt"
  outputlist
}