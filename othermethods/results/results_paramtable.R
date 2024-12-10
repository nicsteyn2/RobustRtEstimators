
# CAUTION: This file makes assumptions about the prior hyperparameters. See readme for more details.

rm(list=ls())


# Pre-allocate priors dataframe
df_priors = data.frame()

# Load posterior values
df_posterior = bind_rows(
  read.csv("othermethods/EpiLPSMAP/params.csv") %>% mutate(method = "EpiLPS (MAP)"),
  read.csv("othermethods/EpiLPSMALA/params.csv")%>% mutate(method = "EpiLPS (MALA)"),
  read.csv("othermethods/EpiFilterSmoothing/params.csv") %>% mutate(method = "EpiFilter (smoothing)"),
  read.csv("othermethods/EpiNow2/params.csv") %>% mutate(method = "EpiNow2"),
  read.csv("othermethods/EpiNow2/params_invgamma.csv") %>% mutate(method = "EpiNow2 (inv-gam)"),
  read.csv("othermethods/rtestim/params.csv") %>% mutate(method = "rtestim")
) %>% mutate(type="posterior")



# --------------- Fetch priors -------------- 

# EpiFilter 
df_priors = rbind(df_priors, data.frame(
  method="EpiFilter (smoothing)", 
  parameter="eta",
  mean=0.5,
  mode=NaN,
  lower=0.025,
  upper=0.975
))


# EpiLPS smoothness parameter lambda
a_delta <- 10
b_delta <- 10
phi <- 2

# Calculate P(lambda | phi, a_delta, b_delta) numerically
plambda = function(lambda, phi, a_delta, b_delta) {
  integrand = function(delta) { 
    dgamma(lambda, shape=phi/2, rate=(phi*delta)/2) * dgamma(delta, shape=a_delta, rate=b_delta)
  }
  return(integrate(integrand, lower=0, upper=Inf)$value)
}

# Find the PDF and CDF on a grid of lambda values
x = seq(0, 1000, length.out=100000)
fx = sapply(x, plambda, phi=phi, a_delta=a_delta, b_delta=b_delta)
Fx = cumsum(fx) / sum(fx)

# Store mean, mode, etc
df_priors = rbind(df_priors, data.frame(
  method="EpiLPS (MAP)", 
  parameter="lambda",
  mean=b_delta/(a_delta-1),
  mode=x[which.max(fx)],
  lower=0,
  upper=x[which.min(abs(Fx-0.95))]
))

df_priors = rbind(df_priors, data.frame(
  method="EpiLPS (MALA)", 
  parameter="lambda",
  mean=b_delta/(a_delta-1),
  mode=x[which.max(fx)],
  lower=0,
  upper=x[which.min(abs(Fx-0.95))]
))


# EpiLPS overdispersion parameter rho*

a_rho <- 1e-04
b_rho <- 1e-04

df_priors = rbind(df_priors, data.frame(
  method="EpiLPS (MAP)",
  parameter="rho",
  mean=a_rho/b_rho,
  mode=ifelse(a_rho > 1, (a_rho-1)/b_rho, 0),
  lower=ifelse(a_rho > 1, qgamma(0.025, shape=a_rho, rate=b_rho), 0),
  upper=ifelse(a_rho > 1, qgamma(0.975, shape=a_rho, rate=b_rho), qgamma(0.95, shape=a_rho, rate=b_rho))
))

df_priors = rbind(df_priors, data.frame(
  method="EpiLPS (MALA)",
  parameter="rho",
  mean=a_rho/b_rho,
  mode=ifelse(a_rho > 1, (a_rho-1)/b_rho, 0),
  lower=ifelse(a_rho > 1, qgamma(0.025, shape=a_rho, rate=b_rho), 0),
  upper=ifelse(a_rho > 1, qgamma(0.975, shape=a_rho, rate=b_rho), qgamma(0.95, shape=a_rho, rate=b_rho))
))

# EpiNow2 gaussian process lengthscale (log-normal prior)
ell_mean = 21
ell_sd = 7

meanlog = log((ell_mean^2)/sqrt(ell_sd^2 + ell_mean^2))
sdlog = sqrt(log(1 + (ell_sd^2)/(ell_mean^2)))

df_priors = rbind(df_priors, data.frame(
  method = "EpiNow2",
  parameter = "ell",
  mean = exp(meanlog + (sdlog^2)/2),
  mode = exp(meanlog - sdlog^2),
  lower = exp(qnorm(0.025, meanlog, sdlog)),
  upper = exp(qnorm(0.975, meanlog, sdlog))
))

# EpiNow2 gaussian process lengthscale (inverse gamma prior)
library(extraDistr)
shape = 1.499007
scale = 0.057277 * 60

df_priors = rbind(df_priors, data.frame(
  method = "EpiNow2 (inv-gam)",
  parameter = "ell",
  mean = scale/(shape-1),
  mode = scale/(shape+1),
  lower = qinvgamma(0.025, shape, scale),
  upper = qinvgamma(0.975, shape, scale)
))

# 
# # EpiNow2 reporting overdispersion
# # If 1/sqrt{phi} ~ HalfNormal(0, 1), then the PDF of phi is:
# 
# f = function(phi) { sqrt(2/pi) * phi^(-3/2) * exp(-1/(2*phi)) }
#
# This is very poorly behaved numerically, so we avoid consider this in detail.
# Code to do so is attached but commented out below.
# 
# x = seq(0, 100000, length.out=10000000)
# fx = sapply(x, f)
# fx[1] = 0
# Fx = cumsum(fx) / sum(fx)
# 
# df_priors = rbind(df_priors, data.frame(
#   method = "EpiNow2",
#   parameter = "phi",
#   mean = sum(x * fx)/sum(fx),
#   mode = x[which.max(fx)],
#   lower = x[which.min(abs(Fx-0.025))],
#   upper = x[which.min(abs(Fx-0.975))]
# ))




# combine it all
df_all = rbind(df_posterior, df_priors %>% mutate(type="prior", simulation="prior"))

smoothing_table = df_all %>%
  filter(parameter %in% c("eta", "lambda", "ell")) %>%
  mutate(value = paste0(signif(mode, digits=3), " (", signif(lower, digits=3), ", ", signif(upper, digits=3), ")")) %>%
  pivot_wider(names_from="simulation", values_from="value", id_cols=c("method", "parameter")) %>%
  dplyr::select(all_of(c("method", "prior", "Random walk", "Sinusoidal", "Step change"))) %>%
  mutate(method = factor(method, levels=c("EpiFilter (smoothing)", "EpiNow2", "EpiNow2 (inv-gam)", "EpiLPS (MAP)", "EpiLPS (MALA)", "rtestim"))) %>%
  arrange(method) %>%
  mutate(row = sprintf("%s & %s & %s & %s & %s \\\\",
                       method,
                       prior,
                       `Random walk`,
                       `Sinusoidal`,
                       `Step change`)) %>%
  pull(row) %>%
  cat(sep = "\n")


noise_table = df_all %>%
  filter(parameter %in% c("rho", "phi")) %>%
  mutate(value = paste0(signif(mode, digits=3), " (", signif(lower, digits=3), ", ", signif(upper, digits=3), ")")) %>%
  pivot_wider(names_from="simulation", values_from="value", id_cols=c("method", "parameter")) %>%
  dplyr::select(all_of(c("method", "prior", "Random walk", "Sinusoidal", "Step change"))) %>%
  mutate(method = factor(method, levels=c("EpiFilter (smoothing)", "EpiNow2", "EpiNow2 (inv-gam)", "EpiLPS (MAP)", "EpiLPS (MALA)", "rtestim"))) %>%
  arrange(method) %>%
  mutate(row = sprintf("%s & %s & %s & %s & %s \\\\",
                       method,
                       prior,
                       `Random walk`,
                       `Sinusoidal`,
                       `Step change`)) %>%
  pull(row) %>%
  cat(sep = "\n")







