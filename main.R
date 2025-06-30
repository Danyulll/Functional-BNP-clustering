##################################
#### Gaussian Process Example ####
##################################

# Example: clustering two distinct Gaussian‐process–generated functional datasets
# ------------------------------------------------------------------------------

# 1. Load libraries and your FBNP implementation
library(fda)
library(MASS)        # for mvrnorm()
library(fdaoutlier)  # for smoothing helper
library(coda)        # for MCMC diagnostics

source("Tools/FBNP/FBNP.R")
source("Tools/FBNP/FBNP_hyper.R")
source("Tools/hyperparameters.R")
source("Tools/Smoothing.R")

set.seed(12345)

# squared‐exponential kernel
se_kernel <- function(x1, x2, l, sigma_f) {
  dist2 <- outer(x1, x2, FUN = function(a,b) (a - b)^2)
  sigma_f^2 * exp(-0.5 * dist2 / l^2)
}

# time grid
T <- 512
time.grid <- seq(0, 1, length.out = T)

# Group 1: zero‐mean, very short length‐scale, low amplitude, + nugget noise
n1      <- 80
ℓ1      <- 0.02           # very rough
σf1     <- 0.5            # low amplitude
Σ1      <- se_kernel(time.grid, time.grid, l = ℓ1, sigma_f = σf1)
X1_core <- mvrnorm(n = n1, mu = rep(0, T), Sigma = Σ1)
# add independent noise so they really look “jittery”
noise1  <- matrix(rnorm(n1 * T, sd = 0.1), nrow = n1, ncol = T)
X1      <- X1_core + noise1

# Group 2: sinusoidal mean, very long length‐scale, high amplitude, no extra noise
n2      <- 20
ℓ2      <- 0.4            # very smooth
σf2     <- 3              # high amplitude
mu2     <- sin(2 * pi * time.grid)
Σ2      <- se_kernel(time.grid, time.grid, l = ℓ2, sigma_f = σf2)
X2      <- mvrnorm(n = n2, mu = mu2, Sigma = Σ2)

# Combine and rescale
X    <- rbind(X1, X2)
X    <- X / max(abs(X))
true_labels <- c(rep(1, n1), rep(2, n2))

known <- true_labels
known[which(known == 2)] <- NA
known[40:length(known)] <- NA

# 4. Smooth curves via basis expansion
smoothing_list <- smoothing(
  X           = X,
  step        = 1,
  nbasis      = 20,
  spline_order= 4
)
basis.t  <- t(eval.basis(smoothing_list$time.grid, smoothing_list$basis))
X_smooth <- smoothing_list$beta %*% basis.t

par(mfrow = c(1,2))
matplot(time.grid, t(X),        type = "l", col = true_labels,
        main = "Original")
matplot(time.grid, t(X_smooth), type = "l", col = true_labels,
        main = "Smoothed")


mass <- 2.0

hyperparam_base_fixed <- list(
  m0 = rep(0, ncol(smoothing_list$beta)),
  k0 = 0.01,  # Very low confidence
  nu0 = ncol(smoothing_list$beta) + 1,
  Lambda0 = diag(ncol(smoothing_list$beta)) * 10,  # More diffuse
  c = 1.01,   # Very weakly informative
  d = 0.01    # Very small rate
)

# Less informative supervised priors
group1_cov_beta <- cov(smoothing_list$beta[group1_indices, ])
hyperparam_sup_fixed <- list(
  m0 = group1_mean_beta,
  k0 = 2,     # Lower confidence
  nu0 = ncol(smoothing_list$beta) + 2,
  Lambda0 = group1_cov_beta * 5,  # More uncertainty
  c = 2,      # Less informative shape
  d = 1 * group1_mean_var  # Adjusted rate
)

out<- FBNP_slice(
  n_iter = 3000,
  burnin = 1000,
  max_clusters = 50,
  mass_init = 5.0,           # Initial mass value
  mass_prior = c(1, 1),      # Weak Gamma(1,1) prior
  update_mass = TRUE,        # Enable mass updating
  smoothing = smoothing_list,
  hyperparam_base = hyperparam_base,
  hyperparam_sup = hyperparam_sup,
  known_labels = known
)

# Quick check of results
K_mode_fixed <- apply(out$K, 2, function(x) {
  tbl <- table(x)
  as.numeric(names(tbl)[which.max(tbl)])
})

par(mfrow = c(1, 2))
matplot(time.grid, t(X), type = "l", col = true_labels, lty = 1,
        main = "True Clusters")
matplot(time.grid, t(X), type = "l", col = factor(K_mode_fixed), lty = 1,
        main = "Fixed Model Predictions")


out$nclusters |> hist()

out$phi[,1,1] |> as.mcmc() |> traceplot()
out$mu_coef[,1,1] |> as.mcmc() |> traceplot()
out$mu[,1,1] |> as.mcmc() |> traceplot()
out$mass |> as.mcmc() |> traceplot()


out$phi[,1,1] |> as.mcmc() |> acfplot()
out$mu_coef[,1,1] |> as.mcmc() |> acfplot()
out$mu[,1,1] |> as.mcmc() |> acfplot()
out$mass |> as.mcmc() |> acfplot()




###################################
####### FDA Outlier Example #######
###################################


## Model 1

# 1. Load libraries and your FBNP implementation
library(fda)
library(MASS)        # for mvrnorm()
library(fdaoutlier)  # for smoothing helper
library(coda)        # for MCMC diagnostics

source("Tools/FBNP/FBNP.R")
source("Tools/FBNP/FBNP_hyper.R")
source("Tools/hyperparameters.R")
source("Tools/Smoothing.R")

#### DATA #### -------------------------------------------------------------------------------
set.seed(87360935)
library(fdaoutlier)
sim <- simulation_model1(n=200,p=512,outlier_rate = 0.05)
X <- sim$data
labs <- rep(1,nrow(X))
labs[sim$true_outliers] <- 2

rescale <- max(X)
X <- X/rescale 

# cut x-axis
par(mfrow=c(1,1))
matplot(t(X), type='l',col=labs)

smoothing_list <- smoothing(X = X, 
                            step = 1, 
                            nbasis = 20, 
                            spline_order = 4)

basis.t <- t(eval.basis(smoothing_list$time.grid, smoothing_list$basis))
X_smooth <- smoothing_list$beta %*% basis.t

par(mfrow=c(1,2))
matplot(t(X), type='l', lwd=1, lty=1,col=labs,main="Original")
matplot(smoothing_list$time.grid, t(X_smooth), type='l',col=labs,main="Smoothed")

hyperparam_base_fixed <- list(
  m0 = rep(0, ncol(smoothing_list$beta)),
  k0 = 0.01,  # Very low confidence
  nu0 = ncol(smoothing_list$beta) + 1,
  Lambda0 = diag(ncol(smoothing_list$beta)) * 10,  # More diffuse
  c = 1.01,   # Very weakly informative
  d = 0.01    # Very small rate
)

# Less informative supervised priors
group1_cov_beta <- cov(smoothing_list$beta[group1_indices, ])
hyperparam_sup_fixed <- list(
  m0 = group1_mean_beta,
  k0 = 2,     # Lower confidence
  nu0 = ncol(smoothing_list$beta) + 2,
  Lambda0 = group1_cov_beta * 5,  # More uncertainty
  c = 2,      # Less informative shape
  d = 1 * group1_mean_var  # Adjusted rate
)

known <- labs
known[which(known == 2)] <- NA
known[50:length(known)] <- NA
# known <- rep(NA, length(labs))

out<- FBNP_hyper(
  n_iter = 8000,
  burnin = 5000,
  M = 10,
  mass_init = 100.0,           # Initial mass value
  mass_prior = c(1, 1),      # Weak Gamma(1,1) prior
  update_mass = TRUE,        # Enable mass updating
  smoothing = smoothing_list,
  hyperparam_base = hyperparam_base,
  hyperparam_sup = hyperparam_sup,
  known_labels = known
)



K_mode_fixed <- apply(out$K, 2, function(x) {
  tbl <- table(x)
  as.numeric(names(tbl)[which.max(tbl)])
})

par(mfrow = c(1, 2))
matplot(time.grid, t(X), type = "l", col = true_labels, lty = 1,
        main = "True Clusters")
matplot(time.grid, t(X), type = "l", col = factor(K_mode_fixed), lty = 1,
        main = "Fixed Model Predictions")


out$nclusters |> hist()

out$phi[,1,1] |> as.mcmc() |> traceplot()
out$mu_coef[,1,1] |> as.mcmc() |> traceplot()
out$mu[,1,1] |> as.mcmc() |> traceplot()
out$mass |> as.mcmc() |> traceplot()


out$phi[,1,1] |> as.mcmc() |> acfplot()
out$mu_coef[,1,1] |> as.mcmc() |> acfplot()
out$mu[,1,1] |> as.mcmc() |> acfplot()
out$mass |> as.mcmc() |> acfplot()




