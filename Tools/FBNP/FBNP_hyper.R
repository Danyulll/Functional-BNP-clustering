library(invgamma)        # for invgamma::invgamma::rinvgamma()
library(LaplacesDemon)   # for rinvwishart()
library(MASS)            # for mvrnorm()
library(pbmcapply)       # for progressBar()

FBNP_hyper <- function(
    n_iter,
    burnin      = 0,
    M,
    mass_init   = 1,          # Initial value for mass
    mass_prior  = c(1, 1),    # Gamma prior for mass: c(shape, rate)
    update_mass = TRUE,       # Whether to update mass parameter
    smoothing,
    hyperparam_base,
    hyperparam_sup,
    known_labels = NULL
) {
  #### 0. Preparations --------------------------------------------------------------------------
  X         <- smoothing$X      # n×T raw data
  basis     <- smoothing$basis
  beta      <- smoothing$beta   # n×L basis coefficients
  tgrid     <- smoothing$time.grid
  n         <- nrow(X)
  T_time    <- ncol(X)
  L         <- basis$nbasis

  # evaluate basis on grid and smooth
  Bmat <- t(eval.basis(tgrid, basis))  # L×T
  X    <- beta %*% Bmat                # n×T

  # which clusters are supervised?
  if (is.null(known_labels)) {
    known_labels <- rep(NA_integer_, n)
  }
  if (length(known_labels) != n)
    stop("`known_labels` must be length n.")
  sup_clusters <- sort(unique(na.omit(known_labels)))

  # helper: pick correct NIW hyperparams for cluster j
  get_hp <- function(j) {
    if (j %in% sup_clusters) {
      hyperparam_sup
    } else {
      hyperparam_base
    }
  }

  #### 1. Hyperpriors for φ (shared) and mass --------------------------------------------------
  c <- hyperparam_base$c
  d <- hyperparam_base$d
  # (we keep φ‐hyperpriors the same for all clusters)

  # Mass parameter initialization
  mass <- mass_init
  mass_shape <- mass_prior[1]
  mass_rate <- mass_prior[2]

  #### 2. Initialize latent variables -----------------------------------------------------------
  # (a) per‐time variances φ_j,t
  phi      <- matrix(invgamma::rinvgamma(M * T_time, shape = c, rate = d),
                     nrow = M, ncol = T_time)

  # (b) NIW for μ‐coefficients
  mu_coef  <- matrix(0, nrow = M, ncol = L)
  for (j in seq_len(M)) {
    hp    <- get_hp(j)
    # reconstruct delta0 = Lambda0 * (nu0−L+1)
    delta0_j <- hp$Lambda0 * (hp$nu0 - L + 1)
    Λ0       <- LaplacesDemon::rinvwishart(hp$nu0, delta0_j)
    m0       <- mvrnorm(1, mu = hp$m0, Sigma = Λ0 / hp$k0)
    mu_coef[j, ] <- mvrnorm(1, mu = m0,  Sigma = Λ0)
  }
  mu <- mu_coef %*% Bmat   # M×T_time mean‐functions

  #### 3. Storage ------------------------------------------------------------------------------
  nsamp       <- n_iter - burnin
  K           <- matrix(0L, nrow = nsamp, ncol = n)
  logL        <- numeric(nsamp)
  COUNTER     <- integer(nsamp)
  nclusters   <- integer(nsamp)
  phi_chain   <- array(  NA_real_, dim = c(nsamp, M, T_time))
  mucoef_chain<- array(  NA_real_, dim = c(nsamp, M, L))
  mu_chain    <- array(  NA_real_, dim = c(nsamp, M, T_time))
  V_chain     <- matrix(  NA_real_, nrow = nsamp, ncol = M)
  p_chain     <- matrix(  NA_real_, nrow = nsamp, ncol = M)
  mass_chain  <- numeric(nsamp)     # Store mass parameter

  # stick‐breaking vars
  V <- numeric(M)
  p <- rep(1/M, M)

  # initialize clusters (clamp supervised)
  K_curr <- sample.int(M, n, replace = TRUE)
  K_curr[!is.na(known_labels)] <- known_labels[!is.na(known_labels)]
  K_new  <- K_curr

  pb <- progressBar(0, max = n_iter, initial = 0, style = "ETA")

  #### 4. Gibbs sampler ------------------------------------------------------------------------
  for (iter in seq_len(n_iter)) {
    ## 4.1 Update φ and μ ----------------------------------------------------------------------
    for (j in seq_len(M)) {
      idx_j <- which(K_curr == j)
      r     <- length(idx_j)

      # (a) φ_j,t
      if (r > 0) {
        c_r <- c + 0.5 * r
        for (t in seq_len(T_time)) {
          d_r      <- d + sum((X[idx_j, t] - mu[j, t])^2)/2
          phi[j, t]<- invgamma::rinvgamma(1, shape = c_r, rate = d_r)
        }
      } else {
        # empty cluster: draw fresh
        phi[j, ] <- invgamma::rinvgamma(T_time, shape = c, rate = d)
      }

      # (b) μ_coef[j,]
      hp    <- get_hp(j)
      # posterior NIW parameters:
      nu_n  <- hp$nu0 + 1
      k_n   <- hp$k0  + 1
      θ_n   <- (hp$k0*hp$m0 + mu_coef[j, ]) / k_n
      Δ_n   <- hp$Lambda0 + (hp$k0/k_n)*tcrossprod(mu_coef[j, ] - hp$m0)
      Λ0    <- LaplacesDemon::rinvwishart(nu_n, Δ_n)
      Λ0_inv<- solve(Λ0)

      # "magic" from the data:
      Mdat  <- Reduce(`+`, lapply(seq_len(T_time),
                                  function(t) tcrossprod(Bmat[, t]) / phi[j, t]))
      post_cov  <- solve(Λ0_inv + r * Mdat)

      # FIX: Handle the dimension issue when r == 1
      if (r > 0) {
        # Ensure beta[idx_j, ] remains a matrix even when r == 1
        beta_j <- beta[idx_j, , drop = FALSE]  # Keep as matrix
        if (r == 1) {
          # For single observation, beta_j is 1×L matrix
          data_term <- Mdat %*% as.vector(beta_j)
        } else {
          # For multiple observations, sum the columns
          data_term <- Mdat %*% colSums(beta_j)
        }
      } else {
        # Empty cluster case
        data_term <- rep(0, L)
      }

      post_mean <- post_cov %*% (
        Λ0_inv %*% (hp$k0 * hp$m0 / k_n) + data_term
      )

      mu_coef[j, ] <- mvrnorm(1, mu = post_mean, Sigma = post_cov)
      mu[j, ]      <- as.vector(mu_coef[j, ] %*% Bmat)
    }

    ## 4.2 Update cluster assignments K_new -----------------------------------------------
    logL_it <- 0
    for (i in seq_len(n)) {
      if (!is.na(known_labels[i])) {
        K_new[i] <- known_labels[i]
      } else {
        logp <- sapply(seq_len(M), function(j) {
          log(p[j]) +
            sum(-0.5*log(2*pi*phi[j,]) - (X[i,] - mu[j,])^2/(2*phi[j,]))
        })
        logp <- logp - max(logp)
        pr   <- exp(logp); pr <- pr/sum(pr)
        K_new[i] <- sample.int(M, 1, prob = pr)
      }
      logL_it <- logL_it +
        sum(-0.5*log(2*pi*phi[K_new[i],]) -
              (X[i,] - mu[K_new[i],])^2/(2*phi[K_new[i],]))
      if (iter > burnin) K[iter-burnin, i] <- K_new[i]
    }
    new_c <- length(unique(K_new)) - length(intersect(unique(K_new), unique(K_curr)))
    K_curr <- K_new

    ## 4.3 Update stick‐breaking weights ---------------------------------------------------
    V[1] <- rbeta(1, 1 + sum(K_curr==1), mass + sum(K_curr>1))
    p[1] <- V[1]
    for (l in 2:(M-1)) {
      V[l] <- rbeta(1,
                    1 + sum(K_curr==l),
                    mass + sum(K_curr>l))
      p[l] <- V[l] * prod(1 - V[1:(l-1)])
    }
    V[M] <- 1; p[M] <- prod(1 - V[1:(M-1)])

    ## 4.4 Update mass parameter using Escobar-West method ---------------------------------
    if (update_mass) {
      # Current number of occupied clusters
      k_curr <- length(unique(K_curr))

      # Escobar-West auxiliary variable update
      # Step 1: Sample auxiliary variable η ~ Beta(mass + 1, n)
      eta <- rbeta(1, mass + 1, n)

      # Step 2: Compute mixture weights for mass update
      # π_η = (mass_shape + k_curr - 1) / (n * (mass_rate - log(eta)))
      # 1 - π_η for the other component
      log_eta <- log(eta)

      # Compute log probabilities to avoid numerical issues
      log_prob1 <- log(mass_shape + k_curr - 1) - log(n) - log(mass_rate - log_eta)
      log_prob2 <- log(mass_shape + k_curr) - log(n) - log(mass_rate - log_eta)

      # Normalize probabilities
      log_probs <- c(log_prob1, log_prob2)
      log_probs <- log_probs - max(log_probs)  # numerical stability
      probs <- exp(log_probs)
      probs <- probs / sum(probs)

      # Step 3: Sample which component for mass posterior
      component <- sample(1:2, 1, prob = probs)

      # Step 4: Update mass parameter
      if (component == 1) {
        # mass ~ Gamma(mass_shape + k_curr - 1, mass_rate - log(η))
        new_shape <- mass_shape + k_curr - 1
        new_rate <- mass_rate - log_eta
      } else {
        # mass ~ Gamma(mass_shape + k_curr, mass_rate - log(η))
        new_shape <- mass_shape + k_curr
        new_rate <- mass_rate - log_eta
      }

      # Ensure positive rate parameter
      if (new_rate > 0) {
        mass <- rgamma(1, shape = new_shape, rate = new_rate)
      } else {
        # Fallback: don't update mass if rate becomes non-positive
        # This can happen with extreme parameter values
        warning("Mass update skipped due to non-positive rate parameter")
      }

      # Ensure mass stays within reasonable bounds
      mass <- pmax(mass, 0.01)  # minimum value
      mass <- pmin(mass, 100)   # maximum value (optional)
    }

    ## 4.5 Store post‐burn‐in ------------------------------------------------------------
    if (iter > burnin) {
      ii             <- iter - burnin
      logL[ii]       <- logL_it
      COUNTER[ii]    <- new_c
      nclusters[ii]  <- length(unique(K_curr))
      phi_chain[ii,,]    <- phi
      mucoef_chain[ii,,] <- mu_coef
      mu_chain[ii,,]     <- mu
      V_chain[ii,]       <- V
      p_chain[ii,]       <- p
      mass_chain[ii]     <- mass
    }

    setTxtProgressBar(pb, iter)
  }
  close(pb)

  #### 5. Return everything ----------------------------------------------------------------------
  list(
    K        = K,
    logL     = logL,
    counter  = COUNTER,
    nclusters= nclusters,
    phi      = phi_chain,
    mu_coef  = mucoef_chain,
    mu       = mu_chain,
    V        = V_chain,
    p        = p_chain,
    mass     = mass_chain,    # Added mass parameter chain
    algo     = list(n_iter = n_iter, burnin = burnin, M = M,
                    mass_init = mass_init, mass_prior = mass_prior,
                    update_mass = update_mass)
  )
}


# library(invgamma)        # for invgamma::invgamma::rinvgamma()
# library(LaplacesDemon)   # for rinvwishart()
# library(MASS)            # for mvrnorm()
# library(pbmcapply)       # for progressBar()
# 
# FBNP_slice <- function(
#     n_iter,
#     burnin      = 0,
#     mass_init   = 1,          # Initial value for mass
#     mass_prior  = c(1, 1),    # Gamma prior for mass: c(shape, rate)
#     update_mass = TRUE,       # Whether to update mass parameter
#     smoothing,
#     hyperparam_base,
#     hyperparam_sup,
#     known_labels = NULL,
#     max_clusters = 50         # Maximum number of clusters to consider (for memory allocation)
# ) {
#   #### 0. Preparations --------------------------------------------------------------------------
#   X         <- smoothing$X      # n×T raw data
#   basis     <- smoothing$basis
#   beta      <- smoothing$beta   # n×L basis coefficients
#   tgrid     <- smoothing$time.grid
#   n         <- nrow(X)
#   T_time    <- ncol(X)
#   L         <- basis$nbasis
#   
#   # evaluate basis on grid and smooth
#   Bmat <- t(eval.basis(tgrid, basis))  # L×T
#   X    <- beta %*% Bmat                # n×T
#   
#   # which clusters are supervised?
#   if (is.null(known_labels)) {
#     known_labels <- rep(NA_integer_, n)
#   }
#   if (length(known_labels) != n)
#     stop("`known_labels` must be length n.")
#   sup_clusters <- sort(unique(na.omit(known_labels)))
#   
#   # helper: pick correct NIW hyperparams for cluster j
#   get_hp <- function(j) {
#     if (j %in% sup_clusters) {
#       hyperparam_sup
#     } else {
#       hyperparam_base
#     }
#   }
#   
#   #### 1. Hyperpriors for φ (shared) and mass --------------------------------------------------
#   c <- hyperparam_base$c
#   d <- hyperparam_base$d
#   
#   # Mass parameter initialization
#   mass <- mass_init
#   mass_shape <- mass_prior[1]
#   mass_rate <- mass_prior[2]
#   
#   #### 2. Initialize latent variables -----------------------------------------------------------
#   # Start with a reasonable number of clusters
#   M_current <- min(n, 10)  # Start with at most 10 clusters
#   
#   # (a) per‐time variances φ_j,t
#   phi      <- matrix(invgamma::rinvgamma(max_clusters * T_time, shape = c, rate = d),
#                      nrow = max_clusters, ncol = T_time)
#   
#   # (b) NIW for μ‐coefficients
#   mu_coef  <- matrix(0, nrow = max_clusters, ncol = L)
#   for (j in seq_len(max_clusters)) {
#     hp    <- get_hp(j)
#     # reconstruct delta0 = Lambda0 * (nu0−L+1)
#     delta0_j <- hp$Lambda0 * (hp$nu0 - L + 1)
#     Λ0       <- LaplacesDemon::rinvwishart(hp$nu0, delta0_j)
#     m0       <- mvrnorm(1, mu = hp$m0, Sigma = Λ0 / hp$k0)
#     mu_coef[j, ] <- mvrnorm(1, mu = m0,  Sigma = Λ0)
#   }
#   mu <- mu_coef %*% Bmat   # max_clusters×T_time mean‐functions
#   
#   #### 3. Storage ------------------------------------------------------------------------------
#   nsamp       <- n_iter - burnin
#   K           <- matrix(0L, nrow = nsamp, ncol = n)
#   logL        <- numeric(nsamp)
#   COUNTER     <- integer(nsamp)
#   nclusters   <- integer(nsamp)
#   phi_chain   <- array(  NA_real_, dim = c(nsamp, max_clusters, T_time))
#   mucoef_chain<- array(  NA_real_, dim = c(nsamp, max_clusters, L))
#   mu_chain    <- array(  NA_real_, dim = c(nsamp, max_clusters, T_time))
#   mass_chain  <- numeric(nsamp)
#   slice_vars  <- matrix(NA_real_, nrow = nsamp, ncol = n)  # Store slice variables
#   
#   # Initialize cluster assignments (clamp supervised)
#   K_curr <- sample.int(M_current, n, replace = TRUE)
#   K_curr[!is.na(known_labels)] <- known_labels[!is.na(known_labels)]
#   K_new  <- K_curr
#   
#   # Initialize slice variables
#   u <- runif(n)
#   
#   pb <- progressBar(0, max = n_iter, initial = 0, style = "ETA")
#   
#   #### 4. Helper functions for slice sampler ---------------------------------------------------
#   
#   # Store stick-breaking variables persistently
#   V_stick <- numeric(max_clusters)
#   p_weights <- numeric(max_clusters)
#   
#   # Function to update stick-breaking weights based on cluster assignments
#   update_stick_weights <- function(K_curr, mass, M_active) {
#     # Count observations in each cluster
#     cluster_counts <- table(factor(K_curr, levels = 1:M_active))
#     n_after <- numeric(M_active)
#     
#     # Calculate n_j^{>} = number of observations in clusters > j
#     for (j in 1:M_active) {
#       n_after[j] <- sum(cluster_counts[(j+1):M_active])
#     }
#     
#     # Update V_j ~ Beta(1 + n_j, mass + n_j^{>})
#     for (j in 1:(M_active-1)) {
#       V_stick[j] <- rbeta(1, 1 + cluster_counts[j], mass + n_after[j])
#     }
#     V_stick[M_active] <- 1  # Last component
#     
#     # Compute weights p_j = V_j * prod(1 - V_k) for k < j
#     p_weights[1] <- V_stick[1]
#     if (M_active > 1) {
#       for (j in 2:M_active) {
#         p_weights[j] <- V_stick[j] * prod(1 - V_stick[1:(j-1)])
#       }
#     }
#     
#     return(p_weights[1:M_active])
#   }
#   
#   # Function to find required number of clusters given slice variables
#   find_required_clusters <- function(u, mass) {
#     # We need enough clusters so that sum(p_1, ..., p_K) >= max(u)
#     # Start with current maximum cluster index
#     min_clusters <- max(K_curr, na.rm = TRUE)
#     
#     # Ensure we have at least as many clusters as currently occupied
#     return(max(min_clusters, ceiling(log(max(u)) / log(mass/(mass+1))) + 3))
#   }
#   
#   #### 5. Gibbs sampler ------------------------------------------------------------------------
#   for (iter in seq_len(n_iter)) {
#     
#     ## 5.1 Determine active clusters and update stick-breaking weights -------------------------
#     # Start with at least the maximum currently used cluster
#     M_active <- max(max(K_curr), 2)  # At least 2 clusters
#     M_active <- min(M_active + 2, max_clusters)  # Add buffer and respect limit
#     
#     # Update stick-breaking weights
#     p_current <- update_stick_weights(K_curr, mass, M_active)
#     
#     ## 5.2 Update slice variables u_i ----------------------------------------------------------
#     for (i in seq_len(n)) {
#       # Slice variable must be below the probability of current cluster assignment
#       cluster_idx <- K_curr[i]
#       if (cluster_idx <= length(p_current) && p_current[cluster_idx] > 0) {
#         u[i] <- runif(1, 0, p_current[cluster_idx])
#       } else {
#         u[i] <- runif(1, 0, 0.001)  # Small positive value for safety
#       }
#     }
#     
#     ## 5.3 Expand clusters if needed based on slice variables ----------------------------------
#     # Check if we need more clusters to satisfy slice constraints
#     while (M_active < max_clusters && max(u) > sum(p_current)) {
#       M_active <- M_active + 1
#       p_current <- update_stick_weights(K_curr, mass, M_active)
#     }
#     
#     ## 5.3 Generate new parameters for newly created clusters ------------------------------------
#     if (M_active > M_current) {
#       for (j in (M_current + 1):M_active) {
#         # Generate new phi values
#         phi[j, ] <- invgamma::rinvgamma(T_time, shape = c, rate = d)
#         
#         # Generate new mu_coef values
#         hp <- get_hp(j)
#         delta0_j <- hp$Lambda0 * (hp$nu0 - L + 1)
#         Λ0 <- LaplacesDemon::rinvwishart(hp$nu0, delta0_j)
#         m0 <- mvrnorm(1, mu = hp$m0, Sigma = Λ0 / hp$k0)
#         mu_coef[j, ] <- mvrnorm(1, mu = m0, Sigma = Λ0)
#         mu[j, ] <- as.vector(mu_coef[j, ] %*% Bmat)
#       }
#     }
#     M_current <- M_active
#     
#     ## 5.3 Update φ and μ for active clusters --------------------------------------------------
#     for (j in seq_len(M_active)) {
#       idx_j <- which(K_curr == j)
#       r     <- length(idx_j)
#       
#       # (a) φ_j,t
#       if (r > 0) {
#         c_r <- c + 0.5 * r
#         for (t in seq_len(T_time)) {
#           d_r      <- d + sum((X[idx_j, t] - mu[j, t])^2)/2
#           phi[j, t]<- invgamma::rinvgamma(1, shape = c_r, rate = d_r)
#         }
#       } else {
#         # empty cluster: draw fresh
#         phi[j, ] <- invgamma::rinvgamma(T_time, shape = c, rate = d)
#       }
#       
#       # (b) μ_coef[j,]
#       hp    <- get_hp(j)
#       # posterior NIW parameters:
#       nu_n  <- hp$nu0 + 1
#       k_n   <- hp$k0  + 1
#       θ_n   <- (hp$k0*hp$m0 + mu_coef[j, ]) / k_n
#       Δ_n   <- hp$Lambda0 + (hp$k0/k_n)*tcrossprod(mu_coef[j, ] - hp$m0)
#       Λ0    <- LaplacesDemon::rinvwishart(nu_n, Δ_n)
#       Λ0_inv<- solve(Λ0)
#       
#       # "magic" from the data:
#       Mdat  <- Reduce(`+`, lapply(seq_len(T_time),
#                                   function(t) tcrossprod(Bmat[, t]) / phi[j, t]))
#       post_cov  <- solve(Λ0_inv + r * Mdat)
#       
#       # Handle the dimension issue when r == 1
#       if (r > 0) {
#         beta_j <- beta[idx_j, , drop = FALSE]
#         if (r == 1) {
#           data_term <- Mdat %*% as.vector(beta_j)
#         } else {
#           data_term <- Mdat %*% colSums(beta_j)
#         }
#       } else {
#         data_term <- rep(0, L)
#       }
#       
#       post_mean <- post_cov %*% (
#         Λ0_inv %*% (hp$k0 * hp$m0 / k_n) + data_term
#       )
#       
#       mu_coef[j, ] <- mvrnorm(1, mu = post_mean, Sigma = post_cov)
#       mu[j, ]      <- as.vector(mu_coef[j, ] %*% Bmat)
#     }
#     
#     ## 5.5 Update cluster assignments K_new with slice constraints -----------------------------
#     logL_it <- 0
#     
#     for (i in seq_len(n)) {
#       if (!is.na(known_labels[i])) {
#         K_new[i] <- known_labels[i]
#       } else {
#         # Only consider clusters j where p_j >= u_i (slice constraint)
#         valid_clusters <- which(p_current >= u[i] & p_current > 0)
#         
#         if (length(valid_clusters) == 0) {
#           # Fallback: assign to cluster with highest probability
#           valid_clusters <- which.max(p_current)
#         }
#         
#         # Compute log probabilities for valid clusters only
#         logp <- numeric(length(valid_clusters))
#         for (idx in seq_along(valid_clusters)) {
#           j <- valid_clusters[idx]
#           logp[idx] <- log(p_current[j]) +
#             sum(-0.5*log(2*pi*phi[j,]) - (X[i,] - mu[j,])^2/(2*phi[j,]))
#         }
#         
#         # Normalize and sample
#         logp <- logp - max(logp)
#         pr   <- exp(logp)
#         pr   <- pr / sum(pr)
#         
#         # Sample from valid clusters
#         if (length(valid_clusters) == 1) {
#           K_new[i] <- valid_clusters[1]
#         } else {
#           K_new[i] <- sample(valid_clusters, 1, prob = pr)
#         }
#       }
#       
#       # Calculate log-likelihood contribution
#       cluster_i <- K_new[i]
#       logL_it <- logL_it +
#         sum(-0.5*log(2*pi*phi[cluster_i,]) -
#               (X[i,] - mu[cluster_i,])^2/(2*phi[cluster_i,]))
#       
#       if (iter > burnin) K[iter-burnin, i] <- K_new[i]
#     }
#     
#     new_c <- length(unique(K_new)) - length(intersect(unique(K_new), unique(K_curr)))
#     K_curr <- K_new
#     
#     ## 5.6 Update mass parameter using Escobar-West method -------------------------------------
#     if (update_mass) {
#       k_curr <- length(unique(K_curr))
#       
#       # Escobar-West auxiliary variable update
#       eta <- rbeta(1, mass + 1, n)
#       log_eta <- log(eta)
#       
#       # Compute log probabilities
#       log_prob1 <- log(mass_shape + k_curr - 1) - log(n) - log(mass_rate - log_eta)
#       log_prob2 <- log(mass_shape + k_curr) - log(n) - log(mass_rate - log_eta)
#       
#       # Normalize probabilities
#       log_probs <- c(log_prob1, log_prob2)
#       log_probs <- log_probs - max(log_probs)
#       probs <- exp(log_probs)
#       probs <- probs / sum(probs)
#       
#       # Sample component
#       component <- sample(1:2, 1, prob = probs)
#       
#       # Update mass parameter
#       if (component == 1) {
#         new_shape <- mass_shape + k_curr - 1
#         new_rate <- mass_rate - log_eta
#       } else {
#         new_shape <- mass_shape + k_curr
#         new_rate <- mass_rate - log_eta
#       }
#       
#       if (new_rate > 0) {
#         mass <- rgamma(1, shape = new_shape, rate = new_rate)
#       } else {
#         warning("Mass update skipped due to non-positive rate parameter")
#       }
#       
#       # Ensure mass stays within reasonable bounds
#       mass <- pmax(mass, 0.01)
#       mass <- pmin(mass, 100)
#     }
#     
#     ## 5.7 Store post‐burn‐in -----------------------------------------------------------------
#     if (iter > burnin) {
#       ii             <- iter - burnin
#       logL[ii]       <- logL_it
#       COUNTER[ii]    <- new_c
#       nclusters[ii]  <- length(unique(K_curr))
#       phi_chain[ii,,]    <- phi
#       mucoef_chain[ii,,] <- mu_coef
#       mu_chain[ii,,]     <- mu
#       mass_chain[ii]     <- mass
#       slice_vars[ii, ]   <- u
#     }
#     
#     setTxtProgressBar(pb, iter)
#   }
#   close(pb)
#   
#   #### 6. Return everything -------------------------------------------------------------------
#   list(
#     K        = K,
#     logL     = logL,
#     counter  = COUNTER,
#     nclusters= nclusters,
#     phi      = phi_chain,
#     mu_coef  = mucoef_chain,
#     mu       = mu_chain,
#     mass     = mass_chain,
#     slice_vars = slice_vars,  # Added slice variables
#     algo     = list(n_iter = n_iter, burnin = burnin, max_clusters = max_clusters,
#                     mass_init = mass_init, mass_prior = mass_prior, 
#                     update_mass = update_mass, sampler = "walker_slice")
#   )
# }

