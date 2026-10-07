# ABC_fit.R
suppressPackageStartupMessages({
  library(zellkonverter)
  library(Matrix)
  library(SingleCellExperiment)
  library(MASS)
})

options(error = function() {
  traceback(2)
  quit(save = "no", status = 1)
})

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

progress_log <- "logs/ABC_fit.log"
if (file.exists(progress_log)) invisible(file.remove(progress_log))

log_text <- function(txt) {
  cat(sprintf("%s\n", txt), file = progress_log, append = TRUE)
}

log_iter <- function(stage, iter, accepted = NA, extra = "") {
  msg <- sprintf("[%s] iteration=%d", stage, iter)
  if (!is.na(accepted)) {
    msg <- paste0(msg, sprintf(" accepted=%d", accepted))
  }
  if (nzchar(extra)) {
    msg <- paste0(msg, " ", extra)
  }
  log_text(msg)
}

say <- function(txt) {
  cat(txt, "\n")
  log_text(txt)
}

say("ABC_fit.R started")

# Load data
norman <- readH5AD(
  "~/projects/scRNA_ABC/data/Norman_2019.h5ad",
  reader = "R",
  use_hdf5 = FALSE
)

perturb <- colData(norman)$perturbation_name
sce_control <- norman[, perturb == "control"]

X_real <- assay(sce_control, "counts")

cells_use <- seq_len(ncol(X_real))
genes_use <- seq_len(nrow(X_real))

X_obs <- t(X_real[genes_use, cells_use])
N_use <- nrow(X_obs)
J_use <- ncol(X_obs)

# K = number of unique single-gene perturbations
pert_all <- as.character(colData(norman)$perturbation_name)
pert_all <- pert_all[!is.na(pert_all) & pert_all != "control"]
pert_single <- pert_all[!grepl("\\+", pert_all)]
K_use <- length(unique(pert_single))

say(sprintf(
  "Loaded X_obs: N=%d, J=%d, K=%d",
  N_use, J_use, K_use
))

# 1) Summary statistics for ABC
summaries <- function(X) {

  x_all <- as.vector(X)
  x_nonzero <- x_all[x_all > 0 & !is.na(x_all)]

  if (length(x_nonzero) == 0 || any(!is.finite(x_all))) {
    return(rep(NA, 54))
  }

  C_i <- Matrix::rowSums(X) / ncol(X)
  D_i <- Matrix::rowSums(X > 0) / ncol(X)
  G_j <- Matrix::colSums(X) / nrow(X)
  M_j <- Matrix::colSums(X > 0) / nrow(X)

  probs_use <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  c(
    quantile(x_all, probs = probs_use, na.rm = TRUE, names = FALSE),
    quantile(x_nonzero, probs = probs_use, na.rm = TRUE, names = FALSE),
    quantile(as.numeric(C_i), probs = probs_use, na.rm = TRUE, names = FALSE),
    quantile(as.numeric(D_i), probs = probs_use, na.rm = TRUE, names = FALSE),
    quantile(as.numeric(G_j), probs = probs_use, na.rm = TRUE, names = FALSE),
    quantile(as.numeric(M_j), probs = probs_use, na.rm = TRUE, names = FALSE)
  )
}

# Summary of X_obs
say("About to compute S_obs")
S_obs <- summaries(X_obs)
say("Finished computing S_obs")
saveRDS(X_obs, file = "results/ABC_X_obs.rds", compress = TRUE)
say("Saved results/ABC_X_obs.rds")
rm(X_obs, X_real, sce_control, norman, perturb, pert_all, pert_single)

# 2) Approximate Bayesian Computation distance
distance_abc <- function(S_sim, S_obs) {
  d <- sqrt(sum((S_sim - S_obs)^2))
  
  # If distance is NA, return Infinity to force rejection
  if (!is.finite(d)) {
    return(Inf)
  }
  
  d
}

# 3) Prior ranges
prior_range <- list(
  phi = c(0.05, 0.95),
  gamma = c(1.5, 25),
  rho = c(10, 60),
  alpha = c(2, 20),
  sigma2 = c(0.05, 0.15),
  eta = c(0.1, 5),
  zeta = c(0.1, 10),
  chi = c(0.1, 5),
  omega = c(0.1, 10)
)

# Sample lambda = (phi, gamma, rho, alpha, sigma2, xi, kappa)
sample_lambda <- function() {
  list(
    # Uniform prior
    phi = runif(1, min = prior_range$phi[1], max = prior_range$phi[2]),
    
    # Log-uniform priors
    gamma = exp(runif(1, min = log(prior_range$gamma[1]), max = log(prior_range$gamma[2]))),
    rho = exp(runif(1, min = log(prior_range$rho[1]), max = log(prior_range$rho[2]))),
    alpha = exp(runif(1, min = log(prior_range$alpha[1]), max = log(prior_range$alpha[2]))),
    sigma2 = exp(runif(1, min = log(prior_range$sigma2[1]), max = log(prior_range$sigma2[2]))),
    eta = exp(runif(1, min = log(prior_range$eta[1]), max = log(prior_range$eta[2]))),
    zeta = exp(runif(1, min = log(prior_range$zeta[1]), max = log(prior_range$zeta[2]))),
    chi = exp(runif(1, min = log(prior_range$chi[1]), max = log(prior_range$chi[2]))),
    omega = exp(runif(1, min = log(prior_range$omega[1]),max = log(prior_range$omega[2])))
  )
}

# 4) Prior density on the transformed parameter scale
prior_density <- function(lambda) {
  
  lambda_values <- unlist(lambda, use.names = FALSE)
  
  if (any(!is.finite(lambda_values))) {
    return(0)
  }
  
  inside_support <- (
    lambda$phi >= prior_range$phi[1] &&
      lambda$phi <= prior_range$phi[2] &&
      lambda$gamma >= prior_range$gamma[1] &&
      lambda$gamma <= prior_range$gamma[2] &&
      lambda$rho >= prior_range$rho[1] &&
      lambda$rho <= prior_range$rho[2] &&
      lambda$alpha >= prior_range$alpha[1] &&
      lambda$alpha <= prior_range$alpha[2] &&
      lambda$sigma2 >= prior_range$sigma2[1] &&
      lambda$sigma2 <= prior_range$sigma2[2] &&
      lambda$eta >= prior_range$eta[1] &&
      lambda$eta <= prior_range$eta[2] &&
      lambda$zeta >= prior_range$zeta[1] &&
      lambda$zeta <= prior_range$zeta[2] &&
      lambda$chi >= prior_range$chi[1] &&
      lambda$chi <= prior_range$chi[2] &&
      lambda$omega >= prior_range$omega[1] &&
      lambda$omega <= prior_range$omega[2]
  )
  
  if (!inside_support) {
    return(0)
  }
  
  d_phi <- 1 / (prior_range$phi[2] - prior_range$phi[1])
  d_log_gamma <- 1 / (log(prior_range$gamma[2]) - log(prior_range$gamma[1]))
  d_log_rho <- 1 / (log(prior_range$rho[2]) - log(prior_range$rho[1]))
  d_log_alpha <- 1 / (log(prior_range$alpha[2]) - log(prior_range$alpha[1]))
  d_log_sigma2 <- 1 / (log(prior_range$sigma2[2]) - log(prior_range$sigma2[1]))
  d_log_eta <- 1 / (log(prior_range$eta[2]) - log(prior_range$eta[1]))
  d_log_zeta <- 1 / (log(prior_range$zeta[2]) - log(prior_range$zeta[1]))
  d_log_chi <- 1 / (log(prior_range$chi[2]) - log(prior_range$chi[1]))
  d_log_omega <- 1 / (log(prior_range$omega[2]) - log(prior_range$omega[1]))
  
  d_phi * d_log_gamma * d_log_rho * d_log_alpha * d_log_sigma2 * d_log_eta * d_log_zeta *
    d_log_chi * d_log_omega
}

# check hyperparameter vector inside support of prior
in_prior_support <- function(theta) {
  prior_density(theta) > 0
}

# 5) Simulator
simulate_X_from_lambda <- function(theta, N, J, K) {
  
  # hyperparameters
  phi <- theta$phi
  gamma <- theta$gamma
  rho <- theta$rho
  alpha <- theta$alpha
  sigma2 <- theta$sigma2
  eta <- theta$eta
  zeta <- theta$zeta
  chi <- theta$chi
  omega <- theta$omega
  
  # (a) TF activities: a_ik ~ Gamma(alpha, alpha)
  a_ik <- matrix(rgamma(N * K, shape = alpha, rate = alpha), nrow = N, ncol = K)
  
  # (b) TF -> gene effects: beta_kj ~ Normal(0, sigma2)
  beta_kj <- matrix(rnorm(K * J, mean = 0, sd = sqrt(sigma2)), nrow = K, ncol = J)
  
  # (c) delta_kj ~ Bernoulli(psi_j * psi_k)
  psi_j <- rbeta(J, shape1 = eta, shape2 = zeta)
  psi_k <- rbeta(K, shape1 = chi, shape2 = omega)
  prob_kj <- outer(psi_k, psi_j, "*")
  delta_kj <- matrix(rbinom(K * J, size = 1, prob = as.vector(prob_kj)), 
                     nrow = K, ncol = J)
  beta_kj <- beta_kj * delta_kj
  
  # (d) Mean expression: m_ij = exp(sum_k a_ik * beta_kj)
  m_ij <- exp(a_ik %*% beta_kj)
  
  # (e) Gene-specific on-probabilities
  shape1 <- phi
  shape2 <- (1 - phi)
  pi_j <- rbeta(J, shape1 = shape1, shape2 = shape2)
  
  # On/off states: v_ij ~ Bernoulli(pi_j)
  v_ij <- matrix(
    rbinom(N * J, size = 1, prob = rep(pi_j, each = N)),
    nrow = N, ncol = J
  )
  
  # (f) Biological expression
  z_ij <- matrix(0, nrow = N, ncol = J)
  idx_on <- which(v_ij == 1, arr.ind = TRUE)
  
  if (nrow(idx_on) > 0) {
    shape_vec <- m_ij[idx_on] * gamma
    z_ij[idx_on] <- rgamma(nrow(idx_on), shape = shape_vec, rate = gamma)
  }
  
  # (g) Technical multiplicative noise: t_ij ~ Gamma(rho, rho)
  t_ij <- matrix(rgamma(N * J, shape = rho, rate = rho), nrow = N, ncol = J)
  
  # (h) Observation: X_ij | z_ij, t_ij ~ Poisson(t_ij * z_ij)
  lambda_ij <- t_ij * z_ij
  X_ij <- matrix(rpois(N * J, lambda = as.vector(lambda_ij)), nrow = N, ncol = J)
  
  list(
    X = X_ij
  )
}

# Convert lambda to numeric vector
lambda_to_vec <- function(lambda) {
  c(
    phi = lambda$phi,
    log_gamma = log(lambda$gamma),
    log_rho = log(lambda$rho),
    log_alpha = log(lambda$alpha),
    log_sigma2 = log(lambda$sigma2),
    log_eta = log(lambda$eta),
    log_zeta = log(lambda$zeta),
    log_chi = log(lambda$chi),
    log_omega = log(lambda$omega)
  )
}

# Convert back
vec_to_lambda <- function(v) {
  list(
    phi = unname(v["phi"]),
    gamma = exp(unname(v["log_gamma"])),
    rho = exp(unname(v["log_rho"])),
    alpha = exp(unname(v["log_alpha"])),
    sigma2 = exp(unname(v["log_sigma2"])),
    eta = exp(unname(v["log_eta"])),
    zeta = exp(unname(v["log_zeta"])),
    chi = exp(unname(v["log_chi"])),
    omega = exp(unname(v["log_omega"]))
  )
}

# Compute the weighted covariance matrix of the particles
weighted_covariance <- function(X, w) {
  mu <- colSums(X * w)
  Xc <- sweep(X, 2, mu, "-")
  crossprod(Xc * sqrt(w), Xc * sqrt(w))
}

# Stabilize covariance matrix so chol() does not fail
make_positive_definite <- function(Sigma, min_eig = 1e-8) {
  
  Sigma <- as.matrix(Sigma)
  Sigma[!is.finite(Sigma)] <- 0
  
  # Force symmetry
  Sigma <- (Sigma + t(Sigma)) / 2
  
  eig_vals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  min_val <- min(eig_vals, na.rm = TRUE)
  
  if (!is.finite(min_val)) {
    return(diag(min_eig, nrow(Sigma)))
  }
  
  if (min_val < min_eig) {
    diag(Sigma) <- diag(Sigma) + (min_eig - min_val)
  }
  
  # Force symmetry again after diagonal adjustment
  Sigma <- (Sigma + t(Sigma)) / 2
  
  Sigma
}

# Compute the log density of a multivariate Gaussian distribution
mvnorm_logpdf <- function(x, mean, Sigma) {
  
  Sigma <- make_positive_definite(Sigma, min_eig = 1e-8)
  
  p <- length(x)
  
  cholS <- tryCatch(
    chol(Sigma),
    error = function(e) {
      chol(make_positive_definite(Sigma, min_eig = 1e-6))
    }
  )
  
  z <- backsolve(cholS, x - mean, transpose = TRUE)
  
  out <- -0.5 * (
    p * log(2 * pi) +
      2 * sum(log(diag(cholS))) +
      sum(z^2)
  )
  
  if (!is.finite(out)) {
    return(-Inf)
  }
  
  out
}

# Take one old particle and perturb it with a Gaussian random walk
perturb_lambda <- function(lambda, Sigma_kernel) {
  
  Sigma_kernel <- make_positive_definite(Sigma_kernel, min_eig = 1e-8)
  
  v <- lambda_to_vec(lambda)
  v_new <- as.numeric(MASS::mvrnorm(1, mu = v, Sigma = Sigma_kernel))
  names(v_new) <- names(v)
  
  vec_to_lambda(v_new)
}

# Gaussian proposal density from the old particle to the new one
kernel_density <- function(lambda_new, lambda_old, Sigma_kernel) {
  x <- lambda_to_vec(lambda_new)
  m <- lambda_to_vec(lambda_old)
  exp(mvnorm_logpdf(x, m, Sigma_kernel))
}

# Diagnostics
effective_sample_size <- function(w) {
  1 / sum(w^2)
}

# ABC_SMC procedure
abc_smc <- function(S_obs, N, J, K, W, T_pop, alpha_quantile, pilot_size) {
  
  populations <- vector("list", T_pop)
  eps_vec <- numeric(T_pop)
  
  # Population 1: Initial Rejection Sampling
  cat("Population 1\n")
  log_text("Starting population 1 pilot")
  
  # Pilot simulations to determine epsilon_1
  d_store <- numeric(pilot_size)
  
  for (m in seq_len(pilot_size)) {
    
    if (m %% 10 == 0 || m == 1 || m == pilot_size) {
      log_iter("pilot_pop1", m)
    }
    
    th <- sample_lambda()
    out_sim <- simulate_X_from_lambda(th, N, J, K)
    d_now <- distance_abc(summaries(out_sim$X), S_obs)
    
    d_store[m] <- d_now
  }
  
  finite_d_store <- d_store[is.finite(d_store)]
  
  if (length(finite_d_store) == 0L) {
    stop("All pilot simulations produced non-finite ABC distances.")
  }
  
  eps_vec[1] <- as.numeric(
    quantile(
      finite_d_store,
      probs = alpha_quantile,
      names = FALSE
    )
  )
  
  cat("epsilon_1 =", eps_vec[1], "\n")
  log_text(sprintf("epsilon_1 = %f", eps_vec[1]))
  
  # Generate population 1
  particles <- vector("list", W)
  distances <- numeric(W)
  
  n_accept <- 0
  total_prop_1 <- 0
  
  log_text("Starting population 1 rejection ABC")
  
  while (n_accept < W) {
    
    th <- sample_lambda()
    out_sim <- simulate_X_from_lambda(th, N, J, K)
    d_now <- distance_abc(summaries(out_sim$X), S_obs)
    
    total_prop_1 <- total_prop_1 + 1
    
    if (total_prop_1 %% 10 == 0 || total_prop_1 == 1) {
      log_iter(
        "pop1_proposal",
        total_prop_1,
        n_accept
      )
    }
    
    if (is.finite(d_now) && d_now <= eps_vec[1]) {
      n_accept <- n_accept + 1
      particles[[n_accept]] <- th
      distances[n_accept] <- d_now
      
      log_iter(
        "pop1_accept",
        total_prop_1,
        n_accept,
        sprintf("distance=%f", d_now)
      )
    }
  }
  
  weights <- rep(1 / W, W)
  
  populations[[1]] <- list(
    particles = particles,
    weights = weights,
    distances = distances,
    distances_for_epsilon = distances,
    epsilon = eps_vec[1],
    proposals = total_prop_1,
    ESS = effective_sample_size(weights)
  )
  
  cat(
    "accepted =", W,
    " proposals =", total_prop_1,
    " ESS =", effective_sample_size(weights),
    "\n"
  )
  
  log_text(sprintf(
    "Population 1 finished: accepted=%d proposals=%d ESS=%f",
    W,
    total_prop_1,
    effective_sample_size(weights)
  ))
  
  # Populations 2 to T_pop
  if (T_pop >= 2) {
    
    for (t in 2:T_pop) {
      
      cat("Population", t, "\n")
      log_text(sprintf("Starting population %d", t))
      
      prev_particles <- populations[[t - 1]]$particles
      prev_weights <- populations[[t - 1]]$weights
      prev_distances <- populations[[t - 1]]$distances_for_epsilon
      
      eps_vec[t] <- as.numeric(
        quantile(
          prev_distances,
          probs = alpha_quantile,
          na.rm = TRUE
        )
      )
      
      cat(paste0("epsilon_", t, " = "), eps_vec[t], "\n")
      log_text(sprintf("epsilon_%d = %f", t, eps_vec[t]))
      
      # Construct the Gaussian perturbation kernel
      prev_mat <- do.call(
        rbind,
        lapply(prev_particles, lambda_to_vec)
      )
      
      Sigma_kernel <- 2 * weighted_covariance(
        prev_mat,
        prev_weights
      )
      
      Sigma_kernel <- make_positive_definite(
        Sigma_kernel,
        min_eig = 1e-8
      )
      
      particles_t <- vector("list", W)
      distances_t <- numeric(W)
      weights_t <- numeric(W)
      
      n_accept <- 0
      total_prop_t <- 0
      
      while (n_accept < W) {
        
        # Select and perturb a previous particle
        idx <- sample(
          seq_len(W),
          size = 1,
          prob = prev_weights
        )
        
        th_star <- perturb_lambda(
          prev_particles[[idx]],
          Sigma_kernel
        )
        
        total_prop_t <- total_prop_t + 1
        
        if (total_prop_t %% 10 == 0 || total_prop_t == 1) {
          log_iter(
            sprintf("pop%d_proposal", t),
            total_prop_t,
            n_accept
          )
        }
        
        if (!in_prior_support(th_star)) {
          next
        }
        
        # Simulate and calculate ABC distance
        out_sim <- simulate_X_from_lambda(
          th_star,
          N,
          J,
          K
        )
        
        d_now <- distance_abc(
          summaries(out_sim$X),
          S_obs
        )
        
        if (is.finite(d_now) && d_now <= eps_vec[t]) {
          
          num <- prior_density(th_star)
          
          denom <- sum(
            prev_weights *
              sapply(prev_particles, function(old) {
                kernel_density(
                  th_star,
                  old,
                  Sigma_kernel
                )
              })
          )
          
          if (
            !is.finite(num) ||
            num <= 0 ||
            !is.finite(denom) ||
            denom <= 0
          ) {
            next
          }
          
          n_accept <- n_accept + 1
          particles_t[[n_accept]] <- th_star
          distances_t[n_accept] <- d_now
          weights_t[n_accept] <- num / denom
          
          log_iter(
            sprintf("pop%d_accept", t),
            total_prop_t,
            n_accept,
            sprintf("distance=%f", d_now)
          )
        }
      }
      
      # Normalize importance weights
      weights_t <- weights_t / sum(weights_t)
      
      ess_t <- effective_sample_size(weights_t)
      
      log_text(sprintf(
        "Population %d pre-resample ESS=%f",
        t,
        ess_t
      ))
      
      # Resample particles and corresponding distances
      accepted_distances_t <- distances_t
      
      if (ess_t < W / 2) {
        
        log_text(sprintf(
          "Population %d resampling triggered",
          t
        ))
        
        idx_resample <- sample(
          seq_len(W),
          size = W,
          replace = TRUE,
          prob = weights_t
        )
        
        particles_t <- particles_t[idx_resample]
        distances_t <- distances_t[idx_resample]
        weights_t <- rep(1 / W, W)
        
        ess_t <- effective_sample_size(weights_t)
      }
      
      populations[[t]] <- list(
        particles = particles_t,
        weights = weights_t,
        distances = distances_t,
        distances_for_epsilon = accepted_distances_t,
        epsilon = eps_vec[t],
        proposals = total_prop_t,
        Sigma_kernel = Sigma_kernel,
        ESS = ess_t
      )
      
      cat(
        "accepted =", W,
        " proposals =", total_prop_t,
        " ESS =", ess_t,
        "\n"
      )
      
      log_text(sprintf(
        "Population %d finished: accepted=%d proposals=%d ESS=%f",
        t,
        W,
        total_prop_t,
        ess_t
      ))
    }
  }
  
  list(
    populations = populations,
    eps_vec = eps_vec
  )
}

# Target accepted samples per population
W <- 200

# Number of ABC-SMC populations
T_pop <- 5

# Quantile used to update tolerances
alpha_quantile <- 0.25

# Pilot size
pilot_size <- 2000

say(sprintf(
  "Running ABC-SMC with W=%d, T_pop=%d, alpha_quantile=%0.2f, pilot_size=%d",
  W, T_pop, alpha_quantile, pilot_size
))

# run ABC-SMC
smc_fit <- abc_smc(
  S_obs = S_obs,
  N = N_use,
  J = J_use,
  K = K_use,
  W = W,
  T_pop = T_pop,
  alpha_quantile = alpha_quantile,
  pilot_size = pilot_size
)

final_population <- smc_fit$populations[[T_pop]]
saveRDS(
  list(
    N_use = N_use,
    J_use = J_use,
    K_use = K_use,
    final_population = final_population
  ),
  file = "results/ABC_fit.rds"
)

say("Saved results/ABC_fit.rds")
say("ABC_fit.R finished")
