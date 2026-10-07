# ABC_sim.R

suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
})

options(error = function() {
  traceback(2)
  quit(save = "no", status = 1)
})

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

progress_log <- "logs/ABC_sim.log"
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

say("ABC_sim.R started")

# 1) Load ABC fit and observed data
fit_obj <- readRDS("results/ABC_fit.rds")
X_obs <- readRDS("results/ABC_X_obs.rds")

N_use <- fit_obj$N_use
J_use <- fit_obj$J_use
K_use <- fit_obj$K_use
last_pop <- fit_obj$final_population

say("Loaded results/ABC_fit.rds")
say("Loaded results/ABC_X_obs.rds")
say(sprintf("N_use = %d", N_use))
say(sprintf("J_use = %d", J_use))
say(sprintf("K_use = %d", K_use))

if (nrow(X_obs) != N_use || ncol(X_obs) != J_use) {
  stop("Dimensions of ABC_X_obs.rds do not match N_use and J_use in ABC_fit.rds.")
}

eps_t <- last_pop$epsilon
acc_t <- length(last_pop$particles)
prop_t <- if (!is.null(last_pop$proposals)) last_pop$proposals else NA
ess_t <- if (!is.null(last_pop$ESS)) last_pop$ESS else NA

say(sprintf(
  "Final population epsilon = %.6f accepted = %d proposals = %s ESS = %s",
  eps_t,
  acc_t,
  ifelse(is.na(prop_t), "NA", as.character(prop_t)),
  ifelse(is.na(ess_t), "NA", format(ess_t, digits = 8))
))

theta_acc <- last_pop$particles
w_acc <- last_pop$weights
d_acc <- last_pop$distances

if (length(theta_acc) == 0L) {
  stop("No posterior particles were found in final_population.")
}

if (length(w_acc) != length(theta_acc)) {
  stop("Length of final_population$weights does not match number of particles.")
}

if (any(!is.finite(w_acc)) || any(w_acc < 0) || sum(w_acc) <= 0) {
  stop("Invalid posterior weights in final_population.")
}

w_acc <- w_acc / sum(w_acc)

# 2) Prior distribution
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

# 3) Simulator
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

# 4) Prior samples
prior_size <- 2000

say(sprintf("Prior size = %d", prior_size))
say(sprintf("Posterior size = %d", length(theta_acc)))

prior_draws <- vector("list", prior_size)

say(sprintf(
  "Drawing prior samples: prior_size = %d",
  prior_size
))

for (m in seq_len(prior_size)) {

  if (m %% 50 == 0 || m == 1 || m == prior_size) {
    log_iter("prior_draw", m)
    cat(sprintf(
      "prior_draw %d / %d\n",
      m,
      prior_size
    ))
  }

  th <- sample_lambda()

  prior_draws[[m]] <- c(
    phi = th$phi,
    gamma = th$gamma,
    rho = th$rho,
    alpha = th$alpha,
    sigma2 = th$sigma2,
    eta = th$eta,
    zeta = th$zeta,
    chi = th$chi,
    omega = th$omega
  )
}

prior_df <- as.data.frame(
  do.call(rbind, prior_draws)
)

rm(prior_draws)
gc()

# 5) Posterior particles for PPC
R_post <- length(theta_acc)

uniform_weights <- max(
  abs(w_acc - rep(1 / R_post, R_post))
) < 1e-12

if (uniform_weights) {

  theta_pp <- theta_acc

  say(
    "Final posterior weights are uniform; using final particles directly for PPC."
  )

} else {

  ppc_idx <- sample(
    seq_len(R_post),
    size = R_post,
    replace = TRUE,
    prob = w_acc
  )

  theta_pp <- theta_acc[ppc_idx]

  say(
    "Final posterior weights are non-uniform; resampled particles using ABC-SMC weights for PPC."
  )
}

post_df <- do.call(
  rbind,
  lapply(theta_pp, function(th) {
    c(
      phi = th$phi,
      gamma = th$gamma,
      rho = th$rho,
      alpha = th$alpha,
      sigma2 = th$sigma2,
      eta = th$eta,
      zeta = th$zeta,
      chi = th$chi,
      omega = th$omega
    )
  })
)

post_df <- as.data.frame(post_df)

# 6) Posterior predictive simulation
X_pp_list <- vector("list", R_post)

say(sprintf(
  "Starting posterior predictive simulation: R_post = %d",
  R_post
))

for (r in seq_len(R_post)) {

  if (r %% 5 == 0 || r == 1 || r == R_post) {
    log_iter("posterior_predictive", r)

    cat(sprintf(
      "posterior_predictive %d / %d\n",
      r,
      R_post
    ))
  }

  X_pp_list[[r]] <- simulate_X_from_lambda(
    theta_pp[[r]],
    N_use,
    J_use,
    K_use
  )

  gc()
}

# 7) Save
saveRDS(
  list(
    prior_df = prior_df,
    post_df = post_df,
    X_obs = X_obs,
    X_pp_list = X_pp_list
  ),
  file = "results/ABC_ppc_cache.rds"
)

say("Saved results/ABC_ppc_cache.rds")
say("ABC_sim.R finished")
