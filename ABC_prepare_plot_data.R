suppressPackageStartupMessages({
  library(Matrix)
})

options(error = function() {
  traceback(2)
  quit(save = "no", status = 1)
})

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

progress_log <- "logs/ABC_prepare_plot_data.log"

if (file.exists(progress_log)) {
  invisible(file.remove(progress_log))
}

log_text <- function(txt) {
  cat(sprintf("%s\n", txt), file = progress_log, append = TRUE)
}

say <- function(txt) {
  cat(txt, "\n")
  log_text(txt)
}

say("ABC_prepare_plot_data.R started")


# ==================================================
# files
# ==================================================

source_cache_file <- "results/ABC_ppc_cache.rds"
plot_cache_file <- "results/ABC_plot_data.rds"

if (!file.exists(source_cache_file)) {
  stop(
    "Could not find ",
    source_cache_file,
    ". Run the ABC posterior/posterior-predictive script first."
  )
}


# ==================================================
# load original ABC cache
# ==================================================

say(sprintf("Loading %s", source_cache_file))

cache <- readRDS(source_cache_file)

prior_df <- cache$prior_df
post_df <- cache$post_df
X_obs <- cache$X_obs
X_pp_list <- cache$X_pp_list


# ==================================================
# extract posterior-predictive matrices
# ==================================================

extract_pp_matrix <- function(obj) {
  
  # Old cache format: object itself is already a matrix
  if (is.matrix(obj) || inherits(obj, "Matrix")) {
    return(obj)
  }
  
  # New cache format: posterior-predictive object is a list
  if (is.list(obj)) {
    
    candidate_names <- c(
      "X",
      "X_sim",
      "counts",
      "X_counts"
    )
    
    for (nm in candidate_names) {
      
      if (!is.null(obj[[nm]]) &&
          (is.matrix(obj[[nm]]) ||
           inherits(obj[[nm]], "Matrix"))) {
        
        return(obj[[nm]])
      }
    }
  }
  
  stop(
    "Could not find a matrix in posterior-predictive object. Class = ",
    paste(class(obj), collapse = ", "),
    "; names = ",
    paste(names(obj), collapse = ", ")
  )
}


say("Extracting posterior-predictive count matrices")

X_pp_list <- lapply(
  X_pp_list,
  extract_pp_matrix
)

say(sprintf(
  "Prior samples = %d",
  nrow(prior_df)
))

say(sprintf(
  "Posterior samples = %d",
  nrow(post_df)
))

say(sprintf(
  "Posterior predictive datasets = %d",
  length(X_pp_list)
))


# ==================================================
# validation
# ==================================================

required_params <- c(
  "phi",
  "gamma",
  "rho",
  "alpha",
  "sigma2",
  "eta",
  "zeta",
  "chi",
  "omega"
)

if (!all(required_params %in% names(prior_df))) {
  stop(
    "prior_df is missing one or more required ABC hyperparameters."
  )
}

if (!all(required_params %in% names(post_df))) {
  stop(
    "post_df is missing one or more required ABC hyperparameters."
  )
}

if (length(X_pp_list) == 0L) {
  stop("X_pp_list is empty.")
}


# ==================================================
# helper functions
# ==================================================

hist_count_only <- function(x, breaks) {
  
  h <- hist(
    x,
    breaks = breaks,
    plot = FALSE,
    right = FALSE
  )
  
  as.numeric(h$counts)
}


# Histogram matrix entries without converting
# sparse matrices to dense matrices
matrix_hist_count <- function(X, breaks) {
  
  if (inherits(X, "sparseMatrix")) {
    
    vals <- X@x
    
    counts <- hist_count_only(
      vals,
      breaks
    )
    
    n_total <- nrow(X) * ncol(X)
    n_zero <- n_total - length(vals)
    
    if (n_zero > 0) {
      
      counts <- counts +
        n_zero * hist_count_only(
          0,
          breaks
        )
    }
    
    return(counts)
  }
  
  hist_count_only(
    as.numeric(X),
    breaks
  )
}


# Extract positive matrix entries without
# converting sparse matrices to dense matrices
positive_matrix_values <- function(X) {
  
  if (inherits(X, "sparseMatrix")) {
    vals <- X@x
  } else {
    vals <- as.numeric(X)
  }
  
  vals <- vals[
    is.finite(vals) &
      vals > 0
  ]
  
  vals
}


# ==================================================
# build plot-ready PPC object from vectors
# ==================================================

build_panel_from_obs_and_simlist <- function(
    obs,
    sim_list,
    title_base,
    xlab_raw,
    xlab_log,
    bins = 100) {
  
  obs <- as.numeric(obs)
  
  obs <- obs[
    is.finite(obs)
  ]
  
  sim_min <- Inf
  sim_max <- -Inf
  
  
  # ------------------------------------------------
  # determine common range
  # ------------------------------------------------
  
  for (i in seq_along(sim_list)) {
    
    if (
      i %% 5 == 0 ||
      i == 1 ||
      i == length(sim_list)
    ) {
      
      say(sprintf(
        "Range pass %d / %d for %s",
        i,
        length(sim_list),
        title_base
      ))
    }
    
    x <- as.numeric(
      sim_list[[i]]
    )
    
    x <- x[
      is.finite(x)
    ]
    
    if (length(x) > 0) {
      
      sim_min <- min(
        sim_min,
        min(x)
      )
      
      sim_max <- max(
        sim_max,
        max(x)
      )
    }
  }
  
  
  x_min <- min(
    c(obs, sim_min),
    na.rm = TRUE
  )
  
  x_max <- max(
    c(obs, sim_max),
    na.rm = TRUE
  )
  
  
  if (!is.finite(x_min) ||
      !is.finite(x_max)) {
    
    stop(
      "Could not determine plotting range for ",
      title_base
    )
  }
  
  
  if (x_min == x_max) {
    x_min <- x_min - 0.5
    x_max <- x_max + 0.5
  }
  
  
  # ------------------------------------------------
  # common histogram bins
  # ------------------------------------------------
  
  br_raw <- seq(
    x_min,
    x_max,
    length.out = bins + 1
  )
  
  br_log <- (
    10^seq(
      log10(1 + x_min),
      log10(1 + x_max),
      length.out = bins + 1
    )
  ) - 1
  
  
  # ------------------------------------------------
  # observed histogram
  # ------------------------------------------------
  
  obs_counts_raw <- hist_count_only(
    obs,
    br_raw
  )
  
  obs_counts_log <- hist_count_only(
    obs,
    br_log
  )
  
  
  # ------------------------------------------------
  # simulated histogram
  # ------------------------------------------------
  
  sim_counts_raw <- rep(
    0,
    length(obs_counts_raw)
  )
  
  sim_counts_log <- rep(
    0,
    length(obs_counts_log)
  )
  
  sim_total_raw <- 0
  sim_total_log <- 0
  
  
  for (i in seq_along(sim_list)) {
    
    if (
      i %% 5 == 0 ||
      i == 1 ||
      i == length(sim_list)
    ) {
      
      say(sprintf(
        "Histogram pass %d / %d for %s",
        i,
        length(sim_list),
        title_base
      ))
    }
    
    
    x <- as.numeric(
      sim_list[[i]]
    )
    
    x <- x[
      is.finite(x)
    ]
    
    
    c_raw <- hist_count_only(
      x,
      br_raw
    )
    
    c_log <- hist_count_only(
      x,
      br_log
    )
    
    
    sim_counts_raw <-
      sim_counts_raw + c_raw
    
    sim_counts_log <-
      sim_counts_log + c_log
    
    sim_total_raw <-
      sim_total_raw + sum(c_raw)
    
    sim_total_log <-
      sim_total_log + sum(c_log)
  }
  
  
  # ------------------------------------------------
  # plot-ready data frames
  # ------------------------------------------------
  
  df_obs_raw <- data.frame(
    xmin = br_raw[
      -length(br_raw)
    ],
    xmax = br_raw[-1],
    prob = obs_counts_raw /
      sum(obs_counts_raw),
    source = factor(
      "Observed",
      levels = c(
        "Observed",
        "Simulated"
      )
    )
  )
  
  
  df_sim_raw <- data.frame(
    xmin = br_raw[
      -length(br_raw)
    ],
    xmax = br_raw[-1],
    prob = sim_counts_raw /
      sim_total_raw,
    source = factor(
      "Simulated",
      levels = c(
        "Observed",
        "Simulated"
      )
    )
  )
  
  
  df_obs_log <- data.frame(
    xmin = br_log[
      -length(br_log)
    ],
    xmax = br_log[-1],
    prob = obs_counts_log /
      sum(obs_counts_log),
    source = factor(
      "Observed",
      levels = c(
        "Observed",
        "Simulated"
      )
    )
  )
  
  
  df_sim_log <- data.frame(
    xmin = br_log[
      -length(br_log)
    ],
    xmax = br_log[-1],
    prob = sim_counts_log /
      sim_total_log,
    source = factor(
      "Simulated",
      levels = c(
        "Observed",
        "Simulated"
      )
    )
  )
  
  
  list(
    df_raw = rbind(
      df_obs_raw,
      df_sim_raw
    ),
    
    df_log = rbind(
      df_obs_log,
      df_sim_log
    ),
    
    x_min_raw = x_min,
    x_max_raw = x_max,
    
    x_min_log = x_min,
    x_max_log = x_max,
    
    xlab_raw = xlab_raw,
    xlab_log = xlab_log,
    
    title_base = title_base
  )
}


# ==================================================
# positive-entry panel builder
# ==================================================

build_positive_panel_from_matrices <- function(
    obs,
    X_list,
    title_base,
    xlab_raw,
    xlab_log,
    bins = 100) {
  
  obs <- as.numeric(obs)
  
  obs <- obs[
    is.finite(obs)
  ]
  
  
  sim_min <- Inf
  sim_max <- -Inf
  
  
  # ------------------------------------------------
  # first pass: determine plotting range
  # ------------------------------------------------
  
  for (i in seq_along(X_list)) {
    
    if (
      i %% 5 == 0 ||
      i == 1 ||
      i == length(X_list)
    ) {
      
      say(sprintf(
        "Positive-range pass %d / %d for %s",
        i,
        length(X_list),
        title_base
      ))
    }
    
    
    x <- positive_matrix_values(
      X_list[[i]]
    )
    
    
    if (length(x) > 0) {
      
      sim_min <- min(
        sim_min,
        min(x)
      )
      
      sim_max <- max(
        sim_max,
        max(x)
      )
    }
  }
  
  
  x_min <- min(
    c(obs, sim_min),
    na.rm = TRUE
  )
  
  x_max <- max(
    c(obs, sim_max),
    na.rm = TRUE
  )
  
  
  if (!is.finite(x_min) ||
      !is.finite(x_max)) {
    
    stop(
      "Could not determine positive-entry range for ",
      title_base
    )
  }
  
  
  if (x_min == x_max) {
    x_min <- x_min - 0.5
    x_max <- x_max + 0.5
  }
  
  
  br_raw <- seq(
    x_min,
    x_max,
    length.out = bins + 1
  )
  
  
  br_log <- (
    10^seq(
      log10(1 + x_min),
      log10(1 + x_max),
      length.out = bins + 1
    )
  ) - 1
  
  
  obs_counts_raw <- hist_count_only(
    obs,
    br_raw
  )
  
  obs_counts_log <- hist_count_only(
    obs,
    br_log
  )
  
  
  sim_counts_raw <- rep(
    0,
    length(obs_counts_raw)
  )
  
  sim_counts_log <- rep(
    0,
    length(obs_counts_log)
  )
  
  sim_total_raw <- 0
  sim_total_log <- 0
  
  
  # ------------------------------------------------
  # second pass: histogram accumulation
  # ------------------------------------------------
  
  for (i in seq_along(X_list)) {
    
    if (
      i %% 5 == 0 ||
      i == 1 ||
      i == length(X_list)
    ) {
      
      say(sprintf(
        "Positive-histogram pass %d / %d for %s",
        i,
        length(X_list),
        title_base
      ))
    }
    
    
    x <- positive_matrix_values(
      X_list[[i]]
    )
    
    
    c_raw <- hist_count_only(
      x,
      br_raw
    )
    
    c_log <- hist_count_only(
      x,
      br_log
    )
    
    
    sim_counts_raw <-
      sim_counts_raw + c_raw
    
    sim_counts_log <-
      sim_counts_log + c_log
    
    sim_total_raw <-
      sim_total_raw + sum(c_raw)
    
    sim_total_log <-
      sim_total_log + sum(c_log)
  }
  
  
  df_obs_raw <- data.frame(
    xmin = br_raw[
      -length(br_raw)
    ],
    xmax = br_raw[-1],
    prob = obs_counts_raw /
      sum(obs_counts_raw),
    source = factor(
      "Observed",
      levels = c(
        "Observed",
        "Simulated"
      )
    )
  )
  
  
  df_sim_raw <- data.frame(
    xmin = br_raw[
      -length(br_raw)
    ],
    xmax = br_raw[-1],
    prob = sim_counts_raw /
      sim_total_raw,
    source = factor(
      "Simulated",
      levels = c(
        "Observed",
        "Simulated"
      )
    )
  )
  
  
  df_obs_log <- data.frame(
    xmin = br_log[
      -length(br_log)
    ],
    xmax = br_log[-1],
    prob = obs_counts_log /
      sum(obs_counts_log),
    source = factor(
      "Observed",
      levels = c(
        "Observed",
        "Simulated"
      )
    )
  )
  
  
  df_sim_log <- data.frame(
    xmin = br_log[
      -length(br_log)
    ],
    xmax = br_log[-1],
    prob = sim_counts_log /
      sim_total_log,
    source = factor(
      "Simulated",
      levels = c(
        "Observed",
        "Simulated"
      )
    )
  )
  
  
  list(
    df_raw = rbind(
      df_obs_raw,
      df_sim_raw
    ),
    
    df_log = rbind(
      df_obs_log,
      df_sim_log
    ),
    
    x_min_raw = x_min,
    x_max_raw = x_max,
    
    x_min_log = x_min,
    x_max_log = x_max,
    
    xlab_raw = xlab_raw,
    xlab_log = xlab_log,
    
    title_base = title_base
  )
}


# ==================================================
# observed summaries
# ==================================================

say("Preparing observed summaries")


# Positive observed counts
x_obs_nonzero <- positive_matrix_values(
  X_obs
)


# Cell-level mean UMI per gene
C_obs <- Matrix::rowSums(
  X_obs
) / ncol(X_obs)


# Cell-level fraction detected genes
D_obs <- Matrix::rowSums(
  sign(X_obs)
) / ncol(X_obs)


# Gene-level mean UMI per cell
G_obs <- Matrix::colSums(
  X_obs
) / nrow(X_obs)


# Gene-level fraction expressing cells
M_obs <- Matrix::colSums(
  sign(X_obs)
) / nrow(X_obs)


# ==================================================
# simulated summaries
# ==================================================

say("Preparing simulated summary lists")


C_sim_list <- vector(
  "list",
  length(X_pp_list)
)

D_sim_list <- vector(
  "list",
  length(X_pp_list)
)

G_sim_list <- vector(
  "list",
  length(X_pp_list)
)

M_sim_list <- vector(
  "list",
  length(X_pp_list)
)


for (i in seq_along(X_pp_list)) {
  
  if (
    i %% 5 == 0 ||
    i == 1 ||
    i == length(X_pp_list)
  ) {
    
    say(sprintf(
      "Summary pass %d / %d",
      i,
      length(X_pp_list)
    ))
  }
  
  
  X <- X_pp_list[[i]]
  
  
  C_sim_list[[i]] <-
    Matrix::rowSums(X) /
    ncol(X)
  
  
  D_sim_list[[i]] <-
    Matrix::rowSums(
      sign(X)
    ) /
    ncol(X)
  
  
  G_sim_list[[i]] <-
    Matrix::colSums(X) /
    nrow(X)
  
  
  M_sim_list[[i]] <-
    Matrix::colSums(
      sign(X)
    ) /
    nrow(X)
}


# ==================================================
# ENTRYWISE ALL
# ==================================================

say("Preparing entrywise all-count PPC")


obs_max <- max(X_obs)

sim_max_all <- 0


for (i in seq_along(X_pp_list)) {
  
  X <- X_pp_list[[i]]
  
  sim_max_all <- max(
    sim_max_all,
    max(X)
  )
}


bins <- 100

x_min_all <- 0

x_max_all <- max(
  obs_max,
  sim_max_all
)


if (x_max_all == x_min_all) {
  x_max_all <- x_min_all + 0.5
}


br_raw_all <- seq(
  x_min_all,
  x_max_all,
  length.out = bins + 1
)


br_log_all <- (
  10^seq(
    log10(1 + x_min_all),
    log10(1 + x_max_all),
    length.out = bins + 1
  )
) - 1


obs_counts_raw_all <-
  matrix_hist_count(
    X_obs,
    br_raw_all
  )


obs_counts_log_all <-
  matrix_hist_count(
    X_obs,
    br_log_all
  )


sim_counts_raw_all <- rep(
  0,
  length(obs_counts_raw_all)
)

sim_counts_log_all <- rep(
  0,
  length(obs_counts_log_all)
)


for (i in seq_along(X_pp_list)) {
  
  if (
    i %% 5 == 0 ||
    i == 1 ||
    i == length(X_pp_list)
  ) {
    
    say(sprintf(
      "Entrywise all histogram pass %d / %d",
      i,
      length(X_pp_list)
    ))
  }
  
  
  X <- X_pp_list[[i]]
  
  
  c_raw <- matrix_hist_count(
    X,
    br_raw_all
  )
  
  
  c_log <- matrix_hist_count(
    X,
    br_log_all
  )
  
  
  sim_counts_raw_all <-
    sim_counts_raw_all + c_raw
  
  sim_counts_log_all <-
    sim_counts_log_all + c_log
}


df_obs_raw_all <- data.frame(
  xmin = br_raw_all[
    -length(br_raw_all)
  ],
  xmax = br_raw_all[-1],
  prob = obs_counts_raw_all /
    sum(obs_counts_raw_all),
  source = factor(
    "Observed",
    levels = c(
      "Observed",
      "Simulated"
    )
  )
)


df_sim_raw_all <- data.frame(
  xmin = br_raw_all[
    -length(br_raw_all)
  ],
  xmax = br_raw_all[-1],
  prob = sim_counts_raw_all /
    sum(sim_counts_raw_all),
  source = factor(
    "Simulated",
    levels = c(
      "Observed",
      "Simulated"
    )
  )
)


df_obs_log_all <- data.frame(
  xmin = br_log_all[
    -length(br_log_all)
  ],
  xmax = br_log_all[-1],
  prob = obs_counts_log_all /
    sum(obs_counts_log_all),
  source = factor(
    "Observed",
    levels = c(
      "Observed",
      "Simulated"
    )
  )
)


df_sim_log_all <- data.frame(
  xmin = br_log_all[
    -length(br_log_all)
  ],
  xmax = br_log_all[-1],
  prob = sim_counts_log_all /
    sum(sim_counts_log_all),
  source = factor(
    "Simulated",
    levels = c(
      "Observed",
      "Simulated"
    )
  )
)


entrywise_all_plot_data <- list(
  
  df_raw = rbind(
    df_obs_raw_all,
    df_sim_raw_all
  ),
  
  df_log = rbind(
    df_obs_log_all,
    df_sim_log_all
  ),
  
  x_min_raw = x_min_all,
  x_max_raw = x_max_all,
  
  x_min_log = x_min_all,
  x_max_log = x_max_all,
  
  xlab_raw = "X_ij",
  xlab_log = "X_ij",
  
  title_base = "Entrywise: counts (all)"
)


# ==================================================
# ENTRYWISE POSITIVE
# ==================================================

say("Preparing entrywise positive-count PPC")


entrywise_positive_plot_data <-
  build_positive_panel_from_matrices(
    
    obs = x_obs_nonzero,
    
    X_list = X_pp_list,
    
    title_base =
      "Entrywise: counts (positive)",
    
    xlab_raw =
      "X_ij | X_ij > 0",
    
    xlab_log =
      "X_ij | X_ij > 0",
    
    bins = 100
  )


# ==================================================
# CELL-LEVEL MEAN
# ==================================================

say("Preparing cell-level mean PPC")


cell_mean_plot_data <-
  build_panel_from_obs_and_simlist(
    
    obs = C_obs,
    
    sim_list = C_sim_list,
    
    title_base =
      "Cell-level: mean UMI per gene",
    
    xlab_raw =
      "(1/J) sum_j X_ij",
    
    xlab_log =
      "(1/J) sum_j X_ij",
    
    bins = 100
  )


# ==================================================
# CELL-LEVEL DETECTION
# ==================================================

say("Preparing cell-level detection PPC")


cell_detect_plot_data <-
  build_panel_from_obs_and_simlist(
    
    obs = D_obs,
    
    sim_list = D_sim_list,
    
    title_base =
      "Cell-level: fraction detected genes",
    
    xlab_raw =
      "(1/J) sum_j 1{X_ij > 0}",
    
    xlab_log =
      "(1/J) sum_j 1{X_ij > 0}",
    
    bins = 100
  )


# ==================================================
# GENE-LEVEL MEAN
# ==================================================

say("Preparing gene-level mean PPC")


gene_mean_plot_data <-
  build_panel_from_obs_and_simlist(
    
    obs = G_obs,
    
    sim_list = G_sim_list,
    
    title_base =
      "Gene-level: mean UMI per cell",
    
    xlab_raw =
      "(1/N) sum_i X_ij",
    
    xlab_log =
      "(1/N) sum_i X_ij",
    
    bins = 100
  )


# ==================================================
# GENE-LEVEL DETECTION
# ==================================================

say("Preparing gene-level detection PPC")


gene_detect_plot_data <-
  build_panel_from_obs_and_simlist(
    
    obs = M_obs,
    
    sim_list = M_sim_list,
    
    title_base =
      "Gene-level: fraction expressing cells",
    
    xlab_raw =
      "(1/N) sum_i 1{X_ij > 0}",
    
    xlab_log =
      "(1/N) sum_i 1{X_ij > 0}",
    
    bins = 100
  )


# ==================================================
# save plot-ready cache
# ==================================================

say("Saving plot-ready cache")


plot_cache <- list(
  
  prior_df = prior_df,
  
  post_df = post_df,
  
  ppc = list(
    
    entrywise_all =
      entrywise_all_plot_data,
    
    entrywise_positive =
      entrywise_positive_plot_data,
    
    cell_mean =
      cell_mean_plot_data,
    
    cell_detect =
      cell_detect_plot_data,
    
    gene_mean =
      gene_mean_plot_data,
    
    gene_detect =
      gene_detect_plot_data
  ),
  
  source_cache_file =
    source_cache_file,
  
  source_cache_mtime =
    file.info(source_cache_file)$mtime
)


saveRDS(
  plot_cache,
  plot_cache_file
)


say(sprintf(
  "Saved plot-ready data to %s",
  plot_cache_file
))


say(
  "You can now redraw all figures using ABC_plot_only.R without extracting posterior-predictive matrices again."
)

say("ABC_prepare_plot_data.R finished")