suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(scales)
  library(io)
})

options(error = function() {
  traceback(2)
  quit(save = "no", status = 1)
})

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

progress_log <- "logs/ABC_plot_only.log"

if (file.exists(progress_log)) {
  invisible(file.remove(progress_log))
}

log_text <- function(txt) {
  cat(
    sprintf("%s\n", txt),
    file = progress_log,
    append = TRUE
  )
}

say <- function(txt) {
  cat(txt, "\n")
  log_text(txt)
}


say("ABC_plot_only.R started")


# ==================================================
# load plot-ready cache
# ==================================================

plot_cache_file <-
  "results/ABC_plot_data.rds"

source_cache_file <-
  "results/ABC_ppc_cache.rds"


if (!file.exists(plot_cache_file)) {
  
  stop(
    "Could not find ",
    plot_cache_file,
    ". Run ABC_prepare_plot_data.R first."
  )
}


say(sprintf(
  "Loading %s",
  plot_cache_file
))


cache <- readRDS(
  plot_cache_file
)


prior_df <- cache$prior_df
post_df <- cache$post_df
ppc <- cache$ppc


# Warn if ABC_ppc_cache.rds has been regenerated
# after the plot cache was prepared
if (
  file.exists(source_cache_file) &&
  !is.null(cache$source_cache_mtime)
) {
  
  current_source_mtime <-
    file.info(source_cache_file)$mtime
  
  if (
    !is.na(current_source_mtime) &&
    current_source_mtime >
    cache$source_cache_mtime
  ) {
    
    warning(
      paste0(
        source_cache_file,
        " is newer than ",
        plot_cache_file,
        ". If the posterior/PPC results changed, ",
        "rerun ABC_prepare_plot_data.R before plotting."
      )
    )
  }
}


# ==================================================
# prior/posterior plotting functions
# ==================================================

make_prior_posterior_hist <- function(
    x,
    breaks,
    xlim,
    title,
    xlab) {

  df <- data.frame(value = x)

  ggplot(df, aes(x = value)) +

    geom_histogram(
      aes(y = after_stat(density)),
      breaks = breaks,
      fill = "grey75",
      color = "black",
      linewidth = 0.45
    ) +

    scale_x_continuous(
      limits = xlim,
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +

    labs(
      title = title,
      x = xlab,
      y = "Density"
    ) +

    # Same plotting format as Benchmark 1
    theme_classic(base_size = 12) +

    theme(
      legend.position = "bottom",
      plot.title = element_text(
        size = 11,
        face = "bold"
      )
    )
}

make_prior_posterior_pair <- function(
    prior_x,
    post_x,
    prior_title,
    post_title,
    xlab,
    n_breaks = 61) {
  
  xlim <- range(
    c(
      prior_x,
      post_x
    ),
    finite = TRUE
  )
  
  
  if (!all(is.finite(xlim))) {
    
    stop(
      sprintf(
        "Non-finite plotting range for %s.",
        xlab
      )
    )
  }
  
  
  if (xlim[1] == xlim[2]) {
    
    xlim <- xlim +
      c(
        -0.5,
        0.5
      )
  }
  
  
  breaks <- seq(
    xlim[1],
    xlim[2],
    length.out = n_breaks
  )
  
  
  list(
    
    prior =
      make_prior_posterior_hist(
        x = prior_x,
        breaks = breaks,
        xlim = xlim,
        title = prior_title,
        xlab = xlab
      ),
    
    posterior =
      make_prior_posterior_hist(
        x = post_x,
        breaks = breaks,
        xlim = xlim,
        title = post_title,
        xlab = xlab
      )
  )
}


save_prior_posterior_page <- function(
    plot_list,
    file_out,
    width = 9,
    height = 10) {
  
  io::qdraw(
    
    gridExtra::grid.arrange(
      grobs = plot_list,
      ncol = 2
    ),
    
    file = file_out,
    device = NULL,
    width = width,
    height = height
  )
}


# ==================================================
# PPC plotting function
# ==================================================

make_two_panel_from_dfs <- function(
    df_raw,
    df_log,
    x_min_raw,
    x_max_raw,
    x_min_log,
    x_max_log,
    xlab_raw,
    xlab_log,
    title_base,
    file_out) {

  # --------------------------------------------------
  # log10(x + 1) transformation
  # --------------------------------------------------

  log10p1_trans <- scales::trans_new(
    name = "log10p1",
    transform = function(x) log10(x + 1),
    inverse = function(x) 10^x - 1,
    domain = c(-1 + 1e-12, Inf)
  )

  # --------------------------------------------------
  # RAW x-axis breaks
  # Keep complete numbers, but use fewer ticks so the
  # labels do not overlap.
  # --------------------------------------------------

  breaks_raw <- seq(
    from = x_min_raw,
    to = x_max_raw,
    length.out = 4
  )

  # --------------------------------------------------
  # LOG x-axis breaks
  # --------------------------------------------------

  log_min <- log10(1 + x_min_log)
  log_max <- log10(1 + x_max_log)

  breaks_log_space <- seq(
    from = log_min,
    to = log_max,
    length.out = 4
  )

  breaks_log_raw <- (10^breaks_log_space) - 1

  # Complete number labels: 5,000,000 rather than 5M
  full_number_labels <- scales::label_number(
    big.mark = ",",
    trim = TRUE
  )

  ymax_raw <- max(df_raw$prob, na.rm = TRUE)
  ymax_log <- max(df_log$prob, na.rm = TRUE)

  # ==================================================
  # RAW PANEL
  # ==================================================

  p_raw <- ggplot(df_raw) +

    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = 0,
        ymax = prob
      ),
      fill = "grey75",
      color = "black",
      linewidth = 0.55
    ) +

    # Same facet style as Benchmark 1.
    # Using facet_wrap also removes the old vertical
    # facet-strip line on the right side of the plots.
    facet_grid(
      rows = vars(source),
      axes = "all_x",
      axis.labels = "all_x"
    ) +

    scale_x_continuous(
      limits = c(x_min_raw, x_max_raw),
      breaks = breaks_raw,
      labels = full_number_labels,
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      breaks = function(lim) pretty(lim, n = 5),
      labels = label_number(accuracy = 0.01),
      expand = c(0, 0)
    ) +

    coord_cartesian(
      ylim = c(0, ymax_raw)
    ) +

    labs(
      title = paste0("RAW ", title_base),
      x = xlab_raw,
      y = "Probability"
    ) +

    # Same plotting format as Benchmark 1
    theme_classic(base_size = 12) +

    theme(
      legend.position = "bottom",
      plot.title = element_text(
        size = 12,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      strip.background = element_blank(),
      strip.text.y.right = element_text(
        angle = 270
      )
    )

  # ==================================================
  # LOG PANEL
  # ==================================================

  p_log <- ggplot(df_log) +

    geom_rect(
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = 0,
        ymax = prob
      ),
      fill = "grey75",
      color = "black",
      linewidth = 0.55
    ) +

    # Same facet style as Benchmark 1.
    facet_grid(
      rows = vars(source),
      axes = "all_x",
      axis.labels = "all_x"
    ) +

    scale_x_continuous(
      limits = c(x_min_log, x_max_log),
      trans = log10p1_trans,
      breaks = breaks_log_raw,
      labels = full_number_labels,
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      breaks = function(lim) pretty(lim, n = 5),
      labels = label_number(accuracy = 0.01),
      expand = c(0, 0)
    ) +

    coord_cartesian(
      ylim = c(0, ymax_log)
    ) +

    labs(
      title = paste0("LOG ", title_base),
      x = xlab_log,
      y = "Probability"
    ) +

    # Same plotting format as Benchmark 1
    theme_classic(base_size = 12) +

    theme(
      legend.position = "bottom",
      plot.title = element_text(
        size = 12,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      strip.background = element_blank(),
      strip.text.y.right = element_text(
        angle = 270
      )
    )

  # Two side-by-side panels. Each half is approximately
  # the same visual scale as the 7 x 5 Benchmark 1 plot,
  # while leaving enough room for full x-axis numbers.
  io::qdraw(
    gridExtra::grid.arrange(
      p_raw,
      p_log,
      ncol = 2
    ),
    file = file_out,
    device = NULL,
    width = 10,
    height = 5
  )
}

# ==================================================
# wrapper for stored PPC object
# ==================================================

draw_ppc <- function(
    obj,
    file_out) {
  
  make_two_panel_from_dfs(
    
    df_raw =
      obj$df_raw,
    
    df_log =
      obj$df_log,
    
    x_min_raw =
      obj$x_min_raw,
    
    x_max_raw =
      obj$x_max_raw,
    
    x_min_log =
      obj$x_min_log,
    
    x_max_log =
      obj$x_max_log,
    
    xlab_raw =
      obj$xlab_raw,
    
    xlab_log =
      obj$xlab_log,
    
    title_base =
      obj$title_base,
    
    file_out =
      file_out
  )
}


# ==================================================
# prior / posterior plots
# ==================================================

say("Making prior/posterior plots")


pair_phi <-
  make_prior_posterior_pair(
    
    prior_x =
      prior_df$phi,
    
    post_x =
      post_df$phi,
    
    prior_title =
      "Prior: gene-on probability (theta)",
    
    post_title =
      "Posterior: theta",
    
    xlab =
      "theta"
  )


pair_gamma <-
  make_prior_posterior_pair(
    
    prior_x =
      log(prior_df$gamma),
    
    post_x =
      log(post_df$gamma),
    
    prior_title =
      "Prior: biological variability log(gamma)",
    
    post_title =
      "Posterior: log(gamma)",
    
    xlab =
      "log(gamma)"
  )


pair_rho <-
  make_prior_posterior_pair(
    
    prior_x =
      log(prior_df$rho),
    
    post_x =
      log(post_df$rho),
    
    prior_title =
      "Prior: technical noise log(rho)",
    
    post_title =
      "Posterior: log(rho)",
    
    xlab =
      "log(rho)"
  )


save_prior_posterior_page(
  
  plot_list = list(
    pair_phi$prior,
    pair_phi$posterior,
    pair_gamma$prior,
    pair_gamma$posterior,
    pair_rho$prior,
    pair_rho$posterior
  ),
  
  file_out =
    "plots/ABC_prior_posterior_page1.png",
  
  width = 9,
  
  height = 10
)


pair_alpha <-
  make_prior_posterior_pair(
    
    prior_x =
      log(prior_df$alpha),
    
    post_x =
      log(post_df$alpha),
    
    prior_title =
      "Prior: TF activity variability log(alpha)",
    
    post_title =
      "Posterior: log(alpha)",
    
    xlab =
      "log(alpha)"
  )


pair_sigma2 <-
  make_prior_posterior_pair(
    
    prior_x =
      log(prior_df$sigma2),
    
    post_x =
      log(post_df$sigma2),
    
    prior_title =
      "Prior: regulatory effect variance log(sigma2)",
    
    post_title =
      "Posterior: log(sigma2)",
    
    xlab =
      "log(sigma2)"
  )


pair_eta <-
  make_prior_posterior_pair(
    
    prior_x =
      log(prior_df$eta),
    
    post_x =
      log(post_df$eta),
    
    prior_title =
      "Prior: target-specific edge propensity log(alpha_phi)",
    
    post_title =
      "Posterior: log(alpha_phi)",
    
    xlab =
      "log(alpha_phi)"
  )


save_prior_posterior_page(
  
  plot_list = list(
    pair_alpha$prior,
    pair_alpha$posterior,
    pair_sigma2$prior,
    pair_sigma2$posterior,
    pair_eta$prior,
    pair_eta$posterior
  ),
  
  file_out =
    "plots/ABC_prior_posterior_page2.png",
  
  width = 9,
  
  height = 10
)


pair_zeta <-
  make_prior_posterior_pair(
    
    prior_x =
      log(prior_df$zeta),
    
    post_x =
      log(post_df$zeta),
    
    prior_title =
      "Prior: target-specific edge propensity log(beta_phi)",
    
    post_title =
      "Posterior: log(beta_phi)",
    
    xlab =
      "log(beta_phi)"
  )


pair_chi <-
  make_prior_posterior_pair(
    
    prior_x =
      log(prior_df$chi),
    
    post_x =
      log(post_df$chi),
    
    prior_title =
      "Prior: regulator-specific edge propensity log(alpha_psi)",
    
    post_title =
      "Posterior: log(alpha_psi)",
    
    xlab =
      "log(alpha_psi)"
  )


pair_omega <-
  make_prior_posterior_pair(
    
    prior_x =
      log(prior_df$omega),
    
    post_x =
      log(post_df$omega),
    
    prior_title =
      "Prior: regulator-specific edge propensity log(beta_psi)",
    
    post_title =
      "Posterior: log(beta_psi)",
    
    xlab =
      "log(beta_psi)"
  )


save_prior_posterior_page(
  
  plot_list = list(
    pair_zeta$prior,
    pair_zeta$posterior,
    pair_chi$prior,
    pair_chi$posterior,
    pair_omega$prior,
    pair_omega$posterior
  ),
  
  file_out =
    "plots/ABC_prior_posterior_page3.png",
  
  width = 9,
  
  height = 10
)


# ==================================================
# PPC plots
# ==================================================

say("Making PPC plots")


say("Drawing entrywise all-count plot")

draw_ppc(
  ppc$entrywise_all,
  "plots/ABC_ppc_entrywise_all.png"
)


say("Drawing entrywise positive-count plot")

draw_ppc(
  ppc$entrywise_positive,
  "plots/ABC_ppc_entrywise_positive.png"
)


say("Drawing cell-level mean plot")

draw_ppc(
  ppc$cell_mean,
  "plots/ABC_ppc_cell_mean.png"
)


say("Drawing cell-level detection plot")

draw_ppc(
  ppc$cell_detect,
  "plots/ABC_ppc_cell_detect.png"
)


say("Drawing gene-level mean plot")

draw_ppc(
  ppc$gene_mean,
  "plots/ABC_ppc_gene_mean.png"
)


say("Drawing gene-level detection plot")

draw_ppc(
  ppc$gene_detect,
  "plots/ABC_ppc_gene_detect.png"
)


say("ABC_plot_only.R finished")