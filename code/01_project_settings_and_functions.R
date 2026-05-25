###############################################################################
# Shared settings and functions for the MOEX liquidity-regime pipeline
#
# The numbered scripts are the public pipeline. This one shared file keeps
# settings, path handling, plotting style, and HMM functions consistent without
# requiring a separate settings folder.
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required.")
}

library(data.table)

split_env_vector <- function(value) {
  out <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
  out[nzchar(out)]
}

resolve_code_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_arg <- args_all[startsWith(args_all, file_arg)]

  if (length(script_arg) > 0L) {
    script_path <- sub(file_arg, "", script_arg[1L], fixed = TRUE)
    script_path <- gsub("~\\+~", " ", script_path)
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }

  current_dir <- normalizePath(getwd(), mustWork = TRUE)
  if (basename(current_dir) == "code") {
    return(current_dir)
  }
  code_candidate <- file.path(current_dir, "code")
  if (dir.exists(code_candidate)) {
    return(normalizePath(code_candidate, mustWork = TRUE))
  }
  current_dir
}

read_project_settings <- function(project_root) {
  # This list is the single source of truth for the empirical design. The later
  # scripts read these values rather than defining their own dates, tickers,
  # folders, response variables, or HMM settings.
  cfg <- list(
    # The baseline universe follows the MOEX Blue Chip sample used in the paper.
    # Tickers can be overridden by MOEX_TICKERS for small smoke tests.
    tickers = c("LKOH", "GAZP", "SBER", "TATN", "T", "NVTK", "GMKN", "ROSN", "YDEX", "PLZL", "SNGS", "CHMF", "NLMK", "MOEX", "HEAD"),
    sample_start_date = "2025-04-01",
    sample_end_date = "2025-06-30",
    period_label = "2025_april_june",
    # Raw order logs are first reconstructed on a 5-second grid. The statistical
    # model then uses a 30-second panel, with 5- and 60-second robustness checks.
    reconstruction_bucket_size_seconds = 5L,
    bucket_size_seconds = 30L,
    robustness_frequencies_seconds = c(5L, 60L),
    # Moscow Exchange main-session timing, in seconds after midnight. The open
    # and close cutoffs define the boundary dummies in the NHMM transition model.
    session_start_seconds = 36000L,
    session_end_seconds = 67200L,
    open_end_seconds = 37800L,
    close_start_seconds = 65400L,
    tod_bin_seconds = 300L,
    # Realized volatility is backward-looking: 12 buckets at 30 seconds is a
    # six-minute window, and the current return is deliberately excluded.
    rv_window_buckets = 12L,
    epsilon = 1e-12,
    # The primary emission is Gaussian. The Student-t setting is used only by
    # the robustness script, where it is treated as a sensitivity check.
    emission_family = "gaussian",
    response_variable = "r_t_winsorized",
    flow_variable = "ofi_depth_scaled",
    n_states = 3L,
    n_bootstrap_reps = 200L,
    bootstrap_maxit = 50L,
    em_tol = 1e-6,
    em_maxit = 200L,
    n_multistart = 4L,
    random_seed = 2025L,
    winsor_lower = 0.005,
    winsor_upper = 0.995,
    student_t_df = 5,
    drop_incomplete_rv = TRUE,
    process_all_selected_dates = TRUE,
    include_weekends = FALSE,
    reprocess_existing_days = FALSE,
    transformation_mode = "pooled_full_sample"
  )

  cfg$project_root <- normalizePath(project_root, mustWork = TRUE)
  cfg$code_dir <- file.path(cfg$project_root, "code")
  cfg$settings_source <- file.path(cfg$code_dir, "01_project_settings_and_functions.R")
  # By default, raw and processed data live beside the project. Environment
  # variables allow the same code to run from another machine or disk without
  # editing the scripts themselves.
  cfg$raw_data_root <- Sys.getenv("MOEX_RAW_DATA_ROOT", unset = file.path(cfg$project_root, "Данные"))
  cfg$results_root <- Sys.getenv("MOEX_RESULTS_ROOT", unset = file.path(cfg$project_root, "results"))
  # Processed data are split into daily ticker files and combined panel files.
  # The helper functions below enforce the file naming convention.
  cfg$daily_processed_root <- file.path(cfg$raw_data_root, "processed data", "daily")
  cfg$combined_data_root <- file.path(cfg$raw_data_root, "processed data", "combined")
  # Numbered result folders mirror the numbered research steps in the paper and
  # keep generated tables/figures out of the project root.
  cfg$results_descriptive_root <- file.path(cfg$results_root, "01_descriptive statistics")
  cfg$results_seasonality_root <- file.path(cfg$results_root, "02_seasonality")
  cfg$results_stationarity_root <- file.path(cfg$results_root, "03_stationarity")
  cfg$results_normality_root <- file.path(cfg$results_root, "04_normality")
  cfg$results_regressors_root <- file.path(cfg$results_root, "05_regressors")
  cfg$results_price_impact_root <- file.path(cfg$results_root, "06_price_impact")
  cfg$results_outliers_root <- file.path(cfg$results_root, "07_outliers")
  cfg$results_hmm_state_probabilities_root <- file.path(cfg$results_root, "08_hmm_state_probabilities")
  cfg$results_hmm_constant_transition_root <- file.path(cfg$results_root, "09_hmm_constant_transition")
  cfg$results_hmm_covariate_transition_root <- file.path(cfg$results_root, "10_hmm_covariate_transition")
  cfg$results_hmm_model_comparison_root <- file.path(cfg$results_root, "11_hmm_model_comparison")
  cfg$results_hmm_state_interpretation_root <- file.path(cfg$results_root, "12_hmm_state_interpretation")
  cfg$results_hmm_hypotheses_root <- file.path(cfg$results_root, "13_hmm_hypotheses")
  cfg$results_hmm_robustness_root <- file.path(cfg$results_root, "14_hmm_robustness")
  cfg$results_hmm_state_sequence_diagnostics_root <- file.path(cfg$results_root, "15_hmm_state_sequence_diagnostics")

  env_tickers <- split_env_vector(Sys.getenv("MOEX_TICKERS", unset = paste(cfg$tickers, collapse = ",")))
  if (length(env_tickers) > 0L) cfg$tickers <- env_tickers

  cfg$tickers <- as.character(cfg$tickers)
  cfg$sample_start_code <- gsub("-", "", cfg$sample_start_date)
  cfg$sample_end_code <- gsub("-", "", cfg$sample_end_date)
  cfg$reconstruction_bucket_size_seconds <- as.integer(cfg$reconstruction_bucket_size_seconds)
  cfg$bucket_size_seconds <- as.integer(cfg$bucket_size_seconds)
  cfg$robustness_frequencies_seconds <- as.integer(unlist(cfg$robustness_frequencies_seconds))
  cfg$session_start_seconds <- as.integer(cfg$session_start_seconds)
  cfg$session_end_seconds <- as.integer(cfg$session_end_seconds)
  cfg$open_end_seconds <- as.integer(cfg$open_end_seconds)
  cfg$close_start_seconds <- as.integer(cfg$close_start_seconds)
  cfg$tod_bin_seconds <- as.integer(cfg$tod_bin_seconds)
  cfg$rv_window_buckets <- as.integer(cfg$rv_window_buckets)
  cfg$n_states <- as.integer(cfg$n_states)
  cfg$n_bootstrap_reps <- as.integer(Sys.getenv("MOEX_NHMM_BOOTSTRAP_REPS", unset = as.character(cfg$n_bootstrap_reps)))
  cfg$bootstrap_maxit <- as.integer(Sys.getenv("MOEX_NHMM_BOOTSTRAP_MAXIT", unset = as.character(cfg$bootstrap_maxit)))
  cfg$em_maxit <- as.integer(Sys.getenv("MOEX_EM_MAXIT", unset = as.character(cfg$em_maxit)))
  cfg$n_multistart <- as.integer(Sys.getenv("MOEX_N_MULTISTART", unset = as.character(cfg$n_multistart)))
  cfg$random_seed <- as.integer(cfg$random_seed)
  cfg$student_t_df <- as.numeric(cfg$student_t_df)
  cfg$epsilon <- as.numeric(cfg$epsilon)
  cfg$winsor_lower <- as.numeric(cfg$winsor_lower)
  cfg$winsor_upper <- as.numeric(cfg$winsor_upper)
  cfg$em_tol <- as.numeric(cfg$em_tol)

  cfg
}

require_dir <- function(path, description = "folder") {
  if (!dir.exists(path)) stop("Required ", description, " not found: ", path)
  invisible(normalizePath(path, mustWork = TRUE))
}

write_table <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path)
  invisible(path)
}

append_table <- function(dt, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, path, append = file.exists(path), col.names = !file.exists(path))
  invisible(path)
}

processed_file_for_date <- function(cfg, ticker, date_string, bucket_size = cfg$reconstruction_bucket_size_seconds) {
  # Daily processed files are intentionally flat inside each ticker folder:
  # one CSV per stock-day, no extra date subfolders. This keeps the output easy
  # to inspect and matches the folder structure described in the paper.
  ticker_key <- tolower(ticker)
  d <- as.Date(date_string, format = "%Y%m%d")
  day_label <- format(d, "%d.%m.%Y")
  file.path(
    cfg$daily_processed_root,
    ticker_key,
    sprintf("processed_%ss_data_%s_%s.csv", bucket_size, ticker_key, day_label)
  )
}

combined_data_file <- function(cfg, frequency = cfg$bucket_size_seconds, trade_buckets_only = FALSE, suffix = NULL) {
  # Combined files encode the frequency and whether only trade buckets are kept.
  # The 5-second all-bucket and 5-second trade-bucket files are both useful:
  # the first is the raw reconstruction panel, while the second is the
  # robustness sample that removes no-trade intervals.
  prefix <- if (length(cfg$tickers) == 1L) tolower(cfg$tickers) else "moex_panel"
  extra <- if (!is.null(suffix)) paste0("_", suffix) else if (trade_buckets_only) "_trade_buckets" else ""
  file.path(
    cfg$combined_data_root,
    sprintf("%s_%s_processed_%ss%s_data.csv", prefix, cfg$period_label, as.integer(frequency), extra)
  )
}

paper_colours <- list(
  calm = "#0047FF",
  normal = "#8A8A8A",
  stressed = "#FF2A2A",
  purple = "#6C2EFF",
  orange = "#FF9A00",
  green = "#00A65A",
  cyan = "#00A6D6",
  light_blue = "#DCE6FF",
  light_red = "#FFE0E0",
  dark = "#111111"
)

paper_theme <- function(base_size = 12) {
  # All generated figures use the same serif theme and color system so that the
  # figures can be included directly in the LaTeX paper without manual restyling.
  ggplot2::theme_bw(base_size = base_size, base_family = "serif") +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(colour = paper_colours$dark),
      axis.text = ggplot2::element_text(colour = "#222222"),
      axis.ticks = ggplot2::element_line(colour = "#222222", linewidth = 0.35),
      panel.border = ggplot2::element_rect(colour = "#222222", fill = NA, linewidth = 0.55),
      panel.grid.major = ggplot2::element_line(linewidth = 0.22, colour = "#DDDDDD"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", colour = "#222222", linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "bold", colour = paper_colours$dark),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

save_paper_plot <- function(plot_object, path, width = 8.8, height = 5.0, dpi = 220) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot_object, width = width, height = height, dpi = dpi, bg = "white")
  invisible(path)
}

clean_numeric <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

safe_sd <- function(x, fallback = 1) {
  x <- clean_numeric(x)
  s <- stats::sd(x)
  if (is.na(s) || s <= 1e-10) fallback else s
}

kurtosis_manual <- function(x) {
  x <- clean_numeric(x)
  if (length(x) < 4L) {
    return(NA_real_)
  }
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  mean(((x - mean(x)) / s)^4)
}

winsorize_by_group <- function(dt, value_col, output_col, group_col, lower, upper) {
  # Winsorization is done within stock, not globally, so one very volatile stock
  # does not set the clipping thresholds for the whole panel.
  bounds <- dt[, .(
    lo = stats::quantile(get(value_col), lower, na.rm = TRUE, names = FALSE),
    hi = stats::quantile(get(value_col), upper, na.rm = TRUE, names = FALSE)
  ), by = group_col]
  dt[bounds, (output_col) := pmax(i.lo, pmin(i.hi, get(value_col))), on = group_col]
  invisible(dt)
}

standardize_by_group <- function(dt, value_col, output_col, by_cols) {
  # Standardization removes each stock's usual intraday level and scale. For
  # transition covariates, this makes a one-unit AME roughly a one-standard-
  # deviation shock within the stock/time-of-day cell.
  dt[, (output_col) := {
    x <- get(value_col)
    s <- stats::sd(x, na.rm = TRUE)
    if (is.na(s) || s <= 0) rep(0, .N) else (x - mean(x, na.rm = TRUE)) / s
  }, by = by_cols]
  invisible(dt)
}

center_by_group <- function(dt, value_col, output_col, by_cols) {
  # Spread and depth levels are centered by time of day to remove deterministic
  # intraday patterns while keeping economically meaningful variation.
  dt[, (output_col) := get(value_col) - mean(get(value_col), na.rm = TRUE), by = by_cols]
  invisible(dt)
}

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

softmax_row <- function(eta) {
  eta <- eta - max(eta)
  out <- exp(eta)
  out / sum(out)
}

stationary_distribution <- function(P) {
  k <- nrow(P)
  A <- t(P) - diag(k)
  A[k, ] <- 1
  b <- c(rep(0, k - 1L), 1)
  out <- tryCatch(as.numeric(solve(A, b)), error = function(e) rep(1 / k, k))
  out <- pmax(out, 1e-8)
  out / sum(out)
}

state_names <- function(k) {
  if (k == 2L) {
    return(c("low_cost", "high_cost"))
  }
  if (k == 3L) {
    return(c("calm", "normal", "stressed"))
  }
  paste0("state_", seq_len(k))
}

state_labels <- function(k) {
  if (k == 2L) {
    return(c("Low-cost", "High-cost / high-impact"))
  }
  if (k == 3L) {
    return(c("Calm", "Normal", "Stressed"))
  }
  paste("State", seq_len(k))
}

state_palette <- function(k) {
  if (k == 2L) {
    return(c("Low-cost" = paper_colours$calm, "High-cost / high-impact" = paper_colours$stressed))
  }
  if (k == 3L) {
    return(c("Calm" = paper_colours$calm, "Normal" = paper_colours$normal, "Stressed" = paper_colours$stressed))
  }
  stats::setNames(grDevices::rainbow(k), paste("State", seq_len(k)))
}

make_start_labels <- function(data, k, seed = 1L, jitter = FALSE) {
  # EM needs starting states. The score below is only an initialization device:
  # high spreads, large absolute returns, low depth, and higher realized
  # volatility are used to create rough starting groups before the likelihood
  # estimation takes over.
  if (k == 1L) {
    return(rep(1L, nrow(data)))
  }
  score <- as.numeric(scale(data$spread_tilde)) +
    0.35 * as.numeric(scale(abs(data$response))) -
    0.25 * as.numeric(scale(data$depth_tilde)) +
    0.20 * as.numeric(scale(data$rv_tilde))
  score[!is.finite(score)] <- 0
  if (jitter) {
    set.seed(seed)
    score <- score + stats::rnorm(length(score), sd = 0.05)
  }
  breaks <- unique(stats::quantile(score, probs = seq(0, 1, length.out = k + 1L), na.rm = TRUE))
  if (length(breaks) <= k) {
    return(sample(seq_len(k), nrow(data), replace = TRUE))
  }
  as.integer(cut(score, breaks = breaks, include.lowest = TRUE, labels = FALSE))
}

start_from_labels <- function(labels, data, k, covariance_mode = "full") {
  # Convert rough starting labels into HMM parameters: one weighted regression
  # for the return equation, one spread intercept, one covariance matrix, and
  # a transition matrix from adjacent labels within each stock-day sequence.
  states <- vector("list", k)
  x_all <- cbind(1, data$flow_scaled)
  for (state in seq_len(k)) {
    x <- data[labels == state]
    if (nrow(x) < max(50L, 4L * k)) x <- data
    X <- cbind(1, x$flow_scaled)
    fit <- tryCatch(stats::lm.fit(X, x$response), error = function(e) NULL)
    beta <- if (is.null(fit)) c(mean(x$response), 0) else as.numeric(fit$coefficients)
    beta[!is.finite(beta)] <- 0
    resid_s <- x$spread_tilde - mean(x$spread_tilde, na.rm = TRUE)
    resid_m <- x$response - as.numeric(X %*% beta)
    rho <- suppressWarnings(stats::cor(resid_s, resid_m, use = "complete.obs"))
    if (!is.finite(rho) || identical(covariance_mode, "diagonal")) rho <- 0
    states[[state]] <- list(
      alpha_s = mean(x$spread_tilde, na.rm = TRUE),
      alpha_m = beta[1L],
      lambda = beta[2L],
      sd_s = max(safe_sd(resid_s), 1e-4),
      sd_m = max(safe_sd(resid_m), 1e-4),
      rho = pmin(pmax(rho, -0.95), 0.95)
    )
  }

  P <- matrix(1, nrow = k, ncol = k)
  for (idx in split(seq_along(labels), data$sequence_id)) {
    y <- labels[idx]
    if (length(y) > 1L) {
      for (tt in seq_len(length(y) - 1L)) P[y[tt], y[tt + 1L]] <- P[y[tt], y[tt + 1L]] + 1
    }
  }
  P <- P / rowSums(P)
  list(states = states, P = P, pi = stationary_distribution(P))
}

emission_logdens <- function(states, data, emission_family = "gaussian", covariance_mode = "full", df = 5) {
  # Emission likelihood for the bivariate observation (spread, return). The
  # return mean includes the state-specific slope on scaled depth-adjusted OFI.
  k <- length(states)
  out <- matrix(NA_real_, nrow = nrow(data), ncol = k)
  for (state in seq_len(k)) {
    p <- states[[state]]
    e_s <- (data$spread_tilde - p$alpha_s) / p$sd_s
    e_m <- (data$response - p$alpha_m - p$lambda * data$flow_scaled) / p$sd_m
    rho <- if (identical(covariance_mode, "diagonal")) 0 else p$rho
    one_minus_rho2 <- max(1 - rho^2, 1e-8)
    quad <- (e_s^2 - 2 * rho * e_s * e_m + e_m^2) / one_minus_rho2
    log_det <- 2 * log(p$sd_s) + 2 * log(p$sd_m) + log(one_minus_rho2)
    if (identical(emission_family, "student_t")) {
      dimension <- 2
      out[, state] <- lgamma((df + dimension) / 2) - lgamma(df / 2) -
        (dimension / 2) * log(df * pi) - 0.5 * log_det -
        ((df + dimension) / 2) * log1p(quad / df)
    } else {
      out[, state] <- -log(2 * pi) - 0.5 * log_det - 0.5 * quad
    }
  }
  out
}

student_t_weights <- function(states, data, covariance_mode = "full", df = 5) {
  k <- length(states)
  out <- matrix(1, nrow = nrow(data), ncol = k)
  for (state in seq_len(k)) {
    p <- states[[state]]
    e_s <- (data$spread_tilde - p$alpha_s) / p$sd_s
    e_m <- (data$response - p$alpha_m - p$lambda * data$flow_scaled) / p$sd_m
    rho <- if (identical(covariance_mode, "diagonal")) 0 else p$rho
    one_minus_rho2 <- max(1 - rho^2, 1e-8)
    quad <- (e_s^2 - 2 * rho * e_s * e_m + e_m^2) / one_minus_rho2
    # In the Student-t robustness check, large Mahalanobis residuals receive
    # lower weights in the emission update rather than forcing a new state.
    out[, state] <- (df + 2) / (df + pmax(quad, 0))
  }
  out
}

forward_backward_constant <- function(states, P, data, indices_by_sequence, emission_family = "gaussian", covariance_mode = "full", df = 5) {
  # E-step for the constant-transition HMM. It returns smoothed state
  # probabilities (gamma) and expected transition counts (xi_sum).
  k <- length(states)
  log_em <- emission_logdens(states, data, emission_family, covariance_mode, df)
  log_P <- log(pmax(P, 1e-12))
  log_pi <- log(stationary_distribution(P))
  gamma_store <- matrix(NA_real_, nrow = nrow(data), ncol = k)
  alpha_store <- matrix(NA_real_, nrow = nrow(data), ncol = k)
  beta_store <- matrix(NA_real_, nrow = nrow(data), ncol = k)
  xi_sum <- matrix(0, nrow = k, ncol = k)
  total_loglik <- 0

  for (idx in indices_by_sequence) {
    n <- length(idx)
    alpha <- matrix(NA_real_, nrow = n, ncol = k)
    beta <- matrix(0, nrow = n, ncol = k)
    scales <- rep(NA_real_, n)

    # Alpha and beta are kept on the log scale and centered by per-time scales.
    # This avoids numerical underflow on long stock-day sequences.
    raw <- log_pi + log_em[idx[1L], ]
    scales[1L] <- logsumexp(raw)
    alpha[1L, ] <- raw - scales[1L]

    if (n > 1L) {
      for (tt in 2:n) {
        raw <- vapply(seq_len(k), function(j) {
          log_em[idx[tt], j] + logsumexp(alpha[tt - 1L, ] + log_P[, j])
        }, numeric(1L))
        scales[tt] <- logsumexp(raw)
        alpha[tt, ] <- raw - scales[tt]
      }
      for (tt in (n - 1L):1L) {
        beta[tt, ] <- vapply(seq_len(k), function(i) {
          logsumexp(log_P[i, ] + log_em[idx[tt + 1L], ] + beta[tt + 1L, ]) - scales[tt + 1L]
        }, numeric(1L))
      }
      for (tt in seq_len(n - 1L)) {
        log_xi <- matrix(NA_real_, nrow = k, ncol = k)
        for (i in seq_len(k)) {
          for (j in seq_len(k)) {
            log_xi[i, j] <- alpha[tt, i] + log_P[i, j] + log_em[idx[tt + 1L], j] +
              beta[tt + 1L, j] - scales[tt + 1L]
          }
        }
        xi_sum <- xi_sum + exp(log_xi)
      }
    }

    log_gamma <- alpha + beta
    norm <- apply(log_gamma, 1L, logsumexp)
    gamma <- exp(log_gamma - norm)
    gamma_store[idx, ] <- gamma
    alpha_store[idx, ] <- alpha
    beta_store[idx, ] <- beta
    total_loglik <- total_loglik + sum(scales)
  }

  list(logLik = total_loglik, gamma = gamma_store, alpha = alpha_store, beta = beta_store, xi_sum = xi_sum, log_emission = log_em)
}

update_emissions <- function(post, data, old_states = NULL, emission_family = "gaussian", covariance_mode = "full", df = 5) {
  # M-step for emission parameters. Posterior probabilities are used as weights,
  # so observations can partly belong to several states rather than being hard-
  # assigned during estimation.
  k <- ncol(post$gamma)
  states <- vector("list", k)
  t_weights <- if (identical(emission_family, "student_t") && !is.null(old_states)) {
    student_t_weights(old_states, data, covariance_mode, df)
  } else {
    matrix(1, nrow = nrow(data), ncol = k)
  }

  X <- cbind(1, data$flow_scaled)
  for (state in seq_len(k)) {
    tau <- post$gamma[, state]
    u <- t_weights[, state]
    w_mu <- tau * u
    denom_mu <- sum(w_mu)
    alpha_s <- sum(w_mu * data$spread_tilde) / denom_mu
    beta_hat <- as.numeric(solve(crossprod(X, X * w_mu) + diag(1e-8, ncol(X)), crossprod(X, data$response * w_mu)))
    resid_s <- data$spread_tilde - alpha_s
    resid_m <- data$response - as.numeric(X %*% beta_hat)

    # For Student-t ECME, the covariance denominator is sum(tau), not sum(tau*u).
    denom_sigma <- if (identical(emission_family, "student_t")) sum(tau) else denom_mu
    var_s <- sum(tau * u * resid_s^2) / denom_sigma
    var_m <- sum(tau * u * resid_m^2) / denom_sigma
    cov_sm <- sum(tau * u * resid_s * resid_m) / denom_sigma
    sd_s <- sqrt(max(var_s, 1e-6))
    sd_m <- sqrt(max(var_m, 1e-6))
    rho <- if (identical(covariance_mode, "diagonal")) 0 else cov_sm / (sd_s * sd_m)

    states[[state]] <- list(
      alpha_s = alpha_s,
      alpha_m = beta_hat[1L],
      lambda = beta_hat[2L],
      sd_s = sd_s,
      sd_m = sd_m,
      rho = pmin(pmax(rho, -0.95), 0.95)
    )
  }

  states
}

fit_constant_hmm <- function(data, k, maxit, tolerance, n_starts = 4L, seed = 1L, emission_family = "gaussian", covariance_mode = "full", df = 5) {
  # Run EM from several starting partitions and keep the highest likelihood
  # solution. This helps guard against local optima in a non-convex HMM problem.
  indices <- split(seq_len(nrow(data)), data$sequence_id)
  best_fit <- NULL
  start_rows <- vector("list", n_starts)
  for (s in seq_len(n_starts)) {
    labels <- make_start_labels(data, k, seed = seed + s, jitter = s > 1L)
    init <- start_from_labels(labels, data, k, covariance_mode)
    states <- init$states
    P <- init$P
    post <- forward_backward_constant(states, P, data, indices, emission_family, covariance_mode, df)
    old_loglik <- post$logLik
    best_iter_loglik <- old_loglik
    best_iter_states <- states
    best_iter_P <- P
    trace <- data.table(iteration = 0L, logLik = old_loglik, improvement = NA_real_, is_negative = FALSE)
    convergence <- 1L
    convergence_message <- "maximum iterations reached"

    for (iter in seq_len(maxit)) {
      states_new <- update_emissions(post, data, states, emission_family, covariance_mode, df)
      P_new <- post$xi_sum + 1e-8
      P_new <- P_new / rowSums(P_new)
      post_new <- forward_backward_constant(states_new, P_new, data, indices, emission_family, covariance_mode, df)
      improvement <- post_new$logLik - old_loglik
      trace <- rbind(trace, data.table(iteration = iter, logLik = post_new$logLik, improvement = improvement, is_negative = improvement < -1e-6))

      states <- states_new
      P <- P_new
      post <- post_new
      if (post_new$logLik > best_iter_loglik) {
        best_iter_loglik <- post_new$logLik
        best_iter_states <- states_new
        best_iter_P <- P_new
      }
      if (is.finite(improvement) && abs(improvement) <= tolerance * (1 + abs(old_loglik))) {
        convergence <- 0L
        convergence_message <- "relative log-likelihood tolerance reached"
        break
      }
      old_loglik <- post_new$logLik
    }

    post_best <- forward_backward_constant(best_iter_states, best_iter_P, data, indices, emission_family, covariance_mode, df)
    fit <- list(
      states = best_iter_states,
      P = best_iter_P,
      pi = stationary_distribution(best_iter_P),
      post = post_best,
      trace = trace,
      logLik = post_best$logLik,
      convergence = convergence,
      message = convergence_message,
      start = s
    )
    start_rows[[s]] <- data.table(start = s, logLik = fit$logLik, iterations = max(trace$iteration), convergence = convergence)
    if (is.null(best_fit) || fit$logLik > best_fit$logLik) best_fit <- fit
  }
  best_fit$start_convergence <- rbindlist(start_rows)
  relabel_fit(best_fit)
}

relabel_fit <- function(fit) {
  k <- length(fit$states)
  ord <- order(vapply(fit$states, function(s) s$alpha_s, numeric(1L)))
  fit$states <- fit$states[ord]
  if (!is.null(fit$P)) fit$P <- fit$P[ord, ord, drop = FALSE]
  if (!is.null(fit$pi)) fit$pi <- fit$pi[ord]
  if (!is.null(fit$post$gamma)) fit$post$gamma <- fit$post$gamma[, ord, drop = FALSE]
  if (!is.null(fit$post$alpha)) fit$post$alpha <- fit$post$alpha[, ord, drop = FALSE]
  if (!is.null(fit$post$beta)) fit$post$beta <- fit$post$beta[, ord, drop = FALSE]
  if (!is.null(fit$post$xi_sum)) fit$post$xi_sum <- fit$post$xi_sum[ord, ord, drop = FALSE]
  if (!is.null(fit$beta_list)) {
    old_beta <- fit$beta_list
    fit$beta_list <- lapply(seq_along(ord), function(new_origin) {
      # Re-express the old destination logits under the new state ordering.
      # Each origin stores K destination linear predictors with destination 1 as
      # the zero baseline. After relabeling, the new destination 1 must again be
      # the zero baseline, so all rows are shifted by that row.
      b <- old_beta[[ord[new_origin]]][ord, , drop = FALSE]
      b - matrix(b[1L, ], nrow = nrow(b), ncol = ncol(b), byrow = TRUE)
    })
  }
  fit$state_order <- ord
  fit
}

transition_probs_for_row <- function(x_row, beta_list) {
  k <- length(beta_list)
  P <- matrix(NA_real_, nrow = k, ncol = k)
  for (i in seq_len(k)) {
    eta <- as.numeric(beta_list[[i]] %*% x_row)
    P[i, ] <- softmax_row(eta)
  }
  P
}

forward_backward_nhmm <- function(states, beta_list, data, indices_by_sequence, design, emission_family = "gaussian", covariance_mode = "full", df = 5) {
  # E-step for the covariate-transition model. The emission parameters are the
  # same type as in the constant HMM, but each transition matrix is now computed
  # from X_t through multinomial logits.
  k <- length(states)
  log_em <- emission_logdens(states, data, emission_family, covariance_mode, df)
  first_indices <- unlist(lapply(indices_by_sequence, `[`, 1L), use.names = FALSE)
  avg_P <- Reduce("+", lapply(first_indices, function(idx) transition_probs_for_row(design[idx, ], beta_list))) / length(first_indices)
  log_pi <- log(stationary_distribution(avg_P))
  gamma_store <- matrix(NA_real_, nrow = nrow(data), ncol = k)
  alpha_store <- matrix(NA_real_, nrow = nrow(data), ncol = k)
  beta_store <- matrix(NA_real_, nrow = nrow(data), ncol = k)
  xi_array <- array(0, dim = c(nrow(data), k, k))
  total_loglik <- 0

  for (idx in indices_by_sequence) {
    n <- length(idx)
    alpha <- matrix(NA_real_, nrow = n, ncol = k)
    beta <- matrix(0, nrow = n, ncol = k)
    scales <- rep(NA_real_, n)
    raw <- log_pi + log_em[idx[1L], ]
    scales[1L] <- logsumexp(raw)
    alpha[1L, ] <- raw - scales[1L]

    if (n > 1L) {
      for (tt in 2:n) {
        log_P <- log(pmax(transition_probs_for_row(design[idx[tt - 1L], ], beta_list), 1e-12))
        raw <- vapply(seq_len(k), function(j) {
          log_em[idx[tt], j] + logsumexp(alpha[tt - 1L, ] + log_P[, j])
        }, numeric(1L))
        scales[tt] <- logsumexp(raw)
        alpha[tt, ] <- raw - scales[tt]
      }
      for (tt in (n - 1L):1L) {
        log_P <- log(pmax(transition_probs_for_row(design[idx[tt], ], beta_list), 1e-12))
        beta[tt, ] <- vapply(seq_len(k), function(i) {
          logsumexp(log_P[i, ] + log_em[idx[tt + 1L], ] + beta[tt + 1L, ]) - scales[tt + 1L]
        }, numeric(1L))
      }
      for (tt in seq_len(n - 1L)) {
        log_P <- log(pmax(transition_probs_for_row(design[idx[tt], ], beta_list), 1e-12))
        log_xi <- matrix(NA_real_, nrow = k, ncol = k)
        for (i in seq_len(k)) {
          for (j in seq_len(k)) {
            log_xi[i, j] <- alpha[tt, i] + log_P[i, j] + log_em[idx[tt + 1L], j] +
              beta[tt + 1L, j] - scales[tt + 1L]
          }
        }
        xi_array[idx[tt], , ] <- exp(log_xi)
      }
    }

    log_gamma <- alpha + beta
    norm <- apply(log_gamma, 1L, logsumexp)
    gamma <- exp(log_gamma - norm)
    gamma_store[idx, ] <- gamma
    alpha_store[idx, ] <- alpha
    beta_store[idx, ] <- beta
    total_loglik <- total_loglik + sum(scales)
  }

  list(logLik = total_loglik, gamma = gamma_store, alpha = alpha_store, beta = beta_store, xi_array = xi_array, log_emission = log_em)
}

initial_beta_from_P <- function(P, design_names) {
  # The NHMM starts from the constant-HMM transition matrix. Intercepts are set
  # so the initial covariate model reproduces those baseline probabilities when
  # covariates are zero.
  k <- nrow(P)
  p <- length(design_names)
  lapply(seq_len(k), function(i) {
    beta <- matrix(0, nrow = k, ncol = p, dimnames = list(paste0("to_", seq_len(k)), design_names))
    for (j in 2:k) beta[j, "(Intercept)"] <- log(pmax(P[i, j], 1e-8) / pmax(P[i, 1L], 1e-8))
    beta
  })
}

fit_multinom_origin <- function(design, xi_i, beta_previous, maxit = 80L) {
  # For one origin state, estimate destination probabilities as a multinomial
  # logit using posterior transition weights as fractional observations.
  k <- ncol(xi_i)
  valid <- rowSums(xi_i) > 1e-12 & complete.cases(design)
  if (sum(valid) < ncol(design) + k + 20L) {
    return(beta_previous)
  }

  y <- xi_i[valid, , drop = FALSE]
  x <- as.data.frame(design[valid, , drop = FALSE])
  names(x) <- make.names(colnames(design), unique = TRUE)
  formula_terms <- paste(names(x), collapse = " + ")
  form <- stats::as.formula(paste("y ~ 0 +", formula_terms))

  fit <- tryCatch(
    suppressWarnings(nnet::multinom(form, data = x, trace = FALSE, maxit = maxit, MaxNWts = 5000)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(beta_previous)
  }

  coefs <- stats::coef(fit)
  if (is.null(dim(coefs))) coefs <- matrix(coefs, nrow = 1L)
  beta <- matrix(0, nrow = k, ncol = ncol(design), dimnames = dimnames(beta_previous))
  row_count <- min(nrow(coefs), k - 1L)
  if (row_count > 0L) beta[seq_len(row_count) + 1L, ] <- coefs[seq_len(row_count), seq_len(ncol(design)), drop = FALSE]
  beta[!is.finite(beta)] <- 0
  beta
}

update_nhmm_transitions <- function(design, xi_array, beta_previous) {
  k <- dim(xi_array)[2L]
  lapply(seq_len(k), function(i) {
    # Each origin state has its own multinomial logit. The posterior transition
    # weights from the E-step act as fractional destination counts.
    fit_multinom_origin(design, xi_array[, i, , drop = TRUE], beta_previous[[i]])
  })
}

fit_nhmm <- function(data, k, baseline_fit, design, maxit, tolerance, seed = 1L, emission_family = "gaussian", covariance_mode = "full", df = 5) {
  # Full EM loop for the non-homogeneous HMM: update emissions, then update the
  # transition logits, then recompute posteriors under time-varying transitions.
  set.seed(seed)
  indices <- split(seq_len(nrow(data)), data$sequence_id)
  states <- baseline_fit$states
  beta_list <- initial_beta_from_P(baseline_fit$P, colnames(design))
  post <- forward_backward_nhmm(states, beta_list, data, indices, design, emission_family, covariance_mode, df)
  old_loglik <- post$logLik
  best_loglik <- old_loglik
  best_states <- states
  best_beta <- beta_list
  trace <- data.table(iteration = 0L, logLik = old_loglik, improvement = NA_real_, is_negative = FALSE)
  convergence <- 1L
  convergence_message <- "maximum iterations reached"

  for (iter in seq_len(maxit)) {
    states_new <- update_emissions(post, data, states, emission_family, covariance_mode, df)
    beta_new <- update_nhmm_transitions(design, post$xi_array, beta_list)
    post_new <- forward_backward_nhmm(states_new, beta_new, data, indices, design, emission_family, covariance_mode, df)
    improvement <- post_new$logLik - old_loglik
    trace <- rbind(trace, data.table(iteration = iter, logLik = post_new$logLik, improvement = improvement, is_negative = improvement < -1e-6))

    states <- states_new
    beta_list <- beta_new
    post <- post_new
    if (post_new$logLik > best_loglik) {
      best_loglik <- post_new$logLik
      best_states <- states_new
      best_beta <- beta_new
    }
    if (is.finite(improvement) && abs(improvement) <= tolerance * (1 + abs(old_loglik))) {
      convergence <- 0L
      convergence_message <- "relative log-likelihood tolerance reached"
      break
    }
    old_loglik <- post_new$logLik
  }

  best_post <- forward_backward_nhmm(best_states, best_beta, data, indices, design, emission_family, covariance_mode, df)
  fit <- list(
    states = best_states,
    beta_list = best_beta,
    post = best_post,
    trace = trace,
    logLik = best_post$logLik,
    convergence = convergence,
    message = convergence_message
  )
  relabel_fit(fit)
}

viterbi_constant <- function(states, P, data, indices_by_sequence, emission_family = "gaussian", covariance_mode = "full", df = 5) {
  # Viterbi decoding gives one most likely state path for tables and figures.
  # It is not used as if the states were perfectly observed in estimation.
  k <- length(states)
  log_em <- emission_logdens(states, data, emission_family, covariance_mode, df)
  log_P <- log(pmax(P, 1e-12))
  log_pi <- log(stationary_distribution(P))
  decoded <- integer(nrow(data))
  for (idx in indices_by_sequence) {
    n <- length(idx)
    delta <- matrix(NA_real_, nrow = n, ncol = k)
    psi <- matrix(1L, nrow = n, ncol = k)
    delta[1L, ] <- log_pi + log_em[idx[1L], ]
    if (n > 1L) {
      for (tt in 2:n) {
        for (j in seq_len(k)) {
          scores <- delta[tt - 1L, ] + log_P[, j]
          psi[tt, j] <- which.max(scores)
          delta[tt, j] <- max(scores) + log_em[idx[tt], j]
        }
      }
    }
    path <- integer(n)
    path[n] <- which.max(delta[n, ])
    if (n > 1L) {
      for (tt in (n - 1L):1L) path[tt] <- psi[tt + 1L, path[tt + 1L]]
    }
    decoded[idx] <- path
  }
  decoded
}

viterbi_nhmm <- function(states, beta_list, data, indices_by_sequence, design, emission_family = "gaussian", covariance_mode = "full", df = 5) {
  k <- length(states)
  log_em <- emission_logdens(states, data, emission_family, covariance_mode, df)
  first_indices <- unlist(lapply(indices_by_sequence, `[`, 1L), use.names = FALSE)
  avg_P <- Reduce("+", lapply(first_indices, function(idx) transition_probs_for_row(design[idx, ], beta_list))) / length(first_indices)
  log_pi <- log(stationary_distribution(avg_P))
  decoded <- integer(nrow(data))
  for (idx in indices_by_sequence) {
    n <- length(idx)
    delta <- matrix(NA_real_, nrow = n, ncol = k)
    psi <- matrix(1L, nrow = n, ncol = k)
    delta[1L, ] <- log_pi + log_em[idx[1L], ]
    if (n > 1L) {
      for (tt in 2:n) {
        log_P <- log(pmax(transition_probs_for_row(design[idx[tt - 1L], ], beta_list), 1e-12))
        for (j in seq_len(k)) {
          scores <- delta[tt - 1L, ] + log_P[, j]
          psi[tt, j] <- which.max(scores)
          delta[tt, j] <- max(scores) + log_em[idx[tt], j]
        }
      }
    }
    path <- integer(n)
    path[n] <- which.max(delta[n, ])
    if (n > 1L) {
      for (tt in (n - 1L):1L) path[tt] <- psi[tt + 1L, path[tt + 1L]]
    }
    decoded[idx] <- path
  }
  decoded
}

parameter_count_emission <- function(k, full_covariance = TRUE) {
  # alpha_s, alpha_m, lambda, sd_s, sd_m per state, plus rho for full covariance.
  k * 5L + if (full_covariance) k else 0L
}

summarize_hmm_fit <- function(fit, data, model_name, k, transition_parameter_count, emission_family = "gaussian", covariance_mode = "full", df = NA_real_) {
  # Compact model-fit row used by both the primary estimation and robustness
  # scripts. BIC comparisons are meaningful only within the same observation set.
  parameters <- parameter_count_emission(k, full_covariance = covariance_mode == "full") + transition_parameter_count
  data.table(
    model = model_name,
    states = k,
    emission_family = emission_family,
    covariance_mode = covariance_mode,
    student_t_df = df,
    logLik = fit$logLik,
    parameters = parameters,
    observations = nrow(data),
    AIC = -2 * fit$logLik + 2 * parameters,
    BIC = -2 * fit$logLik + log(nrow(data)) * parameters,
    convergence_code = fit$convergence,
    convergence_message = fit$message
  )
}
