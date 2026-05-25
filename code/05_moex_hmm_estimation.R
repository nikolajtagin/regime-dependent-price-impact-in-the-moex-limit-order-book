###############################################################################
# MOEX hidden Markov model estimation
#
# What this script does:
#   1. Loads the primary model-ready panel defined in the shared settings.
#   2. Estimates the primary 3-state bivariate Gaussian HMM with winsorized
#      returns and depth-scaled OFI in the emission mean.
#   3. Estimates the K-state covariate-transition NHMM with multinomial logits.
#   4. Saves compact posterior/Viterbi diagnostics, bootstrap diagnostics, state
#      interpretation tables, transition tables, and paper-style figures.
###############################################################################

library(data.table)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for HMM figures.")
}
if (!requireNamespace("nnet", quietly = TRUE)) {
  stop("Package 'nnet' is required for K-state NHMM transition updates.")
}

code_dir <- local({
  args_all <- commandArgs(trailingOnly = FALSE)
  script_arg <- args_all[startsWith(args_all, "--file=")]
  if (length(script_arg) > 0L) {
    dirname(normalizePath(sub("--file=", "", script_arg[1L], fixed = TRUE), mustWork = TRUE))
  } else {
    current_dir <- normalizePath(getwd(), mustWork = TRUE)
    if (basename(current_dir) == "code") current_dir else file.path(current_dir, "code")
  }
})
source(file.path(code_dir, "01_project_settings_and_functions.R"))

project_root <- normalizePath(file.path(code_dir, ".."), mustWork = TRUE)
cfg <- read_project_settings(project_root)
cfg$tickers <- toupper(split_env_vector(Sys.getenv("MOEX_TICKERS", unset = paste(cfg$tickers, collapse = ","))))
set.seed(cfg$random_seed)

data_file <- combined_data_file(cfg, frequency = cfg$bucket_size_seconds)
if (!file.exists(data_file)) stop("Primary combined dataset not found. Run 03_moex_combine_processed_data.R first: ", data_file)

posterior_dir <- cfg$results_hmm_state_probabilities_root
dirs <- list(
  baseline = cfg$results_hmm_constant_transition_root,
  nhmm = cfg$results_hmm_covariate_transition_root,
  comparison = cfg$results_hmm_model_comparison_root,
  hypotheses = cfg$results_hmm_hypotheses_root,
  state_interpretation = cfg$results_hmm_state_interpretation_root,
  state_sequences = cfg$results_hmm_state_sequence_diagnostics_root
)
for (d in c(posterior_dir, unlist(dirs))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1. LOAD AND PREPARE DATA
# =========================

dt <- fread(data_file)
required_columns <- c(
  "sequence_id", "seccode", "trade_date", "calendar_date", "seconds_from_midnight",
  "spread_tilde", cfg$response_variable, cfg$flow_variable,
  "depth", "depth_tilde", "rv_tilde", "abs_ofi_tilde",
  "spread", "RV_t", "number_of_trades"
)
missing_columns <- setdiff(required_columns, names(dt))
if (length(missing_columns) > 0L) {
  stop("Primary dataset is missing required columns: ", paste(missing_columns, collapse = ", "))
}

setorder(dt, seccode, trade_date, seconds_from_midnight)
dt[, response := get(cfg$response_variable)]
dt[, flow_raw := get(cfg$flow_variable)]
# The emission slope is reported per one-standard-deviation OFI shock. Scaling
# the flow once here makes the lambda estimates directly interpretable.
dt[, flow_scale := safe_sd(flow_raw)]
dt[, flow_scaled := flow_raw / flow_scale]
dt[, open_dummy := as.integer(seconds_from_midnight < cfg$open_end_seconds)]
dt[, close_dummy := as.integer(seconds_from_midnight >= cfg$close_start_seconds)]
dt <- dt[complete.cases(dt[, .(spread_tilde, response, flow_scaled, depth_tilde, rv_tilde, abs_ofi_tilde)])]
if (nrow(dt) < 1000L) stop("Too few complete rows for HMM estimation: ", nrow(dt))

k_primary <- cfg$n_states
state_label_vec <- state_labels(k_primary)
state_name_vec <- state_names(k_primary)
indices_by_sequence <- split(seq_len(nrow(dt)), dt$sequence_id)

message("Loaded rows: ", format(nrow(dt), big.mark = ","))
message("Stock-day sequences: ", length(indices_by_sequence))
message("Primary states: ", k_primary)
message("Response: ", cfg$response_variable)
message("Flow variable: ", cfg$flow_variable)

# =========================
# 2. ONE-STATE BENCHMARK
# =========================

fit_one_state_gaussian <- function(data) {
  # The one-state model is a benchmark, not a candidate final model. It asks how
  # much likelihood is gained by allowing hidden liquidity states at all.
  X <- cbind(1, data$flow_scaled)
  beta <- as.numeric(solve(crossprod(X) + diag(1e-8, 2L), crossprod(X, data$response)))
  alpha_s <- mean(data$spread_tilde)
  resid_s <- data$spread_tilde - alpha_s
  resid_m <- data$response - as.numeric(X %*% beta)
  Sigma <- crossprod(cbind(resid_s, resid_m)) / nrow(data)
  inv_sigma <- solve(Sigma)
  log_det <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)
  resid <- cbind(resid_s, resid_m)
  quad <- rowSums((resid %*% inv_sigma) * resid)
  loglik <- sum(-log(2 * pi) - 0.5 * log_det - 0.5 * quad)
  list(
    logLik = loglik,
    parameters = 6L,
    alpha_s = alpha_s,
    alpha_m = beta[1L],
    lambda = beta[2L],
    sigma_s = sqrt(Sigma[1L, 1L]),
    sigma_m = sqrt(Sigma[2L, 2L]),
    rho_sm = Sigma[1L, 2L] / sqrt(Sigma[1L, 1L] * Sigma[2L, 2L])
  )
}

one_state_fit <- fit_one_state_gaussian(dt)
one_state_table <- data.table(
  model = "1-state Gaussian conditional model",
  states = 1L,
  logLik = one_state_fit$logLik,
  parameters = one_state_fit$parameters,
  observations = nrow(dt),
  AIC = -2 * one_state_fit$logLik + 2 * one_state_fit$parameters,
  BIC = -2 * one_state_fit$logLik + log(nrow(dt)) * one_state_fit$parameters,
  convergence_code = 0L,
  convergence_message = "closed-form conditional Gaussian fit"
)

# =========================
# 3. CONSTANT-TRANSITION HMM
# =========================

message("Estimating constant-transition ", k_primary, "-state HMM.")
# The constant-transition HMM gives the baseline state interpretation and the
# expected durations. It also provides stable starting values for the NHMM.
baseline_fit <- fit_constant_hmm(
  data = dt,
  k = k_primary,
  maxit = cfg$em_maxit,
  tolerance = cfg$em_tol,
  n_starts = cfg$n_multistart,
  seed = cfg$random_seed,
  emission_family = cfg$emission_family,
  covariance_mode = "full",
  df = cfg$student_t_df
)
write_table(baseline_fit$start_convergence, file.path(dirs$baseline, "hmm_baseline_start_convergence.csv"))
write_table(baseline_fit$trace, file.path(dirs$baseline, "hmm_baseline_convergence_trace.csv"))
baseline_decoded <- viterbi_constant(baseline_fit$states, baseline_fit$P, dt, indices_by_sequence)
dt[, decoded_state_constant := baseline_decoded]

baseline_fit_table <- summarize_hmm_fit(
  baseline_fit,
  dt,
  sprintf("%d-state Gaussian HMM, constant transitions", k_primary),
  k_primary,
  transition_parameter_count = k_primary * (k_primary - 1L),
  emission_family = cfg$emission_family,
  covariance_mode = "full"
)
write_table(baseline_fit_table, file.path(dirs$baseline, "hmm_baseline_model_fit.csv"))

P_constant_labelled <- baseline_fit$P
dimnames(P_constant_labelled) <- list(from_state = state_label_vec, to_state = state_label_vec)
transition_constant <- as.data.table(as.table(P_constant_labelled))
setnames(transition_constant, c("from_state", "to_state", "probability"))
transition_constant[, `:=`(
  from_state = as.character(from_state),
  to_state = as.character(to_state)
)]
transition_constant[, `:=`(
  from_order = match(from_state, state_label_vec),
  to_order = match(to_state, state_label_vec)
)]
setorder(transition_constant, from_order, to_order)
transition_constant[, `:=`(from_order = NULL, to_order = NULL)]
write_table(transition_constant, file.path(dirs$baseline, "hmm_baseline_transition_matrix.csv"))

duration_table <- data.table(
  state = state_label_vec,
  p_ii = diag(baseline_fit$P),
  expected_duration_buckets = 1 / pmax(1 - diag(baseline_fit$P), 1e-8),
  expected_duration_minutes = (cfg$bucket_size_seconds / 60) / pmax(1 - diag(baseline_fit$P), 1e-8)
)
write_table(duration_table, file.path(dirs$baseline, "hmm_baseline_expected_durations.csv"))

# =========================
# 4. COVARIATE-TRANSITION NHMM
# =========================

transition_formula <- ~ depth_tilde + rv_tilde + abs_ofi_tilde + open_dummy + close_dummy
transition_design <- model.matrix(transition_formula, data = dt)

message("Estimating covariate-transition NHMM.")
# In the NHMM, emissions still describe how each state looks, while transition
# logits describe which observed conditions make moves between states more or
# less likely.
nhmm_fit <- fit_nhmm(
  data = dt,
  k = k_primary,
  baseline_fit = baseline_fit,
  design = transition_design,
  maxit = cfg$em_maxit,
  tolerance = cfg$em_tol,
  seed = cfg$random_seed,
  emission_family = cfg$emission_family,
  covariance_mode = "full",
  df = cfg$student_t_df
)
write_table(nhmm_fit$trace, file.path(dirs$nhmm, "nhmm_convergence_trace_primary.csv"))
nhmm_decoded <- viterbi_nhmm(nhmm_fit$states, nhmm_fit$beta_list, dt, indices_by_sequence, transition_design)
dt[, decoded_state := nhmm_decoded]

nhmm_parameter_count <- parameter_count_emission(k_primary, full_covariance = TRUE) +
  k_primary * (k_primary - 1L) * ncol(transition_design)
nhmm_fit_table <- summarize_hmm_fit(
  nhmm_fit,
  dt,
  sprintf("%d-state Gaussian HMM, covariate transitions", k_primary),
  k_primary,
  transition_parameter_count = k_primary * (k_primary - 1L) * ncol(transition_design),
  emission_family = cfg$emission_family,
  covariance_mode = "full"
)
nhmm_fit_table[, `:=`(
  response_variable = cfg$response_variable,
  flow_variable = cfg$flow_variable,
  flow_scale = dt$flow_scale[1L]
)]
write_table(nhmm_fit_table, file.path(dirs$nhmm, "nhmm_model_fit.csv"))

transition_coefficients <- rbindlist(lapply(seq_len(k_primary), function(origin) {
  beta <- nhmm_fit$beta_list[[origin]]
  out <- melt(
    as.data.table(beta, keep.rownames = "to_state_index"),
    id.vars = "to_state_index",
    variable.name = "term",
    value.name = "estimate"
  )
  out[, `:=`(
    origin_state = state_label_vec[origin],
    to_state = state_label_vec[as.integer(sub("to_", "", to_state_index))]
  )]
  out[, .(origin_state, to_state, term, estimate)]
}), use.names = TRUE, fill = TRUE)
write_table(transition_coefficients, file.path(dirs$nhmm, "nhmm_transition_coefficients.csv"))

transition_rows <- unlist(lapply(indices_by_sequence, function(idx) if (length(idx) <= 1L) integer(0) else idx[-length(idx)]), use.names = FALSE)
average_transition_matrix <- function(beta_list, design_rows) {
  # A time-varying transition model has one matrix for every row. This function
  # averages those matrices over the observed covariate distribution to produce
  # a compact table comparable to the constant-transition matrix.
  k <- length(beta_list)
  out <- matrix(NA_real_, nrow = k, ncol = k)
  for (origin in seq_len(k)) {
    eta <- design_rows %*% t(beta_list[[origin]])
    eta <- sweep(eta, 1L, apply(eta, 1L, max), "-")
    prob <- exp(eta)
    prob <- prob / rowSums(prob)
    out[origin, ] <- colMeans(prob, na.rm = TRUE)
  }
  out
}
avg_P <- average_transition_matrix(nhmm_fit$beta_list, transition_design[transition_rows, , drop = FALSE])
dimnames(avg_P) <- list(from_state = state_label_vec, to_state = state_label_vec)
avg_transition <- as.data.table(as.table(avg_P))
setnames(avg_transition, c("from_state", "to_state", "probability"))
avg_transition[, `:=`(
  from_state = as.character(from_state),
  to_state = as.character(to_state)
)]
avg_transition[, `:=`(
  from_order = match(from_state, state_label_vec),
  to_order = match(to_state, state_label_vec)
)]
setorder(avg_transition, from_order, to_order)
avg_transition[, `:=`(from_order = NULL, to_order = NULL)]
write_table(avg_transition, file.path(dirs$nhmm, "nhmm_average_implied_transition_matrix.csv"))

# =========================
# 5. POSTERIORS, STATE INTERPRETATION, AND RESIDUALS
# =========================

post_cols <- paste0("post_", seq_len(k_primary))
dt[, (post_cols) := as.data.table(nhmm_fit$post$gamma)]
dt[, max_posterior := do.call(pmax, .SD), .SDcols = post_cols]
dt[, stressed_probability := get(post_cols[k_primary])]

# Posterior probabilities and the decoded path are kept in memory for figures,
# state summaries, and diagnostics. They are not written as full CSV files
# because those files are large and do not add information beyond the saved
# compact figures/tables needed for the paper.

state_interpretation <- rbindlist(lapply(seq_len(k_primary), function(state) {
  # Posterior-weighted means are used for state summaries. This avoids treating
  # uncertain observations as if their Viterbi labels were known with certainty.
  w <- nhmm_fit$post$gamma[, state]
  p <- nhmm_fit$states[[state]]
  data.table(
    state = state,
    state_label = state_label_vec[state],
    posterior_share = mean(w),
    alpha_s = p$alpha_s,
    alpha_m = p$alpha_m,
    lambda_scaled_flow = p$lambda,
    sigma_s = p$sd_s,
    sigma_m = p$sd_m,
    rho_sm = p$rho,
    mean_spread_bp = weighted.mean(dt$spread * 10000 / dt$midprice, w, na.rm = TRUE),
    mean_depth = weighted.mean(dt$depth, w, na.rm = TRUE),
    mean_RV = weighted.mean(dt$RV_t, w, na.rm = TRUE),
    mean_abs_ofi = weighted.mean(abs(dt$event_level_ofi), w, na.rm = TRUE)
  )
}), use.names = TRUE, fill = TRUE)
write_table(state_interpretation, file.path(dirs$baseline, "hmm_baseline_state_interpretation.csv"))
write_table(state_interpretation, file.path(dirs$nhmm, "nhmm_state_interpretation.csv"))

dt[, e_s := NA_real_]
dt[, e_m := NA_real_]
for (state in seq_len(k_primary)) {
  # Residuals are standardized with the parameters of the decoded state. They
  # are used only for diagnostics of the conditional Gaussian approximation.
  p <- nhmm_fit$states[[state]]
  dt[decoded_state == state, `:=`(
    e_s = (spread_tilde - p$alpha_s) / p$sd_s,
    e_m = (response - p$alpha_m - p$lambda * flow_scaled) / p$sd_m
  )]
}
residual_kurtosis <- dt[, .(
  kurtosis_r = kurtosis_manual(e_m),
  kurtosis_s = kurtosis_manual(e_s),
  n = .N
), by = .(decoded_state, state_label = state_label_vec[decoded_state])]
write_table(residual_kurtosis, file.path(dirs$nhmm, "state_conditional_residual_kurtosis.csv"))

classification_uncertainty <- data.table(
  average_max_posterior = mean(dt$max_posterior, na.rm = TRUE),
  share_max_posterior_below_070 = mean(dt$max_posterior < 0.70, na.rm = TRUE),
  share_max_posterior_below_080 = mean(dt$max_posterior < 0.80, na.rm = TRUE)
)
write_table(classification_uncertainty, file.path(dirs$nhmm, "classification_uncertainty.csv"))

per_stock <- dt[, .(
  # The primary HMM is pooled, so lambda is common across stocks. This table
  # shows how much the decoded Stressed share and Stressed-state spread still
  # vary across names.
  stressed_share = mean(decoded_state == k_primary),
  lambda_stressed = nhmm_fit$states[[k_primary]]$lambda,
  lambda_calm = nhmm_fit$states[[1L]]$lambda,
  avg_spread_bp_stressed = mean(spread[decoded_state == k_primary] * 10000 / midprice[decoded_state == k_primary], na.rm = TRUE),
  avg_depth_stressed = mean(depth[decoded_state == k_primary], na.rm = TRUE),
  avg_max_posterior = mean(max_posterior, na.rm = TRUE)
), by = seccode]
write_table(per_stock, file.path(dirs$state_sequences, "per_stock_state_summary.csv"))

per_stock[, label_hjust := 0.5]
per_stock[, label_vjust := -0.75]
per_stock[seccode == "TATN", `:=`(label_hjust = 0.05, label_vjust = -1.15)]
per_stock[seccode == "YDEX", `:=`(label_hjust = 0.95, label_vjust = 1.45)]
per_stock_scatter <- ggplot2::ggplot(per_stock, ggplot2::aes(
  x = stressed_share,
  y = avg_spread_bp_stressed,
  label = seccode
)) +
  ggplot2::geom_point(size = 2.4, colour = paper_colours$calm) +
  ggplot2::geom_text(
    ggplot2::aes(hjust = label_hjust, vjust = label_vjust),
    size = 3.1,
    colour = paper_colours$dark
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(x = "Decoded Stressed-state share", y = "Mean spread in Stressed (bp)") +
  paper_theme()
save_paper_plot(
  per_stock_scatter,
  file.path(dirs$state_sequences, "per_stock_stressed_share_spread.png"),
  width = 8.8,
  height = 5.0
)

# =========================
# 6. BOOTSTRAP AND MARGINAL EFFECTS
# =========================

transition_probability_vector <- function(beta_list, design, origin, destination) {
  # Convert a fitted origin-specific multinomial logit into a vector of
  # probabilities for one destination state.
  eta <- design %*% t(beta_list[[origin]])
  eta <- sweep(eta, 1L, apply(eta, 1L, max), "-")
  prob <- exp(eta)
  prob <- prob / rowSums(prob)
  prob[, destination]
}

compute_transition_ame <- function(fit, design, covariate_name, origin, destination) {
  # Average marginal effect: increase one covariate by one unit for every row,
  # recompute transition probabilities, and average the probability change.
  if (!covariate_name %in% colnames(design)) {
    return(NA_real_)
  }
  design_hi <- design
  design_hi[, covariate_name] <- design_hi[, covariate_name] + 1
  p0 <- transition_probability_vector(fit$beta_list, design, origin, destination)
  p1 <- transition_probability_vector(fit$beta_list, design_hi, origin, destination)
  mean(p1 - p0, na.rm = TRUE)
}

ame_terms <- c("depth_tilde", "rv_tilde", "abs_ofi_tilde", "open_dummy", "close_dummy")
ame_table <- rbindlist(lapply(ame_terms, function(term) {
  data.table(
    covariate = term,
    entry_to_stressed_AME = compute_transition_ame(nhmm_fit, transition_design, term, origin = 1L, destination = k_primary),
    exit_to_calm_AME = compute_transition_ame(nhmm_fit, transition_design, term, origin = k_primary, destination = 1L)
  )
}), use.names = TRUE, fill = TRUE)

n_bootstrap_reps <- as.integer(Sys.getenv("MOEX_BOOTSTRAP_REPS", unset = as.character(cfg$n_bootstrap_reps)))
bootstrap_maxit <- as.integer(Sys.getenv("MOEX_BOOTSTRAP_MAXIT", unset = as.character(cfg$bootstrap_maxit)))
bootstrap_results <- list()
if (n_bootstrap_reps > 0L) {
  run_one_bootstrap <- function(rep_id) {
    # Resample whole calendar days. This preserves within-day serial dependence
    # and same-day cross-stock co-movement better than row-by-row resampling.
    set.seed(cfg$random_seed + rep_id)
    sampled_days <- sample(unique(dt$trade_date), replace = TRUE)
    boot_parts <- lapply(seq_along(sampled_days), function(i) {
      x <- copy(dt[trade_date == sampled_days[i]])
      x[, sequence_id := paste(sequence_id, "boot", i, sep = "_")]
      x
    })
    boot_dt <- rbindlist(boot_parts, use.names = TRUE, fill = TRUE)
    boot_design <- model.matrix(transition_formula, data = boot_dt)
    boot_indices <- split(seq_len(nrow(boot_dt)), boot_dt$sequence_id)

    # Reuse the production constant-HMM parameters as the NHMM starting point,
    # but recompute posteriors on the resampled day-cluster panel.
    boot_baseline <- baseline_fit
    boot_baseline$post <- forward_backward_constant(boot_baseline$states, boot_baseline$P, boot_dt, boot_indices)
    boot_fit <- fit_nhmm(
      data = boot_dt,
      k = k_primary,
      baseline_fit = boot_baseline,
      design = boot_design,
      maxit = bootstrap_maxit,
      tolerance = cfg$em_tol,
      seed = cfg$random_seed + rep_id,
      emission_family = cfg$emission_family,
      covariance_mode = "full",
      df = cfg$student_t_df
    )
    rbindlist(lapply(ame_terms, function(term) {
      data.table(
        replication = rep_id,
        covariate = term,
        entry_to_stressed_AME = compute_transition_ame(boot_fit, boot_design, term, origin = 1L, destination = k_primary),
        exit_to_calm_AME = compute_transition_ame(boot_fit, boot_design, term, origin = k_primary, destination = 1L),
        convergence_code = boot_fit$convergence
      )
    }), use.names = TRUE, fill = TRUE)
  }

  message("Running day-cluster bootstrap with ", n_bootstrap_reps, " replications.")
  if (requireNamespace("doParallel", quietly = TRUE) && requireNamespace("foreach", quietly = TRUE)) {
    library(foreach)
    library(doParallel)
    cores <- as.integer(Sys.getenv("MOEX_BOOTSTRAP_WORKERS", unset = as.character(min(max(parallel::detectCores() - 2L, 1L), 4L))))
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    bootstrap_results <- foreach::foreach(
      b = seq_len(n_bootstrap_reps),
      .packages = c("data.table", "nnet"),
      .export = c(
        "dt", "transition_formula", "baseline_fit", "k_primary", "cfg", "ame_terms",
        "run_one_bootstrap",
        "fit_nhmm", "forward_backward_constant", "compute_transition_ame",
        "transition_probability_vector", "transition_probs_for_row",
        "fit_multinom_origin", "update_nhmm_transitions",
        "update_emissions", "forward_backward_nhmm", "emission_logdens",
        "student_t_weights", "safe_sd", "make_start_labels", "start_from_labels",
        "stationary_distribution", "logsumexp", "softmax_row", "relabel_fit",
        "parameter_count_emission"
      ),
      .errorhandling = "pass"
    ) %dopar% {
      run_one_bootstrap(b)
    }
    parallel::stopCluster(cl)
  } else {
    warning("doParallel is unavailable; running bootstrap sequentially.")
    bootstrap_results <- lapply(seq_len(n_bootstrap_reps), run_one_bootstrap)
  }

  boot_dt <- rbindlist(bootstrap_results[!vapply(bootstrap_results, inherits, logical(1L), what = "error")], use.names = TRUE, fill = TRUE)
  write_table(boot_dt, file.path(dirs$nhmm, "nhmm_transition_ame_bootstrap_draws.csv"))
  ame_ci <- boot_dt[, .(
    entry_ci_lower_025 = quantile(entry_to_stressed_AME, 0.025, na.rm = TRUE, names = FALSE),
    entry_ci_upper_975 = quantile(entry_to_stressed_AME, 0.975, na.rm = TRUE, names = FALSE),
    exit_ci_lower_025 = quantile(exit_to_calm_AME, 0.025, na.rm = TRUE, names = FALSE),
    exit_ci_upper_975 = quantile(exit_to_calm_AME, 0.975, na.rm = TRUE, names = FALSE),
    successful_replications = .N
  ), by = covariate]
  ame_table <- merge(ame_table, ame_ci, by = "covariate", all.x = TRUE)
}
write_table(ame_table, file.path(dirs$nhmm, "nhmm_transition_average_marginal_effects.csv"))

# =========================
# 7. MODEL COMPARISON AND HYPOTHESIS TABLE
# =========================

model_comparison <- rbindlist(list(one_state_table, baseline_fit_table, nhmm_fit_table), use.names = TRUE, fill = TRUE)
setorder(model_comparison, BIC)
# BIC comparisons here are safe because all three rows use the same primary
# 30-second sample and the same response/flow variables.
model_comparison[, delta_BIC := BIC - min(BIC, na.rm = TRUE)]
write_table(model_comparison, file.path(dirs$comparison, "hmm_model_comparison.csv"))

hypothesis_summary <- data.table(
  hypothesis = c(
    "H1 entry into stressed liquidity",
    "H2 exit from stressed liquidity",
    "H3 directional session-boundary transition effects",
    "H4 value added of covariate-driven transitions"
  ),
  diagnostic = c(
    "Average marginal effects in the calm-to-stressed transition equation.",
    "Average marginal effects in the stressed-to-calm transition equation.",
    "Average marginal effects of open and close dummies for entry and recovery.",
    "BIC comparison of constant-transition HMM and covariate-transition NHMM, including per-observation difference."
  ),
  result = c(
    paste(sprintf("%s: %.4f", ame_table$covariate, ame_table$entry_to_stressed_AME), collapse = "; "),
    paste(sprintf("%s: %.4f", ame_table$covariate, ame_table$exit_to_calm_AME), collapse = "; "),
    paste(sprintf("%s: entry %.4f, exit %.4f", ame_table[covariate %in% c("open_dummy", "close_dummy")]$covariate, ame_table[covariate %in% c("open_dummy", "close_dummy")]$entry_to_stressed_AME, ame_table[covariate %in% c("open_dummy", "close_dummy")]$exit_to_calm_AME), collapse = "; "),
    sprintf("Constant-transition BIC %.2f; covariate-transition BIC %.2f; difference %.2f; per observation %.5f.", baseline_fit_table$BIC, nhmm_fit_table$BIC, nhmm_fit_table$BIC - baseline_fit_table$BIC, (nhmm_fit_table$BIC - baseline_fit_table$BIC) / nrow(dt))
  )
)
write_table(hypothesis_summary, file.path(dirs$hypotheses, "hypothesis_summary.csv"))

# =========================
# 8. FIGURES
# =========================

plot_transition_heatmap <- function(matrix_dt, output_dir, file_name) {
  p <- ggplot2::ggplot(matrix_dt, ggplot2::aes(x = to_state, y = from_state, fill = probability)) +
    ggplot2::geom_tile(color = "white", linewidth = 1.0) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", probability), color = probability > 0.55), fontface = "bold", size = 4.8) +
    ggplot2::scale_fill_gradientn(colours = c("#D1D1D1", "#8E6BF2", "#1D22FF"), limits = c(0, 1)) +
    ggplot2::scale_color_manual(values = c("FALSE" = paper_colours$dark, "TRUE" = "white"), guide = "none") +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "Next state", y = "Current state", fill = "Probability") +
    paper_theme() +
    ggplot2::theme(legend.position = "right", panel.grid = ggplot2::element_blank())
  save_paper_plot(p, file.path(output_dir, file_name), width = 6.8, height = 5.2)
}
plot_transition_heatmap(transition_constant, dirs$baseline, "baseline_transition_matrix_heatmap.png")
plot_transition_heatmap(avg_transition, dirs$nhmm, "nhmm_average_transition_matrix_heatmap.png")

profile_dt <- dt[, .(mean_stressed_probability = mean(stressed_probability, na.rm = TRUE)), by = .(
  profile_bin = floor(seconds_from_midnight / (15 * 60)) * 15 * 60
)]
profile_dt[, hour := profile_bin / 3600]
profile_plot <- ggplot2::ggplot(profile_dt, ggplot2::aes(x = hour, y = mean_stressed_probability)) +
  ggplot2::geom_area(fill = paper_colours$light_red, alpha = 0.80) +
  ggplot2::geom_line(color = paper_colours$stressed, linewidth = 0.80) +
  ggplot2::geom_point(color = paper_colours$stressed, size = 1.4) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(x = "Hour of day", y = "Mean smoothed probability") +
  paper_theme()
save_paper_plot(profile_plot, file.path(posterior_dir, "stressed_probability_intraday_profile.png"), width = 8.8, height = 5.0)

representative_sequence <- dt[, .N, by = sequence_id][order(-N), sequence_id[1L]]
representative_dt <- dt[sequence_id == representative_sequence]
representative_dt[, hour := seconds_from_midnight / 3600]
representative_dt[, rolling_price_impact := {
  # Right-aligned local slope of return on scaled OFI. It is a visual diagnostic
  # for one representative day, not an additional estimated HMM parameter.
  window <- max(10L, as.integer(round(30 * 60 / cfg$bucket_size_seconds)))
  num <- frollmean(response * flow_scaled, n = window, align = "right", na.rm = TRUE)
  den <- frollmean(flow_scaled^2, n = window, align = "right", na.rm = TRUE)
  fifelse(is.finite(den) & den > 1e-10, num / den, NA_real_)
}]

prob_dt <- representative_dt[, .(
  hour,
  Calm = post_1,
  Normal = if (k_primary >= 3L) post_2 else NA_real_,
  Stressed = get(post_cols[k_primary])
)]
prob_long <- melt(prob_dt, id.vars = "hour", variable.name = "state", value.name = "probability")
prob_long <- prob_long[is.finite(probability)]
prob_long[, panel := "State probabilities"]
var_dt <- representative_dt[, .(
  hour,
  `Transformed spread` = spread_tilde,
  `Rolling price impact` = rolling_price_impact
)]
var_long <- melt(var_dt, id.vars = "hour", variable.name = "panel", value.name = "value")
panel_levels <- c("State probabilities", "Transformed spread", "Rolling price impact")
prob_long[, panel := factor(panel, levels = panel_levels)]
var_long[, panel := factor(panel, levels = panel_levels)]
probability_panel <- ggplot2::ggplot() +
  ggplot2::geom_area(data = prob_long, ggplot2::aes(x = hour, y = probability, fill = state), position = "stack", alpha = 0.98, linewidth = 0) +
  ggplot2::geom_line(data = var_long[is.finite(value)], ggplot2::aes(x = hour, y = value, colour = panel), linewidth = 0.42, show.legend = FALSE) +
  ggplot2::facet_grid(panel ~ ., scales = "free_y", switch = "y") +
  ggplot2::scale_fill_manual(values = state_palette(k_primary)) +
  ggplot2::scale_colour_manual(values = c("Transformed spread" = paper_colours$calm, "Rolling price impact" = paper_colours$purple)) +
  ggplot2::labs(x = "Hour of day", y = NULL) +
  paper_theme() +
  ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", hjust = 0))
save_paper_plot(probability_panel, file.path(posterior_dir, "nhmm_probability_liquidity_panel.png"), width = 9.2, height = 6.4)

viterbi_plot_dt <- representative_dt[, .(hour, decoded_state_label = state_label_vec[decoded_state])]
viterbi_plot_dt[, decoded_state_label := factor(decoded_state_label, levels = state_label_vec)]
viterbi_plot <- ggplot2::ggplot(viterbi_plot_dt, ggplot2::aes(x = hour, y = decoded_state_label, colour = decoded_state_label)) +
  ggplot2::geom_step(linewidth = 0.75) +
  ggplot2::scale_colour_manual(values = state_palette(k_primary)) +
  ggplot2::labs(x = "Hour of day", y = "Decoded state") +
  paper_theme() +
  ggplot2::theme(legend.position = "none")
save_paper_plot(viterbi_plot, file.path(posterior_dir, "representative_viterbi_path.png"), width = 8.8, height = 4.8)

state_summary_long <- melt(
  state_interpretation[, .(state_label, mean_spread_bp, lambda_scaled_flow, sigma_m, rho_sm, mean_RV)],
  id.vars = "state_label",
  variable.name = "metric",
  value.name = "value"
)
state_summary_long[, z_value := as.numeric(scale(value)), by = metric]
# Figure 2 compares metrics with very different units. Standardizing within each
# metric makes the relative state ordering visible in one plot.
metric_labels <- c(
  mean_spread_bp = "Mean spread\n(bp)",
  lambda_scaled_flow = "OFI-return\nslope lambda",
  sigma_m = "Within-state\nreturn s.d.",
  rho_sm = "Spread-return\ncorrelation",
  mean_RV = "Realized\nvariance"
)
state_summary_long[, metric_label := factor(metric_labels[metric], levels = metric_labels)]
state_plot <- ggplot2::ggplot(state_summary_long, ggplot2::aes(x = metric_label, y = z_value, fill = state_label)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.62, colour = paper_colours$dark, linewidth = 0.20) +
  ggplot2::scale_fill_manual(values = state_palette(k_primary)) +
  ggplot2::labs(x = NULL, y = "Standardized value") +
  paper_theme()
save_paper_plot(state_plot, file.path(dirs$state_interpretation, "state_interpretation.png"), width = 8.8, height = 5.0)

qq_dt <- dt[is.finite(e_m), .(e_m, decoded_state_label = state_label_vec[decoded_state])]
if (nrow(qq_dt) > 75000L) {
  set.seed(cfg$random_seed)
  qq_dt <- qq_dt[sample.int(.N, 75000L)]
}
qq_dt[, decoded_state_label := factor(decoded_state_label, levels = state_label_vec)]
qq_plot <- ggplot2::ggplot(qq_dt, ggplot2::aes(sample = e_m)) +
  ggplot2::stat_qq(alpha = 0.28, size = 0.55, colour = paper_colours$calm) +
  ggplot2::stat_qq_line(colour = paper_colours$stressed, linewidth = 0.75) +
  ggplot2::facet_wrap(~decoded_state_label, nrow = 1, scales = "free") +
  ggplot2::labs(x = "Theoretical normal quantiles", y = "Standardized return residual") +
  paper_theme()
save_paper_plot(qq_plot, file.path(dirs$state_interpretation, "state_conditional_return_residual_qq.png"), width = 10.2, height = 4.2)

post_hist <- ggplot2::ggplot(dt, ggplot2::aes(x = max_posterior)) +
  ggplot2::geom_histogram(bins = 60, fill = paper_colours$calm, colour = "white", linewidth = 0.2) +
  ggplot2::labs(x = "Maximum smoothed posterior probability", y = "Count") +
  paper_theme()
save_paper_plot(post_hist, file.path(dirs$state_interpretation, "maximum_posterior_histogram.png"), width = 8.8, height = 5.0)

runs <- dt[order(seccode, trade_date, seconds_from_midnight),
  {
    r <- rle(decoded_state)
    .(lengths = r$lengths, state = r$values)
  },
  by = .(seccode, trade_date)
]
runs[, state_label := state_label_vec[state]]
# Very long runs exist but are rare. Grouping them into a 30+ bin makes the
# short-duration structure readable without creating a second histogram version.
runs[, lengths_display := pmin(lengths, 30L)]
dwell_plot <- ggplot2::ggplot(runs, ggplot2::aes(x = lengths_display)) +
  ggplot2::geom_histogram(
    ggplot2::aes(y = ggplot2::after_stat(density)),
    binwidth = 1,
    boundary = 0.5,
    fill = "#D9DEE7",
    colour = "white"
  ) +
  ggplot2::facet_wrap(~state_label, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20, 25, 30),
    labels = c("1", "5", "10", "15", "20", "25", "30+")
  ) +
  ggplot2::labs(x = sprintf("Dwell time (%s-second buckets; 30+ bin)", cfg$bucket_size_seconds), y = "Density") +
  paper_theme()
save_paper_plot(
  dwell_plot,
  file.path(dirs$state_interpretation, "dwell_time_histograms.png"),
  width = 10.2,
  height = 4.8
)
runs[, lengths_display := NULL]
write_table(runs, file.path(dirs$state_sequences, "dwell_times.csv"))

ame_long <- melt(ame_table, id.vars = "covariate", measure.vars = c("entry_to_stressed_AME", "exit_to_calm_AME"), variable.name = "transition", value.name = "AME")
ame_long[, transition := factor(transition, levels = c("entry_to_stressed_AME", "exit_to_calm_AME"), labels = c("Entry to stressed", "Exit to calm"))]
ame_plot <- ggplot2::ggplot(ame_long, ggplot2::aes(x = covariate, y = AME, fill = transition)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.62, colour = paper_colours$dark, linewidth = 0.2) +
  ggplot2::scale_fill_manual(values = c("Entry to stressed" = paper_colours$stressed, "Exit to calm" = paper_colours$calm)) +
  ggplot2::labs(x = NULL, y = "Average marginal effect") +
  paper_theme()
save_paper_plot(ame_plot, file.path(dirs$nhmm, "transition_average_marginal_effects.png"), width = 8.8, height = 5.0)

model_plot_dt <- copy(model_comparison)
model_plot_dt[, delta_BIC_plot := log10(delta_BIC + 1)]
model_plot_dt[, model := factor(model, levels = rev(model))]
model_plot <- ggplot2::ggplot(model_plot_dt, ggplot2::aes(x = delta_BIC_plot, y = model, fill = states == k_primary)) +
  ggplot2::geom_col(width = 0.62, colour = paper_colours$dark, linewidth = 0.25) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", delta_BIC)), hjust = -0.10, size = 3.3, colour = paper_colours$dark) +
  ggplot2::scale_fill_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = paper_colours$stressed), guide = "none") +
  ggplot2::coord_cartesian(xlim = c(0, max(model_plot_dt$delta_BIC_plot, na.rm = TRUE) * 1.18 + 0.05)) +
  ggplot2::labs(x = expression(log[10] * (1 + Delta * "BIC")), y = NULL) +
  paper_theme()
save_paper_plot(model_plot, file.path(dirs$comparison, "main_model_delta_bic.png"), width = 8.8, height = 4.8)

message("Done.")
message("HMM results folders: ", cfg$results_root)
