###############################################################################
# MOEX HMM robustness models
#
# What this script does:
#   1. Loads the model-ready datasets created by 03_moex_combine_processed_data.R.
#   2. Estimates exactly the five robustness specifications in the final design.
#   3. Saves a compact robustness table and a paper-style BIC comparison figure.
#
# The script does not touch raw order-log files.
###############################################################################

library(data.table)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for robustness figures.")
}

code_dir <- local({
  # Allow the script to be run either with Rscript code/06_...R or from inside
  # the code folder. All relative paths are rebuilt from this detected location.
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
# Keep the same optional ticker override as the main scripts. This is useful for
# quick checks on a smaller ticker set, while the paper values use all 15 stocks.
cfg$tickers <- toupper(split_env_vector(Sys.getenv("MOEX_TICKERS", unset = paste(cfg$tickers, collapse = ","))))
set.seed(cfg$random_seed)

robustness_root <- cfg$results_hmm_robustness_root
robustness_figures_root <- cfg$results_hmm_robustness_root
for (d in c(robustness_root, robustness_figures_root)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Robustness models are intentionally run with the same EM tolerance and the
# same four-start design as the primary model. This makes the comparison cleaner
# than using a fast exploratory grid.
robustness_em_maxit <- as.integer(Sys.getenv("MOEX_ROBUSTNESS_EM_MAXIT", unset = "200"))
robustness_em_tolerance <- as.numeric(Sys.getenv("MOEX_ROBUSTNESS_EM_TOLERANCE", unset = as.character(cfg$em_tol)))
robustness_n_starts <- as.integer(Sys.getenv("MOEX_ROBUSTNESS_N_STARTS", unset = as.character(cfg$n_multistart)))

done_flag <- file.path(robustness_root, "hmm_robustness_model_comparison.csv")
if (file.exists(done_flag)) {
  # Robustness models are slow. If the final comparison table is already present,
  # rerunning the script should be a no-op unless the user deletes the table.
  message("Robustness results already complete, skipping re-estimation: ", done_flag)
  quit(save = "no", status = 0L)
}

load_model_data <- function(path, response_col, flow_col) {
  # Each robustness check may change the sampling frequency or order-flow proxy,
  # but the estimator expects the same canonical column names: response and
  # flow_scaled. This loader performs that standardization in one place.
  if (!file.exists(path)) stop("Required robustness dataset not found: ", path)
  x <- fread(path)
  required <- c(
    "sequence_id", "seccode", "trade_date", "seconds_from_midnight",
    "spread_tilde", response_col, flow_col, "depth", "depth_tilde", "rv_tilde"
  )
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) stop("Dataset is missing required robustness columns: ", paste(missing, collapse = ", "))
  setorder(x, seccode, trade_date, seconds_from_midnight)
  x[, response := get(response_col)]
  x[, flow_raw := get(flow_col)]
  # Scaling is recomputed inside each robustness sample so that the reported
  # lambda is always "basis points per one-s.d. order-flow shock" for that sample.
  x[, flow_scale := safe_sd(flow_raw)]
  x[, flow_scaled := flow_raw / flow_scale]
  x <- x[complete.cases(x[, .(spread_tilde, response, flow_scaled, depth_tilde, rv_tilde)])]
  if (nrow(x) < 1000L) stop("Too few complete rows for robustness model: ", path)
  x
}

fit_robustness_spec <- function(spec) {
  # A robustness specification is deliberately small: it names the dataset,
  # response, order-flow proxy, number of states, and emission family. The actual
  # HMM code is shared with the primary estimation script through the settings
  # and functions file.
  message("Fitting robustness model: ", spec$name)
  data <- load_model_data(spec$data_file, spec$response_variable, spec$flow_variable)
  fit <- fit_constant_hmm(
    data = data,
    k = spec$n_states,
    maxit = robustness_em_maxit,
    tolerance = robustness_em_tolerance,
    n_starts = robustness_n_starts,
    seed = cfg$random_seed + spec$id,
    emission_family = spec$emission_family,
    covariance_mode = "full",
    df = spec$student_t_df
  )
  decoded <- viterbi_constant(
    fit$states,
    fit$P,
    data,
    split(seq_len(nrow(data)), data$sequence_id),
    emission_family = spec$emission_family,
    covariance_mode = "full",
    df = spec$student_t_df
  )
  # Decoded shares are not used for estimation. They make the comparison table
  # readable by showing whether a specification still creates economically
  # recognizable low-, middle-, and high-cost regimes.
  data[, decoded_state := decoded]
  state_share <- data[, .N, by = decoded_state][order(decoded_state)]
  state_share[, share := N / sum(N)]
  lambda_by_state <- vapply(fit$states, function(s) s$lambda, numeric(1L))
  top_lambda_state <- which.max(abs(lambda_by_state))
  summary_row <- summarize_hmm_fit(
    fit,
    data,
    spec$name,
    spec$n_states,
    transition_parameter_count = spec$n_states * (spec$n_states - 1L),
    emission_family = spec$emission_family,
    covariance_mode = "full",
    df = if (identical(spec$emission_family, "student_t")) spec$student_t_df else NA_real_
  )
  summary_row[, `:=`(
    robustness_id = spec$id,
    robustness_check = spec$check,
    frequency_seconds = spec$frequency_seconds,
    response_variable = spec$response_variable,
    flow_variable = spec$flow_variable,
    trade_buckets_only = spec$trade_buckets_only,
    state_shares = paste(sprintf("S%d=%.3f", state_share$decoded_state, state_share$share), collapse = "; "),
    top_lambda_state = top_lambda_state,
    conclusion = spec$conclusion
  )]
  list(fit = fit, summary = summary_row)
}

primary_file <- combined_data_file(cfg, frequency = cfg$bucket_size_seconds)
trade_5s_file <- combined_data_file(cfg, frequency = cfg$reconstruction_bucket_size_seconds, trade_buckets_only = TRUE)
freq_60_file <- combined_data_file(cfg, frequency = 60L)

# These five checks are the final set used in the paper. They are intentionally
# limited to changes that answer a clear question: number of states, tail shape,
# very short sampling, order-flow proxy, and a coarser 60-second grid.
robustness_specs <- list(
  list(
    id = 1L,
    name = "2-state same-emission Gaussian HMM",
    check = "K=3 identification",
    data_file = primary_file,
    frequency_seconds = cfg$bucket_size_seconds,
    n_states = 2L,
    response_variable = cfg$response_variable,
    flow_variable = cfg$flow_variable,
    emission_family = "gaussian",
    student_t_df = NA_real_,
    trade_buckets_only = FALSE,
    conclusion = "Checks whether the three-state structure is needed."
  ),
  list(
    id = 2L,
    name = "3-state Student-t HMM",
    check = "Heavy-tail robustness",
    data_file = primary_file,
    frequency_seconds = cfg$bucket_size_seconds,
    n_states = cfg$n_states,
    response_variable = cfg$response_variable,
    flow_variable = cfg$flow_variable,
    emission_family = "student_t",
    student_t_df = cfg$student_t_df,
    trade_buckets_only = FALSE,
    conclusion = "Checks whether heavy-tailed residuals change state inference."
  ),
  list(
    id = 3L,
    name = "3-state 5-second trade-bucket HMM",
    check = "Ultra-high-frequency robustness",
    data_file = trade_5s_file,
    frequency_seconds = cfg$reconstruction_bucket_size_seconds,
    n_states = cfg$n_states,
    response_variable = cfg$response_variable,
    flow_variable = cfg$flow_variable,
    emission_family = "gaussian",
    student_t_df = NA_real_,
    trade_buckets_only = TRUE,
    conclusion = "Checks whether the regime structure survives at 5-second trade buckets."
  ),
  list(
    id = 4L,
    name = "3-state signed-volume HMM",
    check = "Order-flow variable choice",
    data_file = primary_file,
    frequency_seconds = cfg$bucket_size_seconds,
    n_states = cfg$n_states,
    response_variable = cfg$response_variable,
    flow_variable = "q_t",
    emission_family = "gaussian",
    student_t_df = NA_real_,
    trade_buckets_only = FALSE,
    conclusion = "Checks whether OFI matters relative to signed trade volume."
  ),
  list(
    id = 5L,
    name = "3-state 60-second HMM",
    check = "Sampling-frequency robustness",
    data_file = freq_60_file,
    frequency_seconds = 60L,
    n_states = cfg$n_states,
    response_variable = cfg$response_variable,
    flow_variable = cfg$flow_variable,
    emission_family = "gaussian",
    student_t_df = NA_real_,
    trade_buckets_only = FALSE,
    conclusion = "Checks whether 30 seconds is in the right aggregation range."
  )
)

make_robustness_error_row <- function(spec, msg) {
  # Keep failed specifications visible in the output table instead of silently
  # dropping them. This makes pipeline problems easier to diagnose.
  list(fit = NULL, summary = data.table(
    model = spec$name, states = spec$n_states,
    robustness_id = spec$id, robustness_check = spec$check,
    frequency_seconds = spec$frequency_seconds,
    response_variable = spec$response_variable,
    flow_variable = spec$flow_variable,
    emission_family = spec$emission_family,
    covariance_mode = "full", student_t_df = spec$student_t_df,
    logLik = NA_real_, parameters = NA_integer_,
    observations = NA_integer_, AIC = NA_real_, BIC = NA_real_,
    convergence_code = NA_integer_, convergence_message = msg,
    trade_buckets_only = spec$trade_buckets_only,
    state_shares = NA_character_, top_lambda_state = NA_integer_,
    conclusion = "Model did not complete cleanly."
  ))
}

if (requireNamespace("doParallel", quietly = TRUE) && requireNamespace("foreach", quietly = TRUE)) {
  library(foreach)
  library(doParallel)
  # The five robustness fits are independent, so parallel execution saves time
  # without changing any model output. Each worker reloads the shared functions.
  n_rob_workers <- min(length(robustness_specs), max(parallel::detectCores() - 2L, 1L))
  message("Running ", length(robustness_specs), " robustness models in parallel on ", n_rob_workers, " workers.")
  cl_rob <- parallel::makeCluster(n_rob_workers)
  parallel::clusterExport(cl_rob, varlist = "code_dir", envir = environment())
  parallel::clusterEvalQ(cl_rob, {
    library(data.table)
    source(file.path(code_dir, "01_project_settings_and_functions.R"))
  })
  parallel::clusterExport(cl_rob, varlist = c(
    "load_model_data", "fit_robustness_spec", "make_robustness_error_row",
    "robustness_em_maxit", "robustness_em_tolerance", "robustness_n_starts", "cfg"
  ), envir = environment())
  doParallel::registerDoParallel(cl_rob)
  robustness_list <- foreach::foreach(spec = robustness_specs, .errorhandling = "pass") %dopar% {
    tryCatch(fit_robustness_spec(spec),
      error = function(e) make_robustness_error_row(spec, e$message)
    )
  }
  parallel::stopCluster(cl_rob)
} else {
  warning("doParallel unavailable; running robustness models sequentially.")
  robustness_list <- lapply(robustness_specs, function(spec) {
    tryCatch(fit_robustness_spec(spec), error = function(e) {
      warning("Robustness model failed: ", spec$name, " -- ", e$message)
      make_robustness_error_row(spec, e$message)
    })
  })
}

robustness_comparison <- rbindlist(lapply(robustness_list, `[[`, "summary"), use.names = TRUE, fill = TRUE)
setorder(robustness_comparison, BIC, na.last = TRUE)
# Delta BIC is computed within this robustness table only. The paper warns that
# values from different sampling frequencies are descriptive rather than a strict
# likelihood comparison because the underlying observation counts differ.
robustness_comparison[is.finite(BIC), delta_BIC := BIC - min(BIC, na.rm = TRUE)]
write_table(robustness_comparison, file.path(robustness_root, "hmm_robustness_model_comparison.csv"))

plot_dt <- robustness_comparison[is.finite(delta_BIC)]
if (nrow(plot_dt) > 0L) {
  # Delta BIC differences are large, so the plot uses log10(1 + Delta BIC) for
  # legibility while the text labels keep the original BIC scale.
  plot_dt[, delta_BIC_plot := log10(delta_BIC + 1)]
  plot_dt[, model_label := factor(model, levels = rev(model))]
  plot_dt[, group := fifelse(
    emission_family == "student_t", "Student-t emissions",
    fifelse(
      flow_variable == "q_t", "Signed volume",
      fifelse(
        frequency_seconds == 5L, "5-second trade buckets",
        fifelse(frequency_seconds == 60L, "60-second frequency", "2-state model")
      )
    )
  )]
  robustness_plot <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = delta_BIC_plot, y = model_label, fill = group)) +
    ggplot2::geom_col(width = 0.62, colour = paper_colours$dark, linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", delta_BIC)), hjust = -0.10, size = 3.2, colour = paper_colours$dark) +
    ggplot2::scale_fill_manual(values = c(
      "2-state model" = paper_colours$calm,
      "Student-t emissions" = paper_colours$stressed,
      "5-second trade buckets" = paper_colours$purple,
      "Signed volume" = paper_colours$orange,
      "60-second frequency" = paper_colours$green
    )) +
    ggplot2::coord_cartesian(xlim = c(0, max(plot_dt$delta_BIC_plot, na.rm = TRUE) * 1.18 + 0.05)) +
    ggplot2::labs(x = expression(log[10] * (1 + Delta * "BIC")), y = NULL) +
    paper_theme() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "bottom")
  save_paper_plot(robustness_plot, file.path(robustness_figures_root, "robustness_delta_bic_five_models.png"), width = 9.5, height = 5.6)
}

message("Done.")
message("Robustness folder: ", robustness_root)
