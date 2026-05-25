###############################################################################
# Diagnostics and assumptions for the MOEX model-ready panel
#
# What this script does:
#   1. Loads the primary dataset defined in the shared settings.
#   2. Checks descriptive statistics, intraday profiles, stationarity proxies,
#      conditional-Gaussian diagnostics, regressors, price impact, and outliers.
#   3. Saves compact tables and paper-style figures. It does not generate
#      narrative text and does not estimate the HMM.
###############################################################################

library(data.table)

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

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for diagnostic figures.")
}

data_file <- combined_data_file(cfg, frequency = cfg$bucket_size_seconds)
if (!file.exists(data_file)) stop("Primary combined dataset not found. Run 03_moex_combine_processed_data.R first: ", data_file)
dt <- fread(data_file)
setorder(dt, seccode, trade_date, seconds_from_midnight)

# Result folders are deliberately split by diagnostic theme. This makes it easy
# to audit what supports each assumption in the paper.
dirs <- list(
  descriptive = cfg$results_descriptive_root,
  seasonality_tables = cfg$results_seasonality_root,
  seasonality_figures = cfg$results_seasonality_root,
  stationarity_tables = cfg$results_stationarity_root,
  stationarity_figures = cfg$results_stationarity_root,
  normality_uni = cfg$results_normality_root,
  normality_bi = cfg$results_normality_root,
  regressors_tables = cfg$results_regressors_root,
  regressors_figures = cfg$results_regressors_root,
  price_tables = cfg$results_price_impact_root,
  price_figures = cfg$results_price_impact_root,
  outliers = cfg$results_outliers_root
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

acf_lag_max <- as.integer(Sys.getenv("MOEX_ACF_LAG_MAX", unset = "100"))
max_points_in_scatter <- as.integer(Sys.getenv("MOEX_MAX_SCATTER_POINTS", unset = "50000"))

# =========================
# 1. SMALL DIAGNOSTIC HELPERS
# =========================

sample_for_plot <- function(x, n = max_points_in_scatter) {
  # Scatter/QQ plots do not need almost one million points to show the pattern.
  # Sampling keeps file sizes small and plots readable.
  if (nrow(x) <= n) {
    return(x)
  }
  set.seed(cfg$random_seed)
  x[sample.int(.N, n)]
}

skewness_manual <- function(x) {
  x <- clean_numeric(x)
  if (length(x) < 3L) {
    return(NA_real_)
  }
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  mean((x - mean(x))^3) / s^3
}

kurtosis_manual <- function(x) {
  x <- clean_numeric(x)
  if (length(x) < 4L) {
    return(NA_real_)
  }
  s <- sd(x)
  if (!is.finite(s) || s <= 0) {
    return(NA_real_)
  }
  mean((x - mean(x))^4) / s^4
}

summary_stats <- function(data, variables) {
  # One common summary function keeps table definitions identical across the
  # descriptive, normality, regressor, and outlier diagnostics.
  rbindlist(lapply(variables, function(v) {
    x <- clean_numeric(data[[v]])
    qs <- quantile(x, probs = c(0.001, 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99, 0.999), na.rm = TRUE, names = FALSE)
    data.table(
      variable = v,
      N = length(x),
      missing = sum(is.na(data[[v]])),
      mean = mean(x),
      sd = sd(x),
      min = min(x),
      p0_1 = qs[1L],
      p1 = qs[2L],
      p5 = qs[3L],
      p25 = qs[4L],
      median = qs[5L],
      p75 = qs[6L],
      p95 = qs[7L],
      p99 = qs[8L],
      p99_9 = qs[9L],
      max = max(x),
      skewness = skewness_manual(x),
      kurtosis = kurtosis_manual(x),
      excess_kurtosis = kurtosis_manual(x) - 3
    )
  }), use.names = TRUE, fill = TRUE)
}

mean_acf_by_sequence <- function(data, variable, lag_max) {
  # Autocorrelation is computed within each stock-day, then averaged. This avoids
  # treating overnight gaps or stock boundaries as ordinary adjacent buckets.
  acfs <- data[,
    {
      x <- clean_numeric(get(variable))
      if (length(x) <= lag_max + 1L || sd(x) == 0) {
        list(acf_values = list(rep(NA_real_, lag_max + 1L)))
      } else {
        list(acf_values = list(as.numeric(acf(x, lag.max = lag_max, plot = FALSE, na.action = na.pass)$acf)))
      }
    },
    by = sequence_id
  ]
  mat <- do.call(rbind, acfs$acf_values)
  data.table(
    variable = variable,
    lag = 0:lag_max,
    mean_acf = colMeans(mat, na.rm = TRUE),
    min_acf = apply(mat, 2L, min, na.rm = TRUE),
    max_acf = apply(mat, 2L, max, na.rm = TRUE)
  )
}

plot_line <- function(profile_dt, y, file_name, ylab) {
  plot_dt <- profile_dt[, .(hour, value = get(y))]
  p <- ggplot2::ggplot(plot_dt[is.finite(value)], ggplot2::aes(x = hour, y = value)) +
    ggplot2::geom_line(colour = paper_colours$calm, linewidth = 0.8) +
    ggplot2::labs(x = "Hour of day", y = ylab) +
    paper_theme()
  save_paper_plot(p, file.path(dirs$seasonality_figures, file_name), width = 8.8, height = 5.2)
}

plot_acf <- function(acf_dt, variable) {
  p <- ggplot2::ggplot(acf_dt[lag > 0], ggplot2::aes(x = lag, y = mean_acf)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, y = 0, yend = mean_acf), colour = paper_colours$calm, linewidth = 0.55) +
    ggplot2::labs(x = "Lag", y = "Mean within-sequence autocorrelation") +
    paper_theme()
  save_paper_plot(p, file.path(dirs$stationarity_figures, sprintf("acf_%s.png", variable)), width = 8.8, height = 5.2)
}

normal_density_hist <- function(x, path, xlab) {
  # The histogram is trimmed only for display. The CSV normality statistics are
  # computed on the full variable.
  x <- clean_numeric(x)
  q <- quantile(x, c(0.001, 0.999), na.rm = TRUE, names = FALSE)
  x_plot <- x[x >= q[1L] & x <= q[2L]]
  normal_dt <- data.table(value = seq(min(x_plot), max(x_plot), length.out = 400L))
  normal_dt[, density := dnorm(value, mean = mean(x_plot), sd = sd(x_plot))]
  p <- ggplot2::ggplot(data.table(value = x_plot), ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 80, fill = "#D9DEE7", colour = "white", linewidth = 0.2) +
    ggplot2::geom_line(data = normal_dt, ggplot2::aes(y = density), colour = paper_colours$stressed, linewidth = 0.9) +
    ggplot2::labs(x = xlab, y = "Density") +
    paper_theme()
  save_paper_plot(p, path, width = 8.8, height = 5.2)
}

qq_normal_plot <- function(x, path, ylab) {
  # QQ plots are a visual check of whether a Gaussian emission is a reasonable
  # approximation after aggregation and transformations.
  x <- clean_numeric(x)
  if (length(x) > max_points_in_scatter) {
    set.seed(cfg$random_seed)
    x <- sample(x, max_points_in_scatter)
  }
  qq <- qqnorm(x, plot.it = FALSE)
  qq_dt <- data.table(theoretical = qq$x, sample = qq$y)
  q_sample <- quantile(qq_dt$sample, c(0.25, 0.75), na.rm = TRUE)
  slope <- diff(q_sample) / diff(qnorm(c(0.25, 0.75)))
  intercept <- q_sample[1L] - slope * qnorm(0.25)
  p <- ggplot2::ggplot(qq_dt, ggplot2::aes(x = theoretical, y = sample)) +
    ggplot2::geom_point(alpha = 0.30, size = 0.65, colour = paper_colours$calm) +
    ggplot2::geom_abline(intercept = intercept, slope = slope, colour = paper_colours$stressed, linewidth = 0.8) +
    ggplot2::labs(x = "Theoretical normal quantiles", y = ylab) +
    paper_theme()
  save_paper_plot(p, path, width = 8.8, height = 5.2)
}

# =========================
# 2. DESCRIPTIVE STATISTICS
# =========================

descriptive_variables <- intersect(c(
  "r_t", "r_t_winsorized", "spread_tilde", "spread", "relative_spread",
  "log_relative_spread", "q_t", "event_level_ofi", "ofi_depth_scaled",
  "abs_ofi_tilde", "depth", "depth_tilde", "RV_t", "rv_tilde",
  "number_of_trades", "trade_volume"
), names(dt))
descriptive_statistics <- summary_stats(dt, descriptive_variables)
write_table(descriptive_statistics, file.path(dirs$descriptive, "descriptive_statistics.csv"))

# =========================
# 3. INTRADAY SEASONALITY CHECKS
# =========================

# These profiles are computed before state estimation. They show which patterns
# are removed by de-seasonalization and which regular intraday features remain
# visible in raw variables.
seasonality_variables <- intersect(c(
  "relative_spread", "spread_tilde", "r_t_winsorized", "event_level_ofi",
  "ofi_depth_scaled", "abs_ofi_tilde", "depth", "depth_tilde",
  "RV_t", "rv_tilde", "number_of_trades", "trade_volume"
), names(dt))
intraday_profiles <- dt[, lapply(.SD, mean, na.rm = TRUE), by = tod_bin, .SDcols = seasonality_variables]
intraday_profiles[, hour := tod_bin / 3600]
setcolorder(intraday_profiles, c("tod_bin", "hour", seasonality_variables))
write_table(intraday_profiles, file.path(dirs$seasonality_tables, "intraday_profiles.csv"))

plot_line(intraday_profiles, "relative_spread", "intraday_relative_spread.png", "Relative spread")
plot_line(intraday_profiles, "spread_tilde", "intraday_spread_tilde.png", "Transformed spread")
plot_line(intraday_profiles, "r_t_winsorized", "intraday_winsorized_return.png", "Winsorized return")
plot_line(intraday_profiles, "ofi_depth_scaled", "intraday_ofi_depth_scaled.png", "Depth-scaled OFI")
plot_line(intraday_profiles, "abs_ofi_tilde", "intraday_abs_ofi_tilde.png", "Standardized |OFI|")
plot_line(intraday_profiles, "rv_tilde", "intraday_rv_tilde.png", "Standardized realized volatility")

trade_activity <- intraday_profiles[, .(
  hour,
  `Number of trades` = as.numeric(scale(number_of_trades)),
  `Trade volume` = as.numeric(scale(trade_volume))
)]
trade_activity_long <- melt(trade_activity, id.vars = "hour", variable.name = "series", value.name = "value")
trade_plot <- ggplot2::ggplot(trade_activity_long[is.finite(value)], ggplot2::aes(x = hour, y = value, colour = series)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::scale_colour_manual(values = c("Number of trades" = paper_colours$calm, "Trade volume" = paper_colours$stressed)) +
  ggplot2::labs(x = "Hour of day", y = "Standardized profile") +
  paper_theme()
save_paper_plot(trade_plot, file.path(dirs$seasonality_figures, "intraday_trade_activity.png"), width = 8.8, height = 5.2)

# =========================
# 4. STATIONARITY AND PERSISTENCE DIAGNOSTICS
# =========================

stationarity_variables <- intersect(c(
  "r_t_winsorized", "spread_tilde", "ofi_depth_scaled",
  "abs_ofi_tilde", "depth_tilde", "rv_tilde"
), names(dt))
acf_summary <- rbindlist(lapply(stationarity_variables, function(v) {
  # The first 20 lags are saved to the CSV used in the paper; the figures are
  # still drawn up to acf_lag_max for a fuller diagnostic view.
  out <- mean_acf_by_sequence(dt, v, acf_lag_max)
  plot_acf(out, v)
  out[lag <= 20]
}), use.names = TRUE, fill = TRUE)
write_table(acf_summary, file.path(dirs$stationarity_tables, "acf_summary.csv"))

stationarity_tests <- rbindlist(lapply(stationarity_variables, function(v) {
  x <- clean_numeric(dt[[v]])
  # Unit-root tests are run on a deterministic subsample for speed. With this
  # sample size, the tests are directional diagnostics rather than the main
  # identification argument.
  x <- if (length(x) > 50000L) x[seq(1L, length(x), length.out = 50000L)] else x
  adf <- if (requireNamespace("tseries", quietly = TRUE)) {
    tryCatch(tseries::adf.test(x, alternative = "stationary"), error = function(e) NULL)
  } else {
    NULL
  }
  kpss <- if (requireNamespace("tseries", quietly = TRUE)) {
    tryCatch(tseries::kpss.test(x, null = "Level"), error = function(e) NULL)
  } else {
    NULL
  }
  data.table(
    variable = v,
    test_sample_n = length(x),
    adf_statistic = if (!is.null(adf)) unname(adf$statistic) else NA_real_,
    adf_p_value = if (!is.null(adf)) adf$p.value else NA_real_,
    kpss_statistic = if (!is.null(kpss)) unname(kpss$statistic) else NA_real_,
    kpss_p_value = if (!is.null(kpss)) kpss$p.value else NA_real_
  )
}), use.names = TRUE, fill = TRUE)
write_table(stationarity_tests, file.path(dirs$stationarity_tables, "stationarity_tests.csv"))

# =========================
# 5. NORMALITY DIAGNOSTICS
# =========================

normality_univariate <- summary_stats(dt, c("r_t_winsorized", "spread_tilde"))[, .(
  variable, N, mean, sd, skewness, kurtosis, excess_kurtosis
)]
if (requireNamespace("tseries", quietly = TRUE)) {
  # Jarque-Bera is expected to reject in high-frequency data. It is saved because
  # it documents the gap between unconditional normality and the weaker
  # conditional-Gaussian working assumption used by the HMM.
  normality_univariate[, jarque_bera_p_value := vapply(variable, function(v) {
    out <- tryCatch(tseries::jarque.bera.test(clean_numeric(dt[[v]])), error = function(e) NULL)
    if (is.null(out)) NA_real_ else out$p.value
  }, numeric(1L))]
}
write_table(normality_univariate, file.path(dirs$normality_uni, "normality_univariate.csv"))

normal_density_hist(dt$r_t_winsorized, file.path(dirs$normality_uni, "hist_r_t_winsorized_normal_overlay.png"), "Winsorized return")
normal_density_hist(dt$spread_tilde, file.path(dirs$normality_uni, "hist_spread_tilde_normal_overlay.png"), "Transformed spread")
qq_normal_plot(dt$r_t_winsorized, file.path(dirs$normality_uni, "qq_r_t_winsorized.png"), "Winsorized return quantiles")
qq_normal_plot(dt$spread_tilde, file.path(dirs$normality_uni, "qq_spread_tilde.png"), "Transformed spread quantiles")

emission_cov <- cov(dt[, .(spread_tilde, r_t_winsorized)], use = "complete.obs")
emission_cor <- cor(dt[, .(spread_tilde, r_t_winsorized)], use = "complete.obs")
write_table(data.table(metric = c("covariance", "correlation"), spread_return = c(emission_cov[1, 2], emission_cor[1, 2])), file.path(dirs$normality_bi, "emission_covariance_correlation.csv"))

# The scatter plot is not used as a regression result. It is a quick visual
# diagnostic for the bivariate emission variables before state conditioning.
scatter_dt <- sample_for_plot(dt[is.finite(spread_tilde) & is.finite(r_t_winsorized), .(spread_tilde, r_t_winsorized)])
scatter_plot <- ggplot2::ggplot(scatter_dt, ggplot2::aes(x = spread_tilde, y = r_t_winsorized)) +
  ggplot2::geom_point(alpha = 0.22, size = 0.55, colour = paper_colours$calm) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = paper_colours$stressed, linewidth = 0.8) +
  ggplot2::labs(x = "Transformed spread", y = "Winsorized return") +
  paper_theme()
save_paper_plot(scatter_plot, file.path(dirs$normality_bi, "scatter_emission_variables.png"), width = 8.8, height = 5.2)

# =========================
# 6. REGRESSOR DIAGNOSTICS
# =========================

regressor_variables <- intersect(c("depth_tilde", "rv_tilde", "abs_ofi_tilde", "open_dummy", "close_dummy"), names(dt))
if (!"open_dummy" %in% names(dt)) dt[, open_dummy := as.integer(seconds_from_midnight < cfg$open_end_seconds)]
if (!"close_dummy" %in% names(dt)) dt[, close_dummy := as.integer(seconds_from_midnight >= cfg$close_start_seconds)]
regressor_variables <- intersect(c("depth_tilde", "rv_tilde", "abs_ofi_tilde", "open_dummy", "close_dummy"), names(dt))
write_table(summary_stats(dt, regressor_variables), file.path(dirs$regressors_tables, "transition_regressor_summary.csv"))

correlation_variables <- intersect(c("depth_tilde", "rv_tilde", "abs_ofi_tilde", "number_of_trades", "trade_volume"), names(dt))
cor_mat <- cor(dt[, ..correlation_variables], use = "complete.obs")
write_table(as.data.table(cor_mat, keep.rownames = "variable"), file.path(dirs$regressors_tables, "regressor_correlation_matrix.csv"))

# The heatmap checks whether the transition regressors are visually dominated by
# one common factor. It supports keeping depth, volatility, OFI pressure, and
# activity as distinct diagnostics.
heat_dt <- as.data.table(as.table(cor_mat))
setnames(heat_dt, c("variable_1", "variable_2", "correlation"))
label_map <- c(
  depth_tilde = "Depth\ntransf.",
  rv_tilde = "RV\ntransf.",
  abs_ofi_tilde = "|OFI|\nstd.",
  number_of_trades = "Number\nof trades",
  trade_volume = "Trade\nvolume"
)
heat_dt[, variable_1_label := label_map[as.character(variable_1)]]
heat_dt[, variable_2_label := label_map[as.character(variable_2)]]
heat_dt[is.na(variable_1_label), variable_1_label := as.character(variable_1)]
heat_dt[is.na(variable_2_label), variable_2_label := as.character(variable_2)]
heat_plot <- ggplot2::ggplot(heat_dt, ggplot2::aes(x = variable_2_label, y = variable_1_label, fill = correlation)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.9) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", correlation), color = abs(correlation) > 0.55), size = 4.0, fontface = "bold") +
  ggplot2::scale_fill_gradient2(low = paper_colours$calm, mid = "white", high = paper_colours$stressed, midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  ggplot2::scale_color_manual(values = c("FALSE" = paper_colours$dark, "TRUE" = "white"), guide = "none") +
  ggplot2::coord_fixed() +
  ggplot2::labs(x = NULL, y = NULL) +
  paper_theme(base_size = 13) +
  ggplot2::theme(legend.position = "right", panel.grid = ggplot2::element_blank())
save_paper_plot(heat_plot, file.path(dirs$regressors_figures, "regressor_correlation_heatmap.png"), width = 8.4, height = 6.4)

max_corr <- apply(abs(cor_mat - diag(ncol(cor_mat))), 1L, max, na.rm = TRUE)
write_table(data.table(variable = names(max_corr), max_absolute_pairwise_correlation = as.numeric(max_corr)), file.path(dirs$regressors_tables, "multicollinearity_pairwise.csv"))

# =========================
# 7. PRICE-IMPACT SANITY CHECKS
# =========================

coef_rows_from_lm <- function(fit, se_type, vcov_matrix = NULL) {
  if (!is.null(vcov_matrix) && requireNamespace("lmtest", quietly = TRUE)) {
    ct <- lmtest::coeftest(fit, vcov. = vcov_matrix)
  } else {
    ct <- summary(fit)$coefficients
  }
  data.table(
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    std_error = ct[, "Std. Error"],
    t_statistic = ct[, "t value"],
    p_value = ct[, ncol(ct)],
    se_type = se_type,
    r_squared = summary(fit)$r.squared,
    n = length(fit$residuals)
  )
}

fit_price_impact <- function(flow_col, flow_label) {
  # This is a simple pre-HMM sanity check: signed order-flow pressure should have
  # the expected positive relation with short-horizon returns.
  x <- dt[is.finite(r_t_winsorized) & is.finite(get(flow_col)), .(
    trade_date,
    response = r_t_winsorized, flow = get(flow_col)
  )]
  fit <- lm(response ~ flow, data = x)
  rows <- list(coef_rows_from_lm(fit, "OLS"))
  if (requireNamespace("sandwich", quietly = TRUE) && requireNamespace("lmtest", quietly = TRUE)) {
    rows[[length(rows) + 1L]] <- coef_rows_from_lm(fit, "HC1 robust", sandwich::vcovHC(fit, type = "HC1"))
    if (uniqueN(x$trade_date) >= 2L) {
      rows[[length(rows) + 1L]] <- coef_rows_from_lm(fit, "cluster by trade_date (HC1)", sandwich::vcovCL(fit, cluster = x$trade_date, type = "HC1"))
    }
  }
  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, `:=`(flow_measure = flow_label, flow_column = flow_col)]
  out
}

price_impact_rows <- rbindlist(list(
  fit_price_impact("ofi_depth_scaled", "depth-scaled OFI"),
  fit_price_impact("q_t", "signed trade volume")
), use.names = TRUE, fill = TRUE)
write_table(price_impact_rows, file.path(dirs$price_tables, "price_impact_robustness.csv"))

price_dt <- sample_for_plot(dt[is.finite(r_t_winsorized) & is.finite(ofi_depth_scaled), .(r_t_winsorized, ofi_depth_scaled)])
price_plot <- ggplot2::ggplot(price_dt, ggplot2::aes(x = ofi_depth_scaled, y = r_t_winsorized)) +
  ggplot2::geom_point(alpha = 0.18, size = 0.55, colour = paper_colours$calm) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = paper_colours$stressed, linewidth = 0.85) +
  ggplot2::labs(x = "Depth-scaled OFI", y = "Winsorized return") +
  paper_theme()
save_paper_plot(price_plot, file.path(dirs$price_figures, "ofi_depth_scaled_vs_return.png"), width = 8.8, height = 5.2)

price_fit <- lm(r_t_winsorized ~ ofi_depth_scaled, data = dt[is.finite(r_t_winsorized) & is.finite(ofi_depth_scaled)])
resid_dt <- data.table(residual = residuals(price_fit))
resid_hist <- ggplot2::ggplot(resid_dt, ggplot2::aes(x = residual)) +
  ggplot2::geom_histogram(bins = 100, fill = "#D9DEE7", colour = "white", linewidth = 0.2) +
  ggplot2::labs(x = "Residual", y = "Count") +
  paper_theme()
save_paper_plot(resid_hist, file.path(dirs$price_figures, "price_impact_residual_hist.png"), width = 8.8, height = 5.2)

dt[, price_impact_residual := NA_real_]
dt[is.finite(r_t_winsorized) & is.finite(ofi_depth_scaled), price_impact_residual := residuals(price_fit)]
# Residual ACF checks whether the simple linear price-impact regression leaves
# strong short-run dependence. The HMM later handles richer state dependence.
resid_acf <- mean_acf_by_sequence(dt[is.finite(price_impact_residual)], "price_impact_residual", acf_lag_max)
write_table(resid_acf[lag <= 20], file.path(dirs$price_tables, "price_impact_residual_acf_summary.csv"))
plot_acf(resid_acf, "price_impact_residual")
invisible(file.rename(
  file.path(dirs$stationarity_figures, "acf_price_impact_residual.png"),
  file.path(dirs$price_figures, "price_impact_residual_acf.png")
))

# =========================
# 8. OUTLIER COUNTS
# =========================

outlier_variables <- intersect(c("r_t_winsorized", "spread_tilde", "ofi_depth_scaled", "abs_ofi_tilde", "rv_tilde"), names(dt))
outlier_diagnostics <- rbindlist(lapply(outlier_variables, function(v) {
  # The paper does not delete observations by these percentiles. They are only
  # reported to show which variables are heavy-tailed before HMM estimation.
  x <- clean_numeric(dt[[v]])
  q001 <- quantile(x, 0.001, names = FALSE)
  q999 <- quantile(x, 0.999, names = FALSE)
  q01 <- quantile(x, 0.01, names = FALSE)
  q99 <- quantile(x, 0.99, names = FALSE)
  data.table(
    variable = v,
    N = length(x),
    p0_1 = q001,
    p1 = q01,
    p99 = q99,
    p99_9 = q999,
    outside_p0_1_p99_9 = sum(x < q001 | x > q999),
    outside_p1_p99 = sum(x < q01 | x > q99)
  )
}), use.names = TRUE, fill = TRUE)
write_table(outlier_diagnostics, file.path(dirs$outliers, "outlier_diagnostics.csv"))

message("Done.")
message("Diagnostics folders: ", cfg$results_root)
