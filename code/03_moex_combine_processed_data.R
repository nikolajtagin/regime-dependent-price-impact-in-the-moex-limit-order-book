###############################################################################
# MOEX processed-data combination
#
# What this script does:
#   1. Reads 5-second stock-day files produced by 02_moex_orderlog_processing.R.
#   2. Aggregates them to the model frequencies defined in the shared settings.
#   3. Constructs final model variables: winsorized returns, depth-scaled OFI,
#      intraday-standardized covariates, and stock-day sequence identifiers.
#   4. Saves the primary 30-second dataset, the 5-second trade-bucket robustness
#      dataset, and the 60-second robustness dataset.
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

dir.create(cfg$combined_data_root, recursive = TRUE, showWarnings = FALSE)

# This script starts from daily 5-second files. It does not read the raw MOEX
# logs again. If a daily file is missing, the fix is to rerun script 02.
tickers <- toupper(split_env_vector(Sys.getenv("MOEX_TICKERS", unset = paste(cfg$tickers, collapse = ","))))
cfg$tickers <- tickers
base_frequency <- cfg$reconstruction_bucket_size_seconds
target_frequencies <- unique(c(cfg$bucket_size_seconds, cfg$robustness_frequencies_seconds))

message("Combining processed files for: ", paste(tickers, collapse = ", "))
message("Base frequency: ", base_frequency, " seconds")
message("Target frequencies: ", paste(target_frequencies, collapse = ", "), " seconds")

# =========================
# 1. READ BASE 5-SECOND FILES
# =========================

processed_pattern <- sprintf("^processed_%ss_data_([a-z0-9]+)_([0-9]{2}\\.[0-9]{2}\\.[0-9]{4})\\.csv$", base_frequency)
processed_files <- list.files(cfg$daily_processed_root, pattern = processed_pattern, recursive = TRUE, full.names = TRUE)
if (length(processed_files) == 0L) {
  stop("No processed files found. Run 02_moex_orderlog_processing.R first.")
}

extract_date_code <- function(path) {
  # Daily files use DD.MM.YYYY in the file name because that is the readable
  # style used in the processed data folder. Convert it back to YYYYMMDD for
  # sorting and filtering.
  dotted <- sub("^processed_[0-9]+s_data_[a-z0-9]+_([0-9]{2}\\.[0-9]{2}\\.[0-9]{4})\\.csv$", "\\1", basename(path))
  format(as.Date(dotted, format = "%d.%m.%Y"), "%Y%m%d")
}
extract_ticker_code <- function(path) {
  toupper(sub("^processed_[0-9]+s_data_([a-z0-9]+)_[0-9]{2}\\.[0-9]{2}\\.[0-9]{4}\\.csv$", "\\1", basename(path)))
}

date_codes <- vapply(processed_files, extract_date_code, character(1L))
ticker_codes <- vapply(processed_files, extract_ticker_code, character(1L))
keep <- date_codes >= cfg$sample_start_code & date_codes <= cfg$sample_end_code & ticker_codes %in% tickers
processed_files <- processed_files[keep]
date_codes <- date_codes[keep]
ticker_codes <- ticker_codes[keep]
if (length(processed_files) == 0L) stop("No processed files remain after sample-date and ticker filtering.")

file_order <- order(ticker_codes, date_codes)
processed_files <- processed_files[file_order]
date_codes <- date_codes[file_order]
ticker_codes <- ticker_codes[file_order]

base_dt <- rbindlist(lapply(seq_along(processed_files), function(i) {
  # Add stock/date metadata from the file name. This is safer than relying on
  # the daily files alone, and it makes the combined panel self-contained.
  x <- fread(processed_files[i])
  x[, `:=`(
    seccode = ticker_codes[i],
    trade_date = date_codes[i],
    calendar_date = as.IDate(date_codes[i], format = "%Y%m%d")
  )]
  x
}), use.names = TRUE, fill = TRUE)

setorder(base_dt, seccode, trade_date, seconds_from_midnight)

# =========================
# 2. AGGREGATION AND TRANSFORMATIONS
# =========================

last_non_na <- function(x) {
  # Quote variables are end-of-bucket states, so the correct aggregation rule is
  # the last observed non-missing value inside the larger bucket.
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else x[length(x)]
}

sum_or_zero <- function(x) {
  # Event counts and volume-like variables accumulate over the bucket. Empty
  # buckets should contribute zero, not NA.
  if (all(is.na(x))) 0 else sum(x, na.rm = TRUE)
}

aggregate_to_frequency <- function(x, frequency) {
  # Collapse the 5-second reconstruction to the target model frequency. Quotes
  # use last observation carried within the bucket; flows and counts are summed.
  x <- copy(x)
  x[, bucket_start := floor(seconds_from_midnight / frequency) * frequency]
  x[, bucket_end := bucket_start + frequency]
  x[, seconds_from_midnight := bucket_start]
  key_cols <- c("seccode", "trade_date", "calendar_date", "bucket_start", "bucket_end", "seconds_from_midnight")

  quote_cols <- intersect(c(
    "best_bid", "best_ask", "best_bid_depth", "best_ask_depth",
    "top5_bid_depth", "top5_ask_depth", "top5_depth",
    "queue_imbalance", "top5_queue_imbalance"
  ), names(x))
  count_cols <- intersect(c(
    "number_of_events", "number_of_order_adds", "number_of_cancellations",
    "number_of_trades", "trade_volume",
    "q_t", "q_t_tick", "q_t_buysell", "q_t_buysell_opposite",
    "signed_trade_volume_proxy", "signed_trade_volume_tick",
    "signed_trade_volume_buysell", "signed_trade_volume_buysell_opposite",
    "abs_signed_trade_volume_proxy", "event_level_ofi"
  ), names(x))

  quote_part <- if (length(quote_cols) > 0L) {
    x[, lapply(.SD, last_non_na), by = key_cols, .SDcols = quote_cols]
  } else {
    unique(x[, ..key_cols])
  }

  count_part <- if (length(count_cols) > 0L) {
    x[, lapply(.SD, sum_or_zero), by = key_cols, .SDcols = count_cols]
  } else {
    unique(x[, ..key_cols])
  }

  out <- merge(quote_part, count_part, by = key_cols, all = TRUE, sort = FALSE)

  setorder(out, seccode, trade_date, seconds_from_midnight)
  out[, sequence_id := paste(seccode, trade_date, sep = "_")]
  out[, row_in_sequence := seq_len(.N), by = sequence_id]
  # All model variables are recomputed after aggregation. A 30-second return is
  # therefore a true 30-second midprice return, not a sum of 5-second returns.
  out[, depth := best_bid_depth + best_ask_depth]
  out[, spread := best_ask - best_bid]
  out[, midprice := (best_bid + best_ask) / 2]
  out[, relative_spread := spread / midprice]
  out[, log_relative_spread := fifelse(relative_spread > 0, log(relative_spread), NA_real_)]
  out[, r_t := 10000 * (log(midprice) - shift(log(midprice))), by = sequence_id]
  out[, abs_q_t := abs(q_t)]
  out[, abs_event_level_ofi := abs(event_level_ofi)]
  out[, q_depth_scaled := fifelse(depth > 0, q_t / depth, NA_real_)]
  out[, ofi_depth_scaled := fifelse(depth > 0, event_level_ofi / depth, NA_real_)]
  out[, abs_q_depth_scaled := abs(q_depth_scaled)]
  out[, abs_ofi_depth_scaled := abs(ofi_depth_scaled)]
  out[, crossed_quote_flag := best_ask < best_bid]
  out[, locked_quote_flag := best_ask == best_bid]
  out[, valid_topbook_flag := !is.na(best_bid) & !is.na(best_ask) & best_bid > 0 & best_ask > best_bid]

  out[, RV_t := frollsum(
    shift(r_t^2),
    n = cfg$rv_window_buckets,
    align = "right",
    fill = NA_real_,
    na.rm = TRUE,
    partial = TRUE
  ), by = sequence_id]
  out[, effective_rv_window := frollsum(
    as.integer(!is.na(shift(r_t^2))),
    n = cfg$rv_window_buckets,
    align = "right",
    fill = NA_integer_,
    na.rm = TRUE,
    partial = TRUE
  ), by = sequence_id]
  out[, complete_rv_window_flag := effective_rv_window >= cfg$rv_window_buckets]
  out
}

finalize_model_variables <- function(x, frequency, trade_buckets_only = FALSE) {
  # This function applies the final sample restrictions and transformations used
  # by the HMM. It is called for the primary dataset and all robustness datasets.
  x <- copy(x)
  if (isTRUE(cfg$drop_incomplete_rv)) {
    # The first few buckets of each stock-day do not have a full six-minute
    # realized-volatility history. Dropping them keeps RV_t comparable.
    x <- x[complete_rv_window_flag == TRUE]
  }
  if (trade_buckets_only) {
    # The trade-bucket robustness check removes no-trade intervals. The primary
    # model keeps them because quote persistence is part of liquidity conditions.
    x <- x[number_of_trades > 0]
  }

  x <- x[
    valid_topbook_flag == TRUE &
      crossed_quote_flag == FALSE &
      locked_quote_flag == FALSE &
      is.finite(r_t) &
      is.finite(spread) &
      is.finite(depth) &
      spread > 0 &
      depth > 0
  ]
  if (nrow(x) == 0L) stop("No model rows remain for frequency ", frequency)

  x[, tod_bin := floor(seconds_from_midnight / cfg$tod_bin_seconds) * cfg$tod_bin_seconds]
  winsorize_by_group(x, "r_t", "r_t_winsorized", "seccode", cfg$winsor_lower, cfg$winsor_upper)

  # The paper's transformed variables remove the average intraday pattern by
  # stock and five-minute bin. This prevents the HMM from simply rediscovering
  # the usual open/close seasonality.
  x[, log_depth_for_transform := log(depth)]
  x[, log_RV_for_transform := log(RV_t + cfg$epsilon)]
  seasonality_group <- c("seccode", "tod_bin")
  center_by_group(x, "log_relative_spread", "spread_tilde", seasonality_group)
  center_by_group(x, "log_depth_for_transform", "depth_tilde", seasonality_group)
  standardize_by_group(x, "log_RV_for_transform", "rv_tilde", seasonality_group)
  standardize_by_group(x, "abs_q_t", "abs_q_tilde", seasonality_group)
  standardize_by_group(x, "event_level_ofi", "ofi_tilde", seasonality_group)
  standardize_by_group(x, "abs_event_level_ofi", "abs_ofi_tilde", seasonality_group)
  standardize_by_group(x, "number_of_trades", "trade_count_tilde", seasonality_group)

  x[, frequency_seconds := frequency]
  x[, trade_buckets_only := trade_buckets_only]
  x[, `:=`(log_depth_for_transform = NULL, log_RV_for_transform = NULL)]
  setorder(x, seccode, trade_date, seconds_from_midnight)
  x[, sequence_id := paste(seccode, trade_date, sep = "_")]
  x[, row_in_sequence := seq_len(.N), by = sequence_id]

  required <- c(
    "sequence_id", "seccode", "trade_date", "seconds_from_midnight",
    "spread_tilde", "r_t_winsorized", "ofi_depth_scaled",
    "depth_tilde", "rv_tilde", "abs_ofi_tilde"
  )
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) stop("Missing required columns after transformation: ", paste(missing, collapse = ", "))
  if (sum(!complete.cases(x[, ..required])) > 0L) stop("Missing values remain in key model variables for frequency ", frequency)
  x
}

sample_rows <- list()
for (frequency in target_frequencies) {
  message("Creating ", frequency, "-second dataset.")
  aggregated <- aggregate_to_frequency(base_dt, frequency)
  model_dt <- finalize_model_variables(aggregated, frequency, trade_buckets_only = FALSE)
  out_path <- combined_data_file(cfg, frequency = frequency)
  write_table(model_dt, out_path)

  # Store a compact audit row for each dataset. This is useful for checking that
  # later HMM sample sizes match the data construction step.
  sample_rows[[length(sample_rows) + 1L]] <- data.table(
    dataset = sprintf("%ss_all_buckets", frequency),
    file = out_path,
    observations = nrow(model_dt),
    stocks = uniqueN(model_dt$seccode),
    stock_days = uniqueN(model_dt$sequence_id),
    zero_return_share = mean(model_dt$r_t == 0, na.rm = TRUE),
    trade_bucket_share = mean(model_dt$number_of_trades > 0, na.rm = TRUE)
  )

  if (frequency == base_frequency) {
    # The 5-second trade-bucket dataset is a robustness input only. It is kept
    # separately from the all-bucket 5-second panel.
    trade_dt <- finalize_model_variables(aggregated, frequency, trade_buckets_only = TRUE)
    trade_path <- combined_data_file(cfg, frequency = frequency, trade_buckets_only = TRUE)
    write_table(trade_dt, trade_path)
    sample_rows[[length(sample_rows) + 1L]] <- data.table(
      dataset = sprintf("%ss_trade_buckets", frequency),
      file = trade_path,
      observations = nrow(trade_dt),
      stocks = uniqueN(trade_dt$seccode),
      stock_days = uniqueN(trade_dt$sequence_id),
      zero_return_share = mean(trade_dt$r_t == 0, na.rm = TRUE),
      trade_bucket_share = mean(trade_dt$number_of_trades > 0, na.rm = TRUE)
    )
  }
}

sample_construction <- rbindlist(sample_rows, use.names = TRUE, fill = TRUE)
write_table(sample_construction, file.path(cfg$combined_data_root, "sample_construction.csv"))

message("Done.")
message("Combined datasets folder: ", cfg$combined_data_root)
