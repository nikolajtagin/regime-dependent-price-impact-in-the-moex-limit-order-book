###############################################################################
# MOEX Order Log data preparation
#
# What this script does:
#   1. Takes one large MOEX Full Orders Log file.
#   2. Streams rows for one selected stock without saving a filtered copy.
#   3. Reconstructs the best bid and best ask from the order log.
#   4. Aggregates the data into 5-second buckets.
#   5. Constructs model-ready variables for the hidden Markov model.
#   6. Saves only the processed model-ready data.
#
# Raw input is not changed.
###############################################################################

library(data.table)

# =========================
# 1. SETTINGS TO CHANGE
# =========================

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

split_env_vector <- function(value) {
  out <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
  out[nzchar(out)]
}

require_file <- function(path, description = "file") {
  if (!file.exists(path)) stop("Required ", description, " not found: ", path)
  invisible(normalizePath(path, mustWork = TRUE))
}

require_dir <- function(path, description = "folder") {
  if (!dir.exists(path)) stop("Required ", description, " not found: ", path)
  invisible(normalizePath(path, mustWork = TRUE))
}

code_dir <- resolve_code_dir()
project_root <- normalizePath(file.path(code_dir, ".."), mustWork = TRUE)
source(file.path(code_dir, "01_project_settings_and_functions.R"))
cfg <- read_project_settings(project_root)

raw_data_root <- cfg$raw_data_root
results_root <- cfg$results_root
daily_processed_root <- cfg$daily_processed_root

# Most user-facing switches are environment variables. This lets the same script
# run as a full production batch or as a one-stock, one-day test without editing
# code. The defaults below reproduce the paper sample.
default_tickers <- split_env_vector(Sys.getenv("MOEX_TICKERS", unset = paste(cfg$tickers, collapse = ",")))
bucket_size <- as.integer(Sys.getenv("MOEX_BUCKET_SIZE_SECONDS", unset = as.character(cfg$reconstruction_bucket_size_seconds)))
rv_window <- as.integer(Sys.getenv("MOEX_RV_WINDOW_BUCKETS", unset = as.character(cfg$rv_window_buckets)))
tod_bin_seconds <- as.integer(Sys.getenv("MOEX_TOD_BIN_SECONDS", unset = as.character(cfg$tod_bin_seconds)))
epsilon <- as.numeric(Sys.getenv("MOEX_EPSILON", unset = as.character(cfg$epsilon)))

session_start <- as.integer(Sys.getenv("MOEX_SESSION_START_SECONDS", unset = as.character(cfg$session_start_seconds)))
session_end <- as.integer(Sys.getenv("MOEX_SESSION_END_SECONDS", unset = as.character(cfg$session_end_seconds)))
open_end <- as.integer(Sys.getenv("MOEX_OPEN_END_SECONDS", unset = as.character(cfg$open_end_seconds)))
close_start <- as.integer(Sys.getenv("MOEX_CLOSE_START_SECONDS", unset = as.character(cfg$close_start_seconds)))

raw_file_for_date <- function(date_string) {
  # Raw files may be named either OrderLogYYYYMMDD or DD.MM.YYYY, and they may
  # be plain text or zipped. The first existing candidate is used.
  d <- as.Date(date_string, format = "%Y%m%d")
  dotted <- format(d, "%d.%m.%Y")
  candidates <- file.path(
    raw_data_root,
    c(
      paste0("OrderLog", date_string, ".txt"),
      paste0("OrderLog", date_string, ".zip"),
      paste0(dotted, ".txt"),
      paste0(dotted, ".zip")
    )
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop("Raw MOEX file not found for date ", date_string, " in: ", raw_data_root)
  }
  existing[1L]
}

# Stock tickers to process when process_all_selected_dates is TRUE.
tickers <- toupper(default_tickers)

# Set MOEX_PROCESS_ALL=TRUE for the production panel run. The default setting is
# production-oriented, but the environment flag makes quick smoke tests easy.
process_all_selected_dates <- tolower(Sys.getenv(
  "MOEX_PROCESS_ALL",
  unset = if (isTRUE(cfg$process_all_selected_dates)) "TRUE" else "FALSE"
)) %in% c("true", "1", "yes")

# Leave empty to use all raw files within the date range below.
dates_to_process <- character(0)

# Baseline sample window.
sample_start_date <- Sys.getenv("MOEX_START_DATE", unset = cfg$sample_start_code)
sample_end_date <- Sys.getenv("MOEX_END_DATE", unset = cfg$sample_end_code)

# Weekends are excluded by default because the baseline design is the regular
# main trading session. Set TRUE only if weekend sessions are part of the design.
include_weekends <- tolower(Sys.getenv(
  "MOEX_INCLUDE_WEEKENDS",
  unset = if (isTRUE(cfg$include_weekends)) "TRUE" else "FALSE"
)) %in% c("true", "1", "yes")

# If FALSE, a stock-day is skipped when its processed file already exists.
reprocess_existing_days <- tolower(Sys.getenv(
  "MOEX_REPROCESS_EXISTING",
  unset = if (isTRUE(cfg$reprocess_existing_days)) "TRUE" else "FALSE"
)) %in% c("true", "1", "yes")

# Stock ticker to process in single-date mode.
ticker <- toupper(Sys.getenv("MOEX_TICKER", unset = tickers[1L]))

# Trading date of the raw file, in YYYYMMDD format.
# Change default_trade_date when running this script manually for one day.
# A batch runner may also pass MOEX_TRADE_DATE through the environment.
default_trade_date <- sample_start_date
trade_date <- Sys.getenv("MOEX_TRADE_DATE", unset = default_trade_date)

# For the model file, keep the main observed session only.
# For this MOEX sample, the useful continuous session is 10:00-18:40.
use_main_session_only <- TRUE

# =========================
# 2. FOLDER SETUP
# =========================

require_dir(raw_data_root, "raw MOEX data folder")

discover_raw_dates <- function() {
  # Find every raw file that looks like a MOEX order log and convert its name to
  # YYYYMMDD. The batch loop later filters these dates to the sample window.
  raw_names <- list.files(raw_data_root, pattern = "^(OrderLog[0-9]{8}|[0-9]{2}\\.[0-9]{2}\\.[0-9]{4})\\.(txt|zip)$")
  date_codes <- sub("^OrderLog([0-9]{8})\\.(txt|zip)$", "\\1", raw_names)
  dotted <- grepl("^[0-9]{2}\\.[0-9]{2}\\.[0-9]{4}\\.(txt|zip)$", raw_names)
  date_codes[dotted] <- vapply(strsplit(sub("\\.(txt|zip)$", "", raw_names[dotted]), "\\."), function(parts) {
    paste0(parts[3L], parts[2L], parts[1L])
  }, character(1L))
  unique(sort(date_codes[grepl("^[0-9]{8}$", date_codes)]))
}

if (process_all_selected_dates && !tolower(Sys.getenv("MOEX_RUN_SINGLE_DATE", unset = "FALSE")) %in% c("true", "1", "yes")) {
  raw_dates <- discover_raw_dates()
  selected_dates <- raw_dates[raw_dates >= sample_start_date & raw_dates <= sample_end_date]
  if (!include_weekends) {
    weekdays_flag <- as.POSIXlt(as.Date(selected_dates, "%Y%m%d"))$wday %in% 1:5
    selected_dates <- selected_dates[weekdays_flag]
  }
  if (length(dates_to_process) > 0L) {
    selected_dates <- intersect(selected_dates, dates_to_process)
  }
  if (length(selected_dates) == 0L) stop("No raw dates selected for processing.")

  rscript <- file.path(R.home("bin"), "Rscript")
  script_path <- file.path(code_dir, "02_moex_orderlog_processing.R")
  message("Selected tickers: ", paste(tickers, collapse = ", "))
  message("Selected dates: ", length(selected_dates))

  for (ticker_to_run in tickers) {
    for (date_to_run in selected_dates) {
      target_file <- processed_file_for_date(cfg, ticker_to_run, date_to_run, bucket_size)
      if (file.exists(target_file) && !reprocess_existing_days) {
        message("Skipping ", ticker_to_run, " ", date_to_run, " because processed data already exists.")
        next
      }
      message("Processing ", ticker_to_run, " ", date_to_run)
      # Re-enter this same script in single-date mode. That keeps the batch
      # controller simple and makes failures point to one stock-day at a time.
      spawn_names <- c(
        "MOEX_RUN_SINGLE_DATE", "MOEX_TICKER", "MOEX_TRADE_DATE",
        "MOEX_RAW_DATA_ROOT", "MOEX_RESULTS_ROOT"
      )
      spawn_restore <- Sys.getenv(spawn_names, unset = NA_character_)
      Sys.setenv(
        MOEX_RUN_SINGLE_DATE = "TRUE",
        MOEX_TICKER          = ticker_to_run,
        MOEX_TRADE_DATE      = date_to_run,
        MOEX_RAW_DATA_ROOT   = raw_data_root,
        MOEX_RESULTS_ROOT    = results_root
      )
      status <- system2(rscript, args = script_path)
      for (nm in spawn_names) {
        if (is.na(spawn_restore[nm])) {
          Sys.unsetenv(nm)
        } else {
          v <- list(spawn_restore[nm])
          names(v) <- nm
          do.call(Sys.setenv, v)
        }
      }
      if (!identical(status, 0L)) {
        stop("Daily processing failed for ", ticker_to_run, " ", date_to_run, " with status: ", status)
      }
    }
  }

  message("Batch processing finished.")
  quit(save = "no", status = 0)
}

raw_input <- raw_file_for_date(trade_date)
require_file(raw_input, "raw MOEX file")

data_root <- daily_processed_root
dir.create(data_root, recursive = TRUE, showWarnings = FALSE)

processed_file <- processed_file_for_date(cfg, ticker, trade_date, bucket_size)
dir.create(dirname(processed_file), recursive = TRUE, showWarnings = FALSE)

# =========================
# 3. SMALL HELPER FUNCTIONS
# =========================

# MOEX time is stored as HHMMSSZZZXXX after March 2016.
# Example: 100000117495 means 10:00:00.117495.
parse_moex_time <- function(x) {
  x <- gsub("[^0-9]", "", as.character(x))
  out <- rep(NA_real_, length(x))
  ok <- nzchar(x)

  y <- x[ok]
  y <- ifelse(nchar(y) <= 6L, sprintf("%06s", y), sprintf("%012s", y))
  y <- gsub(" ", "0", y, fixed = TRUE)

  hh <- as.integer(substr(y, 1L, 2L))
  mm <- as.integer(substr(y, 3L, 4L))
  ss <- as.integer(substr(y, 5L, 6L))

  frac_str <- ifelse(nchar(y) > 6L, substr(y, 7L, nchar(y)), "")
  frac <- ifelse(nzchar(frac_str), as.numeric(paste0("0.", frac_str)), 0)

  valid <- !is.na(hh) & !is.na(mm) & !is.na(ss) &
    hh >= 0L & hh <= 23L & mm >= 0L & mm <= 59L & ss >= 0L & ss <= 59L

  vals <- hh * 3600 + mm * 60 + ss + frac
  vals[!valid] <- NA_real_
  out[ok] <- vals
  out
}

# This is used as a key for price-level depth maps.
price_key <- function(price) {
  format(price, scientific = FALSE, trim = TRUE, digits = 15)
}

# =========================
# 4. STREAM RAW FILE TO ONE STOCK
# =========================

message("Step 1: reading ticker ", ticker, " from raw file")

expected_names <- c(
  "NO", "SECCODE", "BUYSELL", "TIME", "ORDERNO",
  "ACTION", "PRICE", "VOLUME", "TRADENO", "TRADEPRICE"
)

build_raw_stream_command <- function(input, ticker) {
  if (!nzchar(Sys.which("awk"))) {
    stop("awk is required for streaming ticker rows without saving a filtered file.")
  }

  # The raw files can be very large. Streaming only one ticker with awk avoids
  # creating temporary filtered files and keeps memory use predictable.
  awk_script <- 'BEGIN {FS=","; OFS=","} NR==1 || $2==ticker {print}'
  awk_part <- sprintf("awk -v ticker=%s %s", shQuote(ticker), shQuote(awk_script))

  if (grepl("\\.zip$", input, ignore.case = TRUE)) {
    if (!nzchar(Sys.which("unzip"))) {
      stop("unzip is required to stream compressed raw MOEX files.")
    }
    sprintf("unzip -p %s | %s", shQuote(input), awk_part)
  } else {
    sprintf("%s %s", awk_part, shQuote(input))
  }
}

dt <- fread(
  cmd = build_raw_stream_command(raw_input, ticker),
  colClasses = list(
    integer = c("NO", "ACTION"),
    character = c("SECCODE", "BUYSELL", "TIME", "ORDERNO", "TRADENO"),
    numeric = c("PRICE", "VOLUME", "TRADEPRICE")
  ),
  na.strings = c("", "NA")
)

if (!identical(names(dt), expected_names)) {
  # Some raw exports arrive without the exact expected header. If the column
  # count still matches, rename defensively; otherwise stop rather than silently
  # reading the wrong fields.
  if (length(names(dt)) == length(expected_names)) {
    setnames(dt, expected_names)
  } else {
    stop("Unexpected column structure in raw MOEX file.")
  }
}

if (nrow(dt) == 0L) stop("No rows found for ticker ", ticker, " on ", trade_date)

dt[, seconds_from_midnight := parse_moex_time(TIME)]
setorder(dt, NO)

# =========================
# 6. RECONSTRUCT TOP OF BOOK
# =========================

message("Step 3: reconstructing top of book and aggregating 5-second buckets")

# The order log contains events for individual orders. To reconstruct the best
# bid and ask, the script has to keep track of which limit orders are currently
# active in the visible book.
#
# The three active_* environments are indexed by ORDERNO:
#   active_side[[ORDERNO]]   = "B" or "S";
#   active_price[[ORDERNO]]  = order limit price;
#   active_volume[[ORDERNO]] = currently remaining visible volume.
#
# Environments are used here because repeated lookup by ORDERNO is much faster
# than repeatedly searching a large data.table inside the event loop.
active_side <- new.env(hash = TRUE, parent = emptyenv())
active_price <- new.env(hash = TRUE, parent = emptyenv())
active_volume <- new.env(hash = TRUE, parent = emptyenv())

# These two environments aggregate active order volume by price level. They make
# it possible to read the depth available at the current best bid and best ask.
bid_depth <- new.env(hash = TRUE, parent = emptyenv())
ask_depth <- new.env(hash = TRUE, parent = emptyenv())

# MOEX trade rows may contain both sides of the same trade. The seen_trades
# environment ensures that each unique TRADENO contributes once to trade volume,
# signed order flow q_t, and the number of trades in the 5-second bucket.
seen_trades <- new.env(hash = TRUE, parent = emptyenv())

# The best quotes are cached. Most events do not change the best quote, so
# caching avoids scanning all price levels after every event. When the whole
# best level disappears, recompute_best() scans the remaining levels.
best_bid <- NA_real_
best_ask <- NA_real_

# Used by trade-sign proxies. The baseline uses a quote-first rule and falls
# back to the tick rule when the trade prints at the midpoint or quotes are not
# available. BUYSELL-based signs are saved separately for validation because
# the full order log often contains two rows per TRADENO, one for each side.
last_trade_price <- NA_real_
last_trade_sign <- 0

level_env <- function(side) {
  if (identical(side, "B")) bid_depth else ask_depth
}

recompute_best <- function(side) {
  env <- level_env(side)
  keys <- ls(env, all.names = TRUE)
  if (length(keys) == 0L) {
    return(NA_real_)
  }
  vals <- as.numeric(keys)
  if (identical(side, "B")) max(vals) else min(vals)
}

# Read total visible volume at one price level. If the price level is absent
# from the reconstructed book, depth is zero.
get_level_depth <- function(side, price) {
  if (is.na(price)) {
    return(0)
  }
  env <- level_env(side)
  key <- price_key(price)
  if (exists(key, envir = env, inherits = FALSE)) get(key, envir = env) else 0
}

top_k_depth <- function(side, k = 5L) {
  env <- level_env(side)
  keys <- ls(env, all.names = TRUE)
  if (length(keys) == 0L) {
    return(0)
  }

  prices <- suppressWarnings(as.numeric(keys))
  valid <- is.finite(prices)
  if (!any(valid)) {
    return(0)
  }

  prices <- prices[valid]
  keys <- keys[valid]
  order_idx <- if (identical(side, "B")) order(prices, decreasing = TRUE) else order(prices)
  top_keys <- keys[order_idx][seq_len(min(k, length(order_idx)))]

  sum(vapply(top_keys, function(key) get(key, envir = env), numeric(1L)), na.rm = TRUE)
}

topbook_snapshot <- function() {
  list(
    bid = best_bid,
    bid_depth = get_level_depth("B", best_bid),
    ask = best_ask,
    ask_depth = get_level_depth("S", best_ask)
  )
}

topbook_ofi_contribution <- function(before, after) {
  if (is.na(before$bid) || is.na(before$ask) || is.na(after$bid) || is.na(after$ask)) {
    return(0)
  }

  # Cont-Kukanov-Stoikov OFI at the best quotes. Improvements on the bid
  # and withdrawals from the ask are positive buy-side pressure; weakening
  # bid quotes and replenishment on the ask are negative sell-side pressure.
  bid_part <- 0
  ask_part <- 0

  if (after$bid >= before$bid) bid_part <- bid_part + after$bid_depth
  if (after$bid <= before$bid) bid_part <- bid_part - before$bid_depth
  if (after$ask <= before$ask) ask_part <- ask_part - after$ask_depth
  if (after$ask >= before$ask) ask_part <- ask_part + before$ask_depth

  bid_part + ask_part
}

# Add or remove volume at a price level and update the cached best bid/ask.
# Positive delta means new visible volume; negative delta means volume removed
# because of cancellation or execution.
add_level_depth <- function(side, price, delta) {
  key <- price_key(price)
  env <- level_env(side)
  old <- if (exists(key, envir = env, inherits = FALSE)) get(key, envir = env) else 0
  new <- old + delta

  if (new <= 0) {
    if (exists(key, envir = env, inherits = FALSE)) rm(list = key, envir = env)
  } else {
    assign(key, new, envir = env)
  }

  if (identical(side, "B")) {
    if (is.na(best_bid) || (delta > 0 && price > best_bid)) best_bid <<- price
    if (!is.na(best_bid) && price == best_bid && new <= 0) best_bid <<- recompute_best("B")
  } else if (identical(side, "S")) {
    if (is.na(best_ask) || (delta > 0 && price < best_ask)) best_ask <<- price
    if (!is.na(best_ask) && price == best_ask && new <= 0) best_ask <<- recompute_best("S")
  }
}

# Remove an order completely from the active book and from its price-level depth.
remove_order <- function(orderno) {
  if (!exists(orderno, envir = active_side, inherits = FALSE)) {
    return(FALSE)
  }
  side <- get(orderno, envir = active_side)
  price <- get(orderno, envir = active_price)
  volume <- get(orderno, envir = active_volume)
  add_level_depth(side, price, -volume)
  rm(list = orderno, envir = active_side)
  rm(list = orderno, envir = active_price)
  rm(list = orderno, envir = active_volume)
  TRUE
}

# Change the remaining volume of an active order and adjust the corresponding
# price-level depth by the volume difference.
set_order_volume <- function(orderno, new_volume) {
  side <- get(orderno, envir = active_side)
  price <- get(orderno, envir = active_price)
  old_volume <- get(orderno, envir = active_volume)
  delta <- new_volume - old_volume
  assign(orderno, new_volume, envir = active_volume)
  add_level_depth(side, price, delta)
}

# Bucket-level variables are accumulated while reading events. A bucket is
# finalized when the next event belongs to a later bucket; therefore the stored
# best bid, best ask, spread, midprice, and depth describe the reconstructed
# book state at the end of that 5-second interval.
bucket_rows <- vector("list", 10000L)
bucket_count <- 0L
current_bucket <- NA_real_

event_count <- add_count <- cancel_count <- trade_count <- 0L
trade_volume <- signed_trade_volume <- abs_signed_trade_volume <- 0
signed_trade_volume_tick <- signed_trade_volume_buysell <- signed_trade_volume_buysell_opposite <- 0
event_level_ofi <- 0
buy_trade_volume <- sell_trade_volume <- 0

ensure_capacity <- function() {
  if (bucket_count >= length(bucket_rows)) length(bucket_rows) <<- length(bucket_rows) * 2L
}

reset_bucket_aggregates <- function() {
  event_count <<- add_count <<- cancel_count <<- trade_count <<- 0L
  trade_volume <<- signed_trade_volume <<- abs_signed_trade_volume <<- 0
  signed_trade_volume_tick <<- signed_trade_volume_buysell <<- signed_trade_volume_buysell_opposite <<- 0
  event_level_ofi <<- 0
  buy_trade_volume <<- sell_trade_volume <<- 0
}

finalize_bucket <- function(bucket_start) {
  ensure_capacity()

  # Depth is measured at the best quoted levels, not across the whole book.
  bid_depth_best <- get_level_depth("B", best_bid)
  ask_depth_best <- get_level_depth("S", best_ask)
  bid_depth_top5 <- top_k_depth("B", 5L)
  ask_depth_top5 <- top_k_depth("S", 5L)
  has_quotes <- !is.na(best_bid) && !is.na(best_ask) && best_bid > 0 && best_ask > 0

  # Crossed and locked quotes are kept as diagnostics. Crossed quotes are not
  # used for the model-ready sample because spread and midpoint are then not
  # economically meaningful for this baseline construction.
  crossed <- has_quotes && best_ask < best_bid
  locked <- has_quotes && best_ask == best_bid
  valid <- has_quotes && !crossed

  spread <- if (has_quotes) best_ask - best_bid else NA_real_
  midprice <- if (has_quotes) (best_ask + best_bid) / 2 else NA_real_
  depth <- if (has_quotes) bid_depth_best + ask_depth_best else NA_real_
  queue_imbalance <- if (has_quotes && depth > 0) (bid_depth_best - ask_depth_best) / depth else NA_real_
  top5_depth <- if (has_quotes) bid_depth_top5 + ask_depth_top5 else NA_real_
  top5_queue_imbalance <- if (has_quotes && top5_depth > 0) (bid_depth_top5 - ask_depth_top5) / top5_depth else NA_real_

  bucket_count <<- bucket_count + 1L
  bucket_rows[[bucket_count]] <<- data.table(
    seccode = ticker,
    bucket_start = bucket_start,
    bucket_end = bucket_start + bucket_size,
    seconds_from_midnight = bucket_start,
    best_bid = best_bid,
    best_ask = best_ask,
    best_bid_depth = bid_depth_best,
    best_ask_depth = ask_depth_best,
    top5_bid_depth = bid_depth_top5,
    top5_ask_depth = ask_depth_top5,
    top5_depth = top5_depth,
    spread = spread,
    midprice = midprice,
    depth = depth,
    number_of_events = event_count,
    number_of_order_adds = add_count,
    number_of_cancellations = cancel_count,
    number_of_trades = trade_count,
    trade_volume = trade_volume,
    signed_trade_volume_proxy = signed_trade_volume,
    signed_trade_volume_tick = signed_trade_volume_tick,
    signed_trade_volume_buysell = signed_trade_volume_buysell,
    signed_trade_volume_buysell_opposite = signed_trade_volume_buysell_opposite,
    abs_signed_trade_volume_proxy = abs_signed_trade_volume,
    buy_trade_volume_proxy = buy_trade_volume,
    sell_trade_volume_proxy = sell_trade_volume,
    event_level_ofi = event_level_ofi,
    queue_imbalance = queue_imbalance,
    top5_queue_imbalance = top5_queue_imbalance,
    crossed_quote_flag = crossed,
    locked_quote_flag = locked,
    valid_topbook_flag = valid
  )
}

unmatched_cancellations <- 0L
unmatched_trades <- 0L
duplicate_placements <- 0L
price_zero_rows <- 0L
cancellation_volume_matches <- 0L
cancellation_volume_mismatches <- 0L
cancellation_volume_missing <- 0L
sign_validation_rows <- vector("list", 10000L)
sign_validation_count <- 0L

time_period_from_seconds <- function(x) {
  if (is.na(x)) {
    return(NA_character_)
  }
  if (x < open_end) {
    return("open")
  }
  if (x >= close_start) {
    return("close")
  }
  "midday"
}

append_sign_validation <- function(bucket_start, quote_direction, tick_direction, buysell_direction, signed_direction) {
  sign_validation_count <<- sign_validation_count + 1L
  if (sign_validation_count > length(sign_validation_rows)) {
    length(sign_validation_rows) <<- length(sign_validation_rows) * 2L
  }

  sign_validation_rows[[sign_validation_count]] <<- data.table(
    bucket_start = bucket_start,
    time_period = time_period_from_seconds(bucket_start),
    quote_direction = quote_direction,
    tick_direction = tick_direction,
    buysell_direction = buysell_direction,
    buysell_opposite_direction = -buysell_direction,
    baseline_direction = signed_direction
  )
}

# Pull columns into vectors because this is faster inside the event loop than
# repeated dt[i, column] access.
NO <- dt$NO
side <- dt$BUYSELL
time_sec <- dt$seconds_from_midnight
orderno <- dt$ORDERNO
action <- dt$ACTION
price <- dt$PRICE
volume <- dt$VOLUME
tradeprice <- dt$TRADEPRICE
tradeno <- dt$TRADENO

for (i in seq_len(nrow(dt))) {
  if (is.na(time_sec[i])) next

  # Assign the event to a 5-second interval. If there are empty intervals between
  # two observed events, they are still written using the latest reconstructed
  # book state and zero event/trade counts.
  b <- floor(time_sec[i] / bucket_size) * bucket_size
  if (is.na(current_bucket)) current_bucket <- b

  while (current_bucket < b) {
    finalize_bucket(current_bucket)
    reset_bucket_aggregates()
    current_bucket <- current_bucket + bucket_size
  }

  event_count <- event_count + 1L
  if (!is.na(price[i]) && price[i] == 0 && !is.na(volume[i]) && volume[i] > 0) {
    price_zero_rows <- price_zero_rows + 1L
  }

  ord <- orderno[i]
  if (is.na(ord) || ord == "") next

  topbook_before_event <- topbook_snapshot()

  if (action[i] == 1L) {
    # ACTION=1 means a new limit order is placed. Market orders have PRICE=0 and
    # are not part of the visible limit order book, so they are not added here.
    add_count <- add_count + 1L
    if (!is.na(price[i]) && !is.na(volume[i]) && price[i] > 0 && volume[i] > 0 &&
      side[i] %in% c("B", "S")) {
      if (exists(ord, envir = active_side, inherits = FALSE)) {
        duplicate_placements <- duplicate_placements + 1L
        remove_order(ord)
      }
      assign(ord, side[i], envir = active_side)
      assign(ord, price[i], envir = active_price)
      assign(ord, volume[i], envir = active_volume)
      add_level_depth(side[i], price[i], volume[i])
    }
  } else if (action[i] == 0L) {
    # ACTION=0 means cancellation. In the MOEX layout, VOLUME is the remaining
    # visible volume of the cancelled order. The reconstructed active order is
    # therefore removed completely, and VOLUME is retained for audit checks.
    cancel_count <- cancel_count + 1L
    found <- exists(ord, envir = active_side, inherits = FALSE)
    if (!found) {
      unmatched_cancellations <- unmatched_cancellations + 1L
    } else {
      old_volume <- get(ord, envir = active_volume)
      if (is.na(volume[i])) {
        cancellation_volume_missing <- cancellation_volume_missing + 1L
      } else if (abs(old_volume - volume[i]) <= 1e-9) {
        cancellation_volume_matches <- cancellation_volume_matches + 1L
      } else {
        cancellation_volume_mismatches <- cancellation_volume_mismatches + 1L
      }
      remove_order(ord)
    }
  } else if (action[i] == 2L) {
    # ACTION=2 means trade. The resting order volume is reduced. MOEX can
    # report both sides of the same trade with the same TRADENO, so bucket-level
    # trade variables are counted only once per trade id.
    v <- if (is.na(volume[i])) 0 else volume[i]
    trade_id <- if (!is.na(tradeno[i]) && tradeno[i] != "") tradeno[i] else paste0("row_", NO[i])

    if (!exists(trade_id, envir = seen_trades, inherits = FALSE)) {
      assign(trade_id, TRUE, envir = seen_trades)
      trade_count <- trade_count + 1L
      trade_volume <- trade_volume + v

      # Baseline signed order-flow proxy:
      #   1. Quote rule first: trades above the midpoint are buyer-initiated,
      #      trades below the midpoint are seller-initiated.
      #   2. If the trade prints at the midpoint or quotes are unavailable, use
      #      the tick rule based on the previous trade price.
      #   3. If the tick is zero, reuse the previous non-zero tick sign.
      #
      # BUYSELL-based signs are saved as robustness variables. In this sample a
      # TRADENO usually has two rows, so BUYSELL is not used as the baseline
      # aggressor-side flag without validation.
      tp <- tradeprice[i]
      quote_direction <- 0
      tick_direction <- 0
      signed_direction <- 0
      if (!is.na(tp)) {
        if (!is.na(best_bid) && !is.na(best_ask)) {
          current_mid <- (best_bid + best_ask) / 2
          quote_direction <- sign(tp - current_mid)
        }

        if (!is.na(last_trade_price)) {
          tick_direction <- sign(tp - last_trade_price)
          if (tick_direction == 0) tick_direction <- last_trade_sign
        }

        signed_direction <- if (quote_direction != 0) quote_direction else tick_direction
        last_trade_price <- tp
        if (tick_direction != 0) last_trade_sign <- tick_direction
      }

      buysell_direction <- if (side[i] == "B") 1 else if (side[i] == "S") -1 else 0
      append_sign_validation(b, quote_direction, tick_direction, buysell_direction, signed_direction)
      signed <- signed_direction * v
      signed_trade_volume <- signed_trade_volume + signed
      signed_trade_volume_tick <- signed_trade_volume_tick + tick_direction * v
      signed_trade_volume_buysell <- signed_trade_volume_buysell + buysell_direction * v
      signed_trade_volume_buysell_opposite <- signed_trade_volume_buysell_opposite - buysell_direction * v
      abs_signed_trade_volume <- abs_signed_trade_volume + abs(signed)
      if (signed_direction > 0) buy_trade_volume <- buy_trade_volume + v
      if (signed_direction < 0) sell_trade_volume <- sell_trade_volume + v
    }

    found <- exists(ord, envir = active_side, inherits = FALSE)
    if (!found) {
      unmatched_trades <- unmatched_trades + 1L
    } else {
      old_volume <- get(ord, envir = active_volume)
      new_volume <- old_volume - v
      if (new_volume > 0) {
        set_order_volume(ord, new_volume)
      } else {
        remove_order(ord)
      }
    }
  }

  topbook_after_event <- topbook_snapshot()
  event_level_ofi <- event_level_ofi + topbook_ofi_contribution(topbook_before_event, topbook_after_event)
}

if (!is.na(current_bucket)) finalize_bucket(current_bucket)
bucket_dt <- rbindlist(bucket_rows[seq_len(bucket_count)], use.names = TRUE, fill = TRUE)

pad_regular_bucket_grid <- function(x) {
  full_grid <- data.table(
    seccode = ticker,
    bucket_start = seq(session_start, session_end - bucket_size, by = bucket_size)
  )
  full_grid[, `:=`(
    bucket_end = bucket_start + bucket_size,
    seconds_from_midnight = bucket_start
  )]

  x <- merge(full_grid, x, by = c("seccode", "bucket_start", "bucket_end", "seconds_from_midnight"), all.x = TRUE)
  setorder(x, bucket_start)

  quote_columns <- c(
    "best_bid", "best_ask", "best_bid_depth", "best_ask_depth",
    "top5_bid_depth", "top5_ask_depth", "top5_depth",
    "spread", "midprice", "depth", "queue_imbalance", "top5_queue_imbalance"
  )
  quote_columns <- intersect(quote_columns, names(x))
  for (col in quote_columns) {
    set(x, j = col, value = nafill(x[[col]], type = "locf"))
  }

  count_columns <- c(
    "number_of_events", "number_of_order_adds", "number_of_cancellations",
    "number_of_trades", "trade_volume", "signed_trade_volume_proxy",
    "signed_trade_volume_tick", "signed_trade_volume_buysell",
    "signed_trade_volume_buysell_opposite", "abs_signed_trade_volume_proxy",
    "buy_trade_volume_proxy", "sell_trade_volume_proxy", "event_level_ofi"
  )
  count_columns <- intersect(count_columns, names(x))
  for (col in count_columns) {
    x[is.na(get(col)), (col) := 0]
  }

  x[, crossed_quote_flag := fifelse(!is.na(best_bid) & !is.na(best_ask), best_ask < best_bid, FALSE)]
  x[, locked_quote_flag := fifelse(!is.na(best_bid) & !is.na(best_ask), best_ask == best_bid, FALSE)]
  x[, valid_topbook_flag := !is.na(best_bid) & !is.na(best_ask) & best_bid > 0 & best_ask > 0 & best_ask > best_bid]

  x
}

bucket_dt <- pad_regular_bucket_grid(bucket_dt)
expected_bucket_count <- as.integer((session_end - session_start) / bucket_size)
if (nrow(bucket_dt) != expected_bucket_count) {
  stop("Padded bucket grid has unexpected length: ", nrow(bucket_dt), " instead of ", expected_bucket_count)
}
if (bucket_dt[, anyDuplicated(bucket_start)] > 0L) {
  stop("Padded bucket grid contains duplicate bucket_start values.")
}

# =========================
# 7. CONSTRUCT MODEL VARIABLES
# =========================

message("Step 4: constructing model-ready variables")

model_base <- copy(bucket_dt)

# The raw file contains events outside the main continuous trading session. For
# the baseline model, only the main observed session is used, because liquidity
# formation outside regular trading can follow a different mechanism.
if (use_main_session_only) {
  model_base <- model_base[seconds_from_midnight >= session_start & seconds_from_midnight < session_end]
}

rows_before_cleaning <- nrow(model_base)

# Recompute the main quote variables from the reconstructed best bid and ask.
# This keeps the model inputs internally consistent even if diagnostic columns
# from bucket_dt are later changed or extended.
model_base[, spread := best_ask - best_bid]
model_base[, midprice := (best_ask + best_bid) / 2]
model_base[, depth := best_bid_depth + best_ask_depth]

# Relative spread is scale-free: a spread of 1 ruble means something different
# for a 10-ruble stock and a 300-ruble stock. The logarithm is used because
# spreads are positive and often right-skewed.
model_base[, relative_spread := spread / midprice]
model_base[, log_relative_spread := NA_real_]
model_base[relative_spread > 0, log_relative_spread := log(relative_spread)]

# Depth is also positive and skewed, so the model uses log depth before removing
# the intraday pattern.
model_base[, log_depth := NA_real_]
model_base[depth > 0, log_depth := log(depth)]

# Midprice returns are measured in basis points. This is the return variable r_t
# in the emission distribution of the HMM.
model_base[, r_t := 10000 * (log(midprice) - shift(log(midprice), 1L))]

# q_t is the signed order-flow proxy entering the price-impact part of the
# emission equation. abs_q_t is used for transition probabilities, where order-
# flow pressure matters by intensity rather than direction.
model_base[, q_t := signed_trade_volume_proxy]
model_base[, q_t_tick := signed_trade_volume_tick]
model_base[, q_t_buysell := signed_trade_volume_buysell]
model_base[, q_t_buysell_opposite := signed_trade_volume_buysell_opposite]
model_base[, abs_q_t := abs(q_t)]
model_base[, abs_event_level_ofi := abs(event_level_ofi)]

# Realized volatility is computed only from lagged returns. The current return
# is not included, so RV_t is predetermined relative to bucket t.
model_base[, RV_t := frollsum(
  shift(r_t^2, 1L),
  n = rv_window,
  align = "right",
  fill = NA_real_,
  na.rm = TRUE,
  partial = TRUE
)]
model_base[, effective_rv_window := frollsum(
  as.integer(!is.na(shift(r_t^2, 1L))),
  n = rv_window,
  align = "right",
  fill = NA_integer_,
  na.rm = TRUE,
  partial = TRUE
)]
model_base[, complete_rv_window_flag := effective_rv_window >= rv_window]

# The time-of-day bin tau_t is used to remove regular intraday patterns. With
# tod_bin_seconds = 300, each bin is five minutes long.
model_base[, tod_bin := floor(seconds_from_midnight / tod_bin_seconds) * tod_bin_seconds]

# Provisional daily de-seasonalized spread and depth:
#   spread_tilde = log(relative spread) - average log(relative spread) in tau_t;
#   depth_tilde  = log(depth)           - average log(depth) in tau_t.
#
# This removes the average intraday shape separately for each time-of-day bin
# inside this daily file. When several days are combined, the final modelling
# dataset recomputes these transformed variables using pooled intraday profiles.
model_base[, spread_tilde := log_relative_spread - mean(log_relative_spread, na.rm = TRUE), by = tod_bin]
model_base[, depth_tilde := log_depth - mean(log_depth, na.rm = TRUE), by = tod_bin]

# Realized volatility is transformed in two steps:
#   1. log(RV_t + epsilon), because realized volatility is non-negative and
#      strongly skewed;
#   2. standardization within the time-of-day bin, so the resulting rv_tilde is
#      measured relative to the usual volatility level at that time of day.
model_base[, log_RV := log(RV_t + epsilon)]
model_base[, rv_tod_mean := mean(log_RV, na.rm = TRUE), by = tod_bin]
model_base[, rv_tod_sd := sd(log_RV, na.rm = TRUE), by = tod_bin]
model_base[, rv_tilde := fifelse(
  !is.na(rv_tod_sd) & rv_tod_sd > 0,
  (log_RV - rv_tod_mean) / rv_tod_sd,
  fifelse(!is.na(log_RV), 0, NA_real_)
)]

# Standardize absolute order-flow pressure within each time-of-day bin. This is
# useful as a transition covariate because raw trading intensity also has an
# intraday pattern. The combined-data script recomputes this standardization on
# the pooled sample used for estimation.
model_base[, abs_q_tilde := {
  s <- sd(abs_q_t, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, .N) else (abs_q_t - mean(abs_q_t, na.rm = TRUE)) / s
}, by = tod_bin]

# Keep only rows where the reconstructed quote state can support the model
# variables. Crossed quotes, locked quotes, and non-positive spread/midprice/depth
# observations are removed from the model-ready sample.
model_base[, clean_model_row := valid_topbook_flag == TRUE &
  best_bid > 0 &
  best_ask > 0 &
  best_ask >= best_bid &
  midprice > 0 &
  relative_spread > 0 &
  depth > 0 &
  !is.na(r_t)]

model_dt <- model_base[clean_model_row == TRUE]

if (nrow(model_dt) == 0L) stop("Model-ready data is empty after cleaning.")
if (model_dt[, any(depth <= 0 | midprice <= 0 | spread <= 0, na.rm = TRUE)]) {
  stop("Model-ready data contains non-positive depth, midprice, or spread.")
}
if (model_dt[, any(best_ask <= best_bid, na.rm = TRUE)]) {
  stop("Model-ready data contains locked or crossed quotes.")
}

keep_cols <- c(
  "seccode", "bucket_start", "bucket_end", "seconds_from_midnight",
  "best_bid", "best_ask", "best_bid_depth", "best_ask_depth",
  "top5_bid_depth", "top5_ask_depth", "top5_depth",
  "spread", "midprice", "depth", "r_t",
  "relative_spread", "log_relative_spread", "spread_tilde",
  "q_t", "q_t_tick", "q_t_buysell", "q_t_buysell_opposite",
  "abs_q_t", "abs_q_tilde", "event_level_ofi", "abs_event_level_ofi",
  "queue_imbalance", "top5_queue_imbalance",
  "RV_t", "effective_rv_window", "complete_rv_window_flag",
  "rv_tilde", "depth_tilde", "tod_bin",
  "valid_topbook_flag", "crossed_quote_flag", "locked_quote_flag",
  "number_of_events", "number_of_trades", "trade_volume",
  "signed_trade_volume_proxy", "signed_trade_volume_tick",
  "signed_trade_volume_buysell", "signed_trade_volume_buysell_opposite",
  "abs_signed_trade_volume_proxy"
)

fwrite(model_dt[, ..keep_cols], processed_file)

message("Done.")
message("Processed file: ", processed_file)
message("Ticker rows read: ", format(nrow(dt), big.mark = ","))
message("Model rows saved: ", format(nrow(model_dt), big.mark = ","))
