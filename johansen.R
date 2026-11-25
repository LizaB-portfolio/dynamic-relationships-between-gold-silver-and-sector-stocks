# ============================================================
# SECTOR STOCKS: Johansen cointegration tests (Case II)
#
# Specification:
#   - Input: log-level datasets only (levels1_c ... levels7_m)
#   - VAR order in levels: 2
#   - VECM lagged differences: 1 (= VAR order - 1)
#   - Deterministic component: restricted constant in the
#     cointegration relation, no unrestricted intercept
#     (urca::ca.jo: ecdet = "const")
#   - Tests: trace and maximum eigenvalue
#
# Run from the project root (the folder containing outputs/):
#   Rscript johansen_case_ii.R
#
# Required package:
#   install.packages("urca")
# ============================================================

# --------------------------- 1. SETTINGS ---------------------
INPUT_DIR <- file.path("outputs", "data_for_models")
OUTPUT_DIR <- file.path("outputs", "econometric_tests")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")

VAR_ORDER <- 2L
VECM_LAGGED_DIFFERENCES <- VAR_ORDER - 1L
ALPHA <- 0.05
DETERMINISTIC_SPEC <- "Case II: restricted constant in cointegration relation"

GROUPS <- list(
  c = c("Gold", "Silver", "Cons_discr", "Financials",
        "Industrials", "Materials"),
  d = c("Gold", "Silver", "Cons_staples", "Health", "Utilities"),
  m = c("Gold", "Silver", "Inf_tech", "Real_estate", "Telecom", "Energy")
)

PERIODS <- 1:7
GROUP_CODES <- names(GROUPS)

# --------------------------- 2. PACKAGES ---------------------
if (!requireNamespace("urca", quietly = TRUE)) {
  stop("Package 'urca' is required. Install it with install.packages('urca').")
}

# --------------------------- 3. DIRECTORIES ------------------
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# --------------------------- 4. INPUT VALIDATION --------------
read_levels_dataset <- function(dataset_name, expected_variables) {
  file_path <- file.path(INPUT_DIR, paste0(dataset_name, ".csv"))
  if (!file.exists(file_path)) {
    stop("Input file not found: ", normalizePath(file_path, mustWork = FALSE))
  }

  df <- read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8"
  )

  if (!"Date" %in% names(df)) {
    stop("Column 'Date' is missing in: ", file_path)
  }

  dates <- as.Date(df$Date)
  if (anyNA(dates)) stop("Could not parse dates in: ", file_path)
  if (anyDuplicated(dates)) stop("Duplicate dates in: ", file_path)
  if (is.unsorted(dates, strictly = TRUE)) {
    stop("Dates must be strictly increasing in: ", file_path)
  }

  missing_variables <- setdiff(expected_variables, names(df))
  extra_variables <- setdiff(names(df), c("Date", expected_variables))
  if (length(missing_variables) > 0L || length(extra_variables) > 0L) {
    stop(
      "Unexpected columns in ", file_path,
      ". Missing: ", paste(missing_variables, collapse = ", "),
      ". Extra: ", paste(extra_variables, collapse = ", ")
    )
  }

  # Preserve the variable order specified in GROUPS.
  x_df <- df[, expected_variables, drop = FALSE]
  x_df[] <- lapply(x_df, function(z) {
    as.numeric(gsub(",", ".", as.character(z), fixed = TRUE))
  })

  x <- as.matrix(x_df)
  storage.mode(x) <- "double"
  rownames(x) <- as.character(dates)

  if (nrow(x) <= ncol(x) + VAR_ORDER + 1L) {
    stop("Too few observations for the system in ", dataset_name,
         ": n = ", nrow(x), ", variables = ", ncol(x))
  }
  if (any(!is.finite(x))) {
    stop("Missing, non-numeric, Inf or -Inf values in: ", file_path)
  }
  if (any(apply(x, 2L, stats::sd) == 0)) {
    stop("At least one variable is constant in: ", file_path)
  }

  list(data = x, dates = dates)
}

# --------------------------- 5. JOHANSEN TEST ----------------
# Extract test statistics and critical values. urca reports the
# hypotheses as row names (e.g. 'r = 0 |' or 'r <= 1 |'). We parse
# those labels and sort by r0 so that the rank decision is made
# from r0 = 0 upwards, regardless of the print order.
extract_test_table <- function(johansen_object, test_type, dataset_name, alpha = ALPHA) {
  test_stats <- as.numeric(johansen_object@teststat)
  critical_values <- as.matrix(johansen_object@cval)
  labels <- rownames(critical_values)

  if (is.null(labels) || length(test_stats) != nrow(critical_values)) {
    stop("Unexpected structure of urca::ca.jo output for ", dataset_name,
         " (", test_type, ").")
  }

  rank_values <- suppressWarnings(as.integer(
    sub(".*r\\s*(?:<=|=)\\s*([0-9]+).*", "\\1", labels, perl = TRUE)
  ))
  if (anyNA(rank_values)) {
    stop("Could not parse null-hypothesis labels from urca output for ",
         dataset_name, " (", test_type, "): ", paste(labels, collapse = "; "))
  }

  if (!"5pct" %in% colnames(critical_values)) {
    stop("The 5% critical values are missing in urca output for ", dataset_name)
  }

  out <- data.frame(
    dataset = dataset_name,
    test_type = test_type,
    rank_null = rank_values,
    null_hypothesis = labels,
    statistic = test_stats,
    critical_10pct = critical_values[, "10pct"],
    critical_5pct = critical_values[, "5pct"],
    critical_1pct = critical_values[, "1pct"],
    stringsAsFactors = FALSE
  )

  out$reject_H0_5pct <- out$statistic > out$critical_5pct
  out <- out[order(out$rank_null), , drop = FALSE]
  rownames(out) <- NULL

  # Sequential rank decision: stop at the first non-rejected H0.
  # If every tested H0 is rejected, the selected rank is k.
  first_not_rejected <- which(!out$reject_H0_5pct)[1L]
  selected_rank <- if (is.na(first_not_rejected)) {
    nrow(critical_values)
  } else {
    out$rank_null[first_not_rejected]
  }

  attr(out, "selected_rank") <- selected_rank
  out
}

run_one_dataset <- function(dataset_name, group_code, period_number) {
  expected_variables <- GROUPS[[group_code]]
  input <- read_levels_dataset(dataset_name, expected_variables)
  x <- input$data

  # K is the VAR order in levels, not the number of differenced lags.
  # ecdet = "const" places the constant inside the cointegration
  # relation (restricted constant; Case II), with no unrestricted
  # short-run intercept.
  johansen_trace <- urca::ca.jo(
    x,
    type = "trace",
    ecdet = "const",
    K = VAR_ORDER,
    spec = "longrun"
  )

  johansen_eigen <- urca::ca.jo(
    x,
    type = "eigen",
    ecdet = "const",
    K = VAR_ORDER,
    spec = "longrun"
  )

  trace_table <- extract_test_table(johansen_trace, "trace", dataset_name)
  eigen_table <- extract_test_table(johansen_eigen, "max_eigen", dataset_name)

  # The number of variables is the maximum possible cointegration rank.
  k <- ncol(x)
  trace_rank <- attr(trace_table, "selected_rank")
  eigen_rank <- attr(eigen_table, "selected_rank")
  if (trace_rank > k) trace_rank <- k
  if (eigen_rank > k) eigen_rank <- k

  metadata <- data.frame(
    dataset = dataset_name,
    period = period_number,
    group = group_code,
    n_obs = nrow(x),
    n_variables = k,
    start_date = as.character(min(input$dates)),
    end_date = as.character(max(input$dates)),
    var_order_levels = VAR_ORDER,
    vecm_lagged_differences = VECM_LAGGED_DIFFERENCES,
    deterministic_specification = DETERMINISTIC_SPEC,
    alpha = ALPHA,
    trace_rank_5pct = trace_rank,
    max_eigen_rank_5pct = eigen_rank,
    ranks_agree = identical(as.integer(trace_rank), as.integer(eigen_rank)),
    stringsAsFactors = FALSE
  )

  list(
    summary = metadata,
    details = rbind(trace_table, eigen_table)
  )
}

# --------------------------- 6. RUN ALL 21 DATASETS -----------
all_results <- list()
result_index <- 1L

for (group_code in GROUP_CODES) {
  for (period_number in PERIODS) {
    dataset_name <- paste0("levels", period_number, "_", group_code)
    message("Running Johansen tests: ", dataset_name)

    all_results[[result_index]] <- run_one_dataset(
      dataset_name = dataset_name,
      group_code = group_code,
      period_number = period_number
    )
    result_index <- result_index + 1L
  }
}

summary_table <- do.call(rbind, lapply(all_results, `[[`, "summary"))
details_table <- do.call(rbind, lapply(all_results, `[[`, "details"))
rownames(summary_table) <- NULL
rownames(details_table) <- NULL

# --------------------------- 7. SAVE RESULTS -----------------
summary_file <- file.path(TABLES_DIR, "johansen_case_ii_summary.csv")
details_file <- file.path(TABLES_DIR, "johansen_case_ii_details.csv")

write.csv(summary_table, summary_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(details_table, details_file, row.names = FALSE, fileEncoding = "UTF-8")

# --------------------------- 8. FINAL CHECKS ------------------
expected_datasets <- length(GROUPS) * length(PERIODS)
if (nrow(summary_table) != expected_datasets) {
  stop("Expected ", expected_datasets, " summary rows, got ", nrow(summary_table), ".")
}

expected_detail_rows <- sum(vapply(GROUPS, length, integer(1)) * length(PERIODS) * 2L)
if (nrow(details_table) != expected_detail_rows) {
  stop("Expected ", expected_detail_rows, " test rows, got ", nrow(details_table), ".")
}

message("\nJohansen tests completed successfully.")
message("Summary: ", normalizePath(summary_file, winslash = "/", mustWork = TRUE))
message("Details: ", normalizePath(details_file, winslash = "/", mustWork = TRUE))
message("Specification: VAR(", VAR_ORDER, ") in levels; ecdet = 'const' (Case II).")
