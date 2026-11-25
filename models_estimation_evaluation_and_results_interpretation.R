# ============================================================
# SECTOR STOCKS: VAR / VECM, diagnostics, IRF and FEVD
#
# Input:
#   outputs/data_for_models/period*_*.csv
#   outputs/data_for_models/levels*_*.csv
#
# Run from the project root (the folder containing outputs/).
# ============================================================


# --------------------------- 0. SETTINGS ---------------------

INPUT_DIR  <- file.path("outputs", "data_for_models")
OUTPUT_DIR <- file.path("outputs", "R_models")
PLOTS_DIR  <- file.path(OUTPUT_DIR, "figures")
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
MODELS_DIR <- file.path(OUTPUT_DIR, "model_objects")

HORIZON <- 20L

BOOT_RUNS <- 500L
BOOT_SEED <- 20261117L
BOOT_CI   <- 0.95
# Bootstrap settings for GFEVD
GFEVD_BOOT_RUNS <- 500L
GFEVD_BOOT_CI <- 0.95

# NULL: automatic rule based on sample size.
# Otherwise, use the specified fixed block length.
GFEVD_BLOCK_LENGTH <- NULL

PORTMANTEAU_LAGS <- 16L
ARCH_LAGS <- 5L

GROUPS <- c("c", "d", "m")
PERIODS <- 1:7

ALL_DATASETS <- unlist(
  lapply(GROUPS, function(g) paste0("period", PERIODS, "_", g))
)

# Only these three models are estimated as VECM.
VECM_DATASETS <- c(
  "levels7_c",
  "levels4_d",
  "levels6_d"
)

# All remaining 18 models are estimated as VAR(1).
VAR_DATASETS <- setdiff(
  ALL_DATASETS,
  sub("^levels", "period", VECM_DATASETS)
)

# VECM specification:
# r = 1
# one lagged difference = VAR(2) in levels
# Case II: restricted constant in the cointegration relation.
VECM_RANK <- 1L
VECM_LAGGED_DIFFERENCES <- 1L
VECM_VAR_ORDER <- VECM_LAGGED_DIFFERENCES + 1L


# --------------------------- 1. PACKAGES ---------------------

required_packages <- c(
  "vars",
  "tsDyn",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(vars)
  library(tsDyn)
  library(ggplot2)
})

# The conversion function is internal to tsDyn.
# Its availability and output structure may depend on package version.
if (!exists(
  "vec2var.tsDyn",
  envir = asNamespace("tsDyn"),
  inherits = FALSE
)) {
  stop(
    "The installed tsDyn version does not contain vec2var.tsDyn(). ",
    "Check the package version and its documentation."
  )
}

for (directory in c(
  OUTPUT_DIR,
  PLOTS_DIR,
  TABLES_DIR,
  MODELS_DIR
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# --------------------------- 2. INPUT CHECKS -----------------

read_model_data <- function(dataset_name) {
  
  file_path <- file.path(
    INPUT_DIR,
    paste0(dataset_name, ".csv")
  )
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  df <- read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  if (!"Date" %in% names(df)) {
    stop("Column Date is missing in: ", file_path)
  }
  
  dates <- as.Date(df$Date)
  
  if (anyNA(dates)) {
    stop("Could not parse Date in: ", file_path)
  }
  
  if (anyDuplicated(dates)) {
    stop("Duplicate dates in: ", file_path)
  }
  
  if (is.unsorted(dates, strictly = TRUE)) {
    stop("Dates are not strictly increasing in: ", file_path)
  }
  
  x <- df[
    ,
    setdiff(names(df), "Date"),
    drop = FALSE
  ]
  
  if (ncol(x) < 2L) {
    stop("At least two endogenous variables are required: ", file_path)
  }
  
  x[] <- lapply(
    x,
    function(z) {
      as.numeric(
        gsub(",", ".", as.character(z), fixed = TRUE)
      )
    }
  )
  
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  rownames(x) <- as.character(dates)
  
  if (any(!is.finite(x))) {
    stop("NA/Inf/non-numeric values in: ", file_path)
  }
  
  if (any(apply(x, 2, sd) == 0)) {
    stop("A constant variable was found in: ", file_path)
  }
  
  x
}


# --------------------------- 3. MODEL ESTIMATION --------------

# VAR(1) is estimated on log returns.
fit_var1 <- function(y) {
  
  vars::VAR(
    as.data.frame(y),
    p = 1,
    type = "const"
  )
}


# VECM is estimated on log price levels.
#
# tsDyn::VECM(lag = 1) means one lagged difference.
# This corresponds to VAR(2) in levels.
#
# Case II:
# include = "none"  -> no unrestricted intercept
# LRinclude = "const" -> constant restricted to the
#                        cointegration relation.

fit_vecm <- function(y) {
  
  tsDyn::VECM(
    data = y,
    lag = VECM_LAGGED_DIFFERENCES,
    r = VECM_RANK,
    include = "none",
    estim = "ML",
    LRinclude = "const"
  )
}


# Extract parameters from a VAR model.
extract_var_parameters <- function(model, y) {
  
  A_raw <- vars::Acoef(model)
  
  if (!is.list(A_raw) || length(A_raw) < 1L) {
    stop("Unexpected object returned by vars::Acoef().")
  }
  
  A_list <- lapply(A_raw, as.matrix)
  
  k <- ncol(y)
  
  if (any(vapply(
    A_list,
    function(A) !all(dim(A) == c(k, k)),
    logical(1)
  ))) {
    stop("VAR coefficient matrices have unexpected dimensions.")
  }
  
  E <- as.matrix(stats::residuals(model))
  
  if (ncol(E) != k) {
    stop("VAR residual matrix has unexpected dimensions.")
  }
  
  # ML covariance estimate, used consistently in IRF calculations.
  Sigma <- crossprod(E) / nrow(E)
  
  p <- length(A_list)
  
  Y_target <- y[(p + 1L):nrow(y), , drop = FALSE]
  
  fitted_without_const <- matrix(
    0,
    nrow = nrow(Y_target),
    ncol = k
  )
  
  for (i in seq_len(p)) {
    
    lagged_y <- y[
      (p + 1L - i):(nrow(y) - i),
      ,
      drop = FALSE
    ]
    
    fitted_without_const <- fitted_without_const +
      lagged_y %*% t(A_list[[i]])
  }
  
  intercept <- colMeans(
    Y_target - fitted_without_const
  )
  
  list(
    A_list = A_list,
    A1 = A_list[[1L]],
    Sigma = Sigma,
    residuals = E,
    intercept = intercept,
    p = p
  )
}


# Extract parameters from a VECM and convert it to its
# equivalent VAR representation in levels.
extract_vecm_parameters <- function(model) {
  
  vec2var_model <- tsDyn:::vec2var.tsDyn(model)
  
  A_list <- lapply(
    vec2var_model$A,
    as.matrix
  )
  
  k <- model$k
  p <- vec2var_model$p
  
  if (
    length(A_list) != p ||
    any(vapply(
      A_list,
      function(A) !all(dim(A) == c(k, k)),
      logical(1)
    ))
  ) {
    stop("Unexpected VAR representation returned for the VECM.")
  }
  
  # Extract deterministic component from the converted object.
  det_mat <- as.matrix(vec2var_model$deterministic)
  
  intercept <- if (
    !is.null(colnames(det_mat)) &&
    "constant" %in% colnames(det_mat)
  ) {
    as.numeric(det_mat[, "constant"])
  } else if (
    !is.null(colnames(det_mat)) &&
    "Intercept" %in% colnames(det_mat)
  ) {
    as.numeric(det_mat[, "Intercept"])
  } else {
    rep(0, k)
  }
  
  E <- as.matrix(stats::residuals(model))
  
  if (ncol(E) != k && nrow(E) == k) {
    E <- t(E)
  }
  
  if (ncol(E) != k) {
    stop("VECM residual matrix has unexpected dimensions.")
  }
  
  Sigma <- crossprod(E) / nrow(E)
  
  alpha <- as.matrix(tsDyn::coefA(model))
  
  beta_full <- as.matrix(
    model$model.specific$beta
  )
  
  beta <- beta_full[
    seq_len(k),
    ,
    drop = FALSE
  ]
  
  beta_const <- beta_full[
    k + 1L,
    ,
    drop = FALSE
  ]
  
  list(
    A_list = A_list,
    A1 = A_list[[1L]],
    Sigma = Sigma,
    residuals = E,
    intercept = intercept,
    p = p,
    alpha = alpha,
    beta = beta,
    beta_const = beta_const,
    vec2var_model = vec2var_model
  )
}


# --------------------------- 4. IRF CALCULATIONS --------------

# Moving-average representation:
#
# Psi_0 = I
# Psi_h = A_1 Psi_(h-1) + ... + A_p Psi_(h-p)
#
# OIRF uses Cholesky decomposition and therefore depends
# on the ordering of variables.
#
# GIRF does not depend on variable ordering.

calc_irfs <- function(
    A_list,
    Sigma,
    variable_names,
    horizon = HORIZON
) {
  
  k <- length(variable_names)
  p <- length(A_list)
  
  if (p < 1L) {
    stop("At least one autoregressive matrix is required.")
  }
  
  Sigma <- (Sigma + t(Sigma)) / 2
  
  eigenvalues <- eigen(
    Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  if (min(eigenvalues) <= 0) {
    stop("Residual covariance matrix is not positive definite.")
  }
  
  chol_lower <- t(chol(Sigma))
  
  out_oirf <- array(
    NA_real_,
    dim = c(horizon + 1L, k, k),
    dimnames = list(
      horizon = as.character(0:horizon),
      response = variable_names,
      impulse = variable_names
    )
  )
  
  out_girf <- out_oirf
  
  Psi <- vector("list", horizon + 1L)
  Psi[[1L]] <- diag(k)
  
  for (h in 0:horizon) {
    
    if (h > 0L) {
      
      Psi_h <- matrix(0, k, k)
      
      for (i in seq_len(min(p, h))) {
        
        Psi_h <- Psi_h +
          A_list[[i]] %*%
          Psi[[h - i + 1L]]
      }
      
      Psi[[h + 1L]] <- Psi_h
    }
    
    for (j in seq_len(k)) {
      
      out_oirf[h + 1L, , j] <-
        as.vector(
          Psi[[h + 1L]] %*% chol_lower[, j]
        )
      
      out_girf[h + 1L, , j] <-
        as.vector(
          Psi[[h + 1L]] %*%
            Sigma[, j] /
            sqrt(Sigma[j, j])
        )
    }
  }
  
  list(
    oirf = out_oirf,
    girf = out_girf
  )
}


# --------------------------- 4A. ORTHOGONALIZED FEVD ----------

calc_ofevd <- function(
    model,
    model_type,
    parameters,
    variable_names,
    horizon = HORIZON
) {
  
  if (model_type == "VAR") {
    
    fevd_model <- model
    
  } else if (model_type == "VECM") {
    
    fevd_model <- parameters$vec2var_model
    
  } else {
    
    stop("Unknown model type: ", model_type)
  }
  
  fevd_result <- vars::fevd(
    fevd_model,
    n.ahead = horizon
  )
  
  k <- length(variable_names)
  
  out <- array(
    NA_real_,
    dim = c(horizon, k, k),
    dimnames = list(
      horizon = as.character(seq_len(horizon)),
      response = variable_names,
      shock = variable_names
    )
  )
  
  for (i in seq_len(k)) {
    
    fevd_matrix <- as.matrix(
      fevd_result[[i]]
    )
    
    if (!all(dim(fevd_matrix) == c(horizon, k))) {
      stop(
        "Unexpected FEVD dimensions for variable: ",
        variable_names[i]
      )
    }
    
    out[, i, ] <- fevd_matrix * 100
  }
  
  row_sums <- apply(out, c(1, 2), sum)
  
  if (
    any(!is.finite(out)) ||
    any(abs(row_sums - 100) > 1e-6)
  ) {
    stop("OFEVD validation failed.")
  }
  
  out
}


# --------------------------- 4B. GENERALIZED FEVD --------------

calc_gfevd <- function(
    A_list,
    Sigma,
    variable_names,
    horizon = HORIZON
) {
  
  k <- length(variable_names)
  p <- length(A_list)
  
  if (any(vapply(
    A_list,
    function(A) !all(dim(A) == c(k, k)),
    logical(1)
  ))) {
    stop("Unexpected dimensions in autoregressive matrices.")
  }
  
  if (
    !all(dim(Sigma) == c(k, k)) ||
    any(!is.finite(Sigma))
  ) {
    stop("Invalid residual covariance matrix.")
  }
  
  Sigma <- (Sigma + t(Sigma)) / 2
  
  eigenvalues <- eigen(
    Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  if (
    min(eigenvalues) <= 0 ||
    any(diag(Sigma) <= 0)
  ) {
    stop("Residual covariance matrix must be positive definite.")
  }
  
  # Psi_0 = I.
  # Psi_h = A_1 Psi_(h-1) + ... + A_p Psi_(h-p).
  
  Psi <- vector("list", horizon)
  Psi[[1L]] <- diag(k)
  
  if (horizon >= 2L) {
    
    for (h in 2:horizon) {
      
      Psi_h <- matrix(0, k, k)
      
      for (i in seq_len(min(p, h - 1L))) {
        
        Psi_h <- Psi_h +
          A_list[[i]] %*%
          Psi[[h - i]]
      }
      
      Psi[[h]] <- Psi_h
    }
  }
  
  cumulative_girf <- matrix(0, k, k)
  
  out <- array(
    NA_real_,
    dim = c(horizon, k, k),
    dimnames = list(
      horizon = as.character(seq_len(horizon)),
      response = variable_names,
      shock = variable_names
    )
  )
  
  for (h in seq_len(horizon)) {
    
    # Generalized impact:
    # Psi_h Sigma[,j] / sqrt(Sigma[j,j])
    
    generalized_impact <- sweep(
      Psi[[h]] %*% Sigma,
      MARGIN = 2,
      STATS = sqrt(diag(Sigma)),
      FUN = "/"
    )
    
    cumulative_girf <- cumulative_girf +
      generalized_impact^2
    
    denominator <- rowSums(cumulative_girf)
    
    if (
      any(!is.finite(denominator)) ||
      any(denominator <= 0)
    ) {
      stop("Invalid denominator in generalized FEVD.")
    }
    
    out[h, , ] <- sweep(
      cumulative_girf,
      MARGIN = 1,
      STATS = denominator,
      FUN = "/"
    ) * 100
  }
  
  row_sums <- apply(out, c(1, 2), sum)
  
  if (
    any(!is.finite(out)) ||
    any(abs(row_sums - 100) > 1e-6)
  ) {
    stop("GFEVD validation failed.")
  }
  
  out
}


# --------------------------- 4C. GFEVD BLOCK BOOTSTRAP ---------

# Circular moving block bootstrap of multivariate residuals.
#
# Entire residual vectors are resampled together to preserve
# contemporaneous cross-equation dependence.
#
# Consecutive observations within each block retain their
# temporal ordering and local dependence structure.

sample_residual_blocks <- function(
    residuals,
    sample_size,
    block_length
) {
  
  E <- as.matrix(residuals)
  
  if (
    nrow(E) < 2L ||
    ncol(E) < 1L ||
    any(!is.finite(E))
  ) {
    stop("Invalid residual matrix.")
  }
  
  if (
    sample_size < 1L ||
    block_length < 1L ||
    block_length > nrow(E)
  ) {
    stop("Invalid bootstrap block dimensions.")
  }
  
  n_residuals <- nrow(E)
  
  # Center residuals equation by equation.
  E <- sweep(
    E,
    2L,
    colMeans(E),
    FUN = "-"
  )
  
  # Number of blocks needed to construct the sample.
  n_blocks <- ceiling(
    sample_size / block_length
  )
  
  # Random starting positions for circular blocks.
  starts <- sample.int(
    n_residuals,
    size = n_blocks,
    replace = TRUE
  )
  
  # Generate circular block indices.
  block_indices <- unlist(
    lapply(
      starts,
      function(start) {
        (
          (start - 1L + seq_len(block_length) - 1L) %%
            n_residuals
        ) + 1L
      }
    ),
    use.names = FALSE
  )
  
  # Retain exactly the required number of observations.
  selected_indices <- block_indices[
    seq_len(sample_size)
  ]
  
  E[selected_indices, , drop = FALSE]
}


# Residual-based moving block bootstrap for GFEVD.
#
# The same model specification is retained in every replication:
# - VAR(1) for stationary log returns;
# - VECM with the original cointegration rank and deterministic
#   specification for cointegrated log price levels.
#
# Models are re-estimated in every bootstrap replication.

bootstrap_gfevd <- function(
    model_type,
    y,
    point_parameters,
    runs = GFEVD_BOOT_RUNS,
    confidence_level = GFEVD_BOOT_CI,
    block_length = GFEVD_BLOCK_LENGTH,
    seed = BOOT_SEED
) {
  
  set.seed(seed)
  
  y <- as.matrix(y)
  
  n <- nrow(y)
  k <- ncol(y)
  p <- point_parameters$p
  
  if (
    n <= p ||
    k < 2L
  ) {
    stop("Insufficient observations or variables.")
  }
  
  if (!model_type %in% c("VAR", "VECM")) {
    stop("Unknown model type: ", model_type)
  }
  
  if (
    runs < 1L ||
    confidence_level <= 0 ||
    confidence_level >= 1
  ) {
    stop("Invalid bootstrap settings.")
  }
  
  E <- as.matrix(
    point_parameters$residuals
  )
  
  # The recursive simulation requires n-p residual vectors.
  n_required <- n - p
  
  if (nrow(E) < n_required) {
    stop(
      "Insufficient residual observations: ",
      nrow(E),
      " available; ",
      n_required,
      " required."
    )
  }
  
  # Keep the same residual alignment as in the original
  # recursive bootstrap implementation.
  E <- E[
    seq_len(n_required),
    ,
    drop = FALSE
  ]
  
  if (ncol(E) != k || any(!is.finite(E))) {
    stop("Residual dimensions do not match the model.")
  }
  
  # Automatic initial block-length rule.
  if (is.null(block_length)) {
    
    block_length <- ceiling(
      n_required^(1 / 3)
    )
    
  }
  
  block_length <- as.integer(block_length)
  
  if (
    block_length < 1L ||
    block_length > n_required
  ) {
    stop("Invalid block length.")
  }
  
  boot_gfevd <- vector(
    "list",
    runs
  )
  
  success <- 0L
  failed <- 0L
  
  failure_messages <- character(0)
  
  for (b in seq_len(runs)) {
    
    tryCatch({
      
      # 1. Generate a bootstrap sequence of residual vectors.
      E_star <- sample_residual_blocks(
        residuals = E,
        sample_size = n_required,
        block_length = block_length
      )
      
      # 2. Initialize the artificial time series.
      y_star <- matrix(
        NA_real_,
        nrow = n,
        ncol = k,
        dimnames = list(
          NULL,
          colnames(y)
        )
      )
      
      # Retain the initial p observations.
      y_star[seq_len(p), ] <- y[
        seq_len(p),
        ,
        drop = FALSE
      ]
      
      # 3. Recursively generate bootstrap observations
      #    using the original estimated VAR representation.
      for (tt in seq.int(p + 1L, n)) {
        
        pred <- point_parameters$intercept
        
        for (i in seq_len(p)) {
          
          pred <- pred +
            as.vector(
              y_star[
                tt - i,
                ,
                drop = FALSE
              ] %*%
                t(point_parameters$A_list[[i]])
            )
          
        }
        
        residual_index <- tt - p
        
        y_star[tt, ] <- pred +
          E_star[residual_index, ]
        
      }
      
      if (any(!is.finite(y_star))) {
        stop("Non-finite bootstrap observations.")
      }
      
      # 4. Re-estimate the original model specification.
      if (model_type == "VAR") {
        
        boot_fit <- fit_var1(y_star)
        
        boot_parameters <- extract_var_parameters(
          boot_fit,
          y_star
        )
        
      } else {
        
        boot_fit <- fit_vecm(y_star)
        
        boot_parameters <- extract_vecm_parameters(
          boot_fit
        )
        
      }
      
      # 5. Recalculate GFEVD using bootstrap estimates.
      gfevd_b <- calc_gfevd(
        A_list = boot_parameters$A_list,
        Sigma = boot_parameters$Sigma,
        variable_names = colnames(y),
        horizon = HORIZON
      )
      
      if (any(!is.finite(gfevd_b))) {
        stop("Non-finite GFEVD estimates.")
      }
      
      success <- success + 1L
      
      boot_gfevd[[success]] <- gfevd_b
      
    }, error = function(e) {
      
      failed <<- failed + 1L
      
      # Keep a small sample of error messages for diagnostics.
      if (length(failure_messages) < 10L) {
        
        failure_messages <<- c(
          failure_messages,
          conditionMessage(e)
        )
        
      }
      
    })
    
  }
  
  # 6. Validate the number of successful replications.
  if (success == 0L) {
    stop("All GFEVD bootstrap replications failed.")
  }
  
  if (success < max(50L, floor(runs * 0.5))) {
    
    stop(
      "Too few successful GFEVD bootstrap replications: ",
      success,
      " / ",
      runs,
      ". Check model convergence and residual dimensions."
    )
    
  }
  
  if (success < runs) {
    
    warning(
      "Some GFEVD bootstrap replications failed: ",
      failed,
      " / ",
      runs
    )
    
  }
  
  boot_gfevd <- boot_gfevd[
    seq_len(success)
  ]
  
  # 7. Combine bootstrap estimates into a four-dimensional array.
  #
  # Dimensions:
  # 1 - forecast horizon
  # 2 - response variable
  # 3 - shock variable
  # 4 - bootstrap replication
  
  first_result <- boot_gfevd[[1L]]
  
  boot_array <- array(
    unlist(
      boot_gfevd,
      use.names = FALSE
    ),
    dim = c(
      dim(first_result),
      success
    ),
    dimnames = c(
      dimnames(first_result),
      list(
        replication = as.character(
          seq_len(success)
        )
      )
    )
  )
  
  # 8. Calculate percentile confidence intervals.
  alpha_ci <- (1 - confidence_level) / 2
  
  lower <- apply(
    boot_array,
    MARGIN = 1:3,
    FUN = quantile,
    probs = alpha_ci,
    na.rm = TRUE,
    names = FALSE
  )
  
  upper <- apply(
    boot_array,
    MARGIN = 1:3,
    FUN = quantile,
    probs = 1 - alpha_ci,
    na.rm = TRUE,
    names = FALSE
  )
  
  list(
    lower = lower,
    upper = upper,
    successful_runs = success,
    failed_runs = failed,
    block_length = block_length,
    confidence_level = confidence_level,
    failure_messages = failure_messages
  )
  
}


# --------------------------- 5. WILD RESIDUAL BOOTSTRAP -------

# A common Rademacher multiplier is applied to the residual
# vector at each time point. This retains contemporaneous
# cross-equation dependence.
#
# Bootstrap confidence intervals are calculated for IRFs.
# FEVD confidence intervals are not calculated in this script.

bootstrap_model <- function(
    model_type,
    y,
    fitted_model,
    point_parameters,
    runs = BOOT_RUNS,
    seed = BOOT_SEED
) {
  
  set.seed(seed)
  
  n <- nrow(y)
  k <- ncol(y)
  p <- point_parameters$p
  
  boot_oirf <- vector("list", runs)
  boot_girf <- vector("list", runs)
  
  success <- 0L
  
  E <- point_parameters$residuals
  
  # Center residuals before applying wild multipliers.
  E <- sweep(
    E,
    2L,
    colMeans(E),
    FUN = "-"
  )
  
  for (b in seq_len(runs)) {
    
    tryCatch({
      
      multiplier <- sample(
        c(-1, 1),
        nrow(E),
        replace = TRUE
      )
      
      E_star <- sweep(
        E,
        1L,
        multiplier,
        FUN = "*"
      )
      
      y_star <- matrix(
        NA_real_,
        nrow = n,
        ncol = k,
        dimnames = list(NULL, colnames(y))
      )
      
      # Retain the initial p observations.
      y_star[seq_len(p), ] <- y[
        seq_len(p),
        ,
        drop = FALSE
      ]
      
      # Recursive simulation from the estimated level VAR.
      for (tt in (p + 1L):n) {
        
        pred <- point_parameters$intercept
        
        for (i in seq_len(p)) {
          
          pred <- pred +
            as.vector(
              y_star[tt - i, , drop = FALSE] %*%
                t(point_parameters$A_list[[i]])
            )
        }
        
        residual_index <- tt - p
        
        if (residual_index > nrow(E_star)) {
          stop("Bootstrap residual index exceeds available residuals.")
        }
        
        y_star[tt, ] <- pred + E_star[residual_index, ]
      }
      
      if (model_type == "VAR") {
        
        boot_fit <- fit_var1(y_star)
        
        boot_par <- extract_var_parameters(
          boot_fit,
          y_star
        )
        
      } else {
        
        boot_fit <- fit_vecm(y_star)
        
        boot_par <- extract_vecm_parameters(
          boot_fit
        )
      }
      
      irf_b <- calc_irfs(
        boot_par$A_list,
        boot_par$Sigma,
        colnames(y),
        HORIZON
      )
      
      success <- success + 1L
      
      boot_oirf[[success]] <- irf_b$oirf
      boot_girf[[success]] <- irf_b$girf
      
    }, error = function(e) {
      NULL
    })
  }
  
  if (success == 0L) {
    stop("All bootstrap replications failed.")
  }
  
  if (success < max(50L, floor(runs * 0.5))) {
    warning(
      "Too few successful bootstrap replications: ",
      success, " / ", runs
    )
  }
  
  boot_oirf <- boot_oirf[seq_len(success)]
  boot_girf <- boot_girf[seq_len(success)]
  
  arr_oirf <- array(
    unlist(boot_oirf),
    dim = c(
      dim(boot_oirf[[1L]]),
      success
    )
  )
  
  arr_girf <- array(
    unlist(boot_girf),
    dim = c(
      dim(boot_girf[[1L]]),
      success
    )
  )
  
  alpha_ci <- (1 - BOOT_CI) / 2
  
  list(
    oirf_lower = apply(
      arr_oirf,
      1:3,
      quantile,
      probs = alpha_ci,
      na.rm = TRUE
    ),
    
    oirf_upper = apply(
      arr_oirf,
      1:3,
      quantile,
      probs = 1 - alpha_ci,
      na.rm = TRUE
    ),
    
    girf_lower = apply(
      arr_girf,
      1:3,
      quantile,
      probs = alpha_ci,
      na.rm = TRUE
    ),
    
    girf_upper = apply(
      arr_girf,
      1:3,
      quantile,
      probs = 1 - alpha_ci,
      na.rm = TRUE
    ),
    
    successful_runs = success
  )
}


# --------------------------- 6. DIAGNOSTICS -------------------

portmanteau_test <- function(
    model,
    model_type,
    lags = PORTMANTEAU_LAGS
) {
  
  if (model_type == "VECM") {
    
    test_object <- tsDyn:::vec2var.tsDyn(model)
    
  } else {
    
    test_object <- model
  }
  
  effective_lags <- min(
    as.integer(lags),
    floor(test_object$obs / 5)
  )
  
  empty_result <- c(
    statistic = NA_real_,
    df = NA_real_,
    p_value = NA_real_
  )
  
  if (effective_lags <= test_object$p) {
    
    return(list(
      asymptotic = empty_result,
      adjusted = empty_result,
      lags = effective_lags,
      level_VAR_order = test_object$p
    ))
  }
  
  find_htest <- function(x) {
    
    if (inherits(x, "htest")) {
      return(x)
    }
    
    if (is.list(x)) {
      
      for (item in x) {
        
        found <- find_htest(item)
        
        if (!is.null(found)) {
          return(found)
        }
      }
    }
    
    NULL
  }
  
  extract_htest <- function(x, test_label) {
    
    htest <- find_htest(x)
    
    if (is.null(htest)) {
      stop(
        "Could not locate htest result for ",
        test_label,
        ". Inspect vars::serial.test() output."
      )
    }
    
    c(
      statistic = as.numeric(htest$statistic)[1L],
      df = as.numeric(htest$parameter)[1L],
      p_value = as.numeric(htest$p.value)[1L]
    )
  }
  
  asymptotic_result <- vars::serial.test(
    test_object,
    lags.pt = effective_lags,
    type = "PT.asymptotic"
  )
  
  adjusted_result <- vars::serial.test(
    test_object,
    lags.pt = effective_lags,
    type = "PT.adjusted"
  )
  
  list(
    asymptotic = extract_htest(
      asymptotic_result,
      "PT.asymptotic"
    ),
    
    adjusted = extract_htest(
      adjusted_result,
      "PT.adjusted"
    ),
    
    lags = effective_lags,
    level_VAR_order = test_object$p
  )
}


# Univariate ARCH-LM test.
arch_lm_test <- function(
    x,
    lags = ARCH_LAGS
) {
  
  x <- as.numeric(x)
  x <- x - mean(x)
  
  n <- length(x)
  
  lags <- min(
    lags,
    floor(n / 10)
  )
  
  if (lags < 1L) {
    return(c(
      statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_
    ))
  }
  
  z <- x^2
  
  yy <- z[(lags + 1L):n]
  
  xx <- sapply(
    seq_len(lags),
    function(j) {
      z[(lags + 1L - j):(n - j)]
    }
  )
  
  fit <- lm.fit(
    cbind(1, xx),
    yy
  )
  
  tss <- sum(
    (yy - mean(yy))^2
  )
  
  r2 <- if (tss == 0) {
    NA_real_
  } else {
    1 - sum(fit$residuals^2) / tss
  }
  
  LM <- length(yy) * r2
  
  c(
    statistic = LM,
    df = lags,
    p_value = pchisq(
      LM,
      df = lags,
      lower.tail = FALSE
    )
  )
}


# Jarque-Bera test.
# Central moments are calculated consistently using denominator n.
jarque_bera_test <- function(x) {
  
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  n <- length(x)
  
  if (n < 3L) {
    return(c(
      statistic = NA_real_,
      df = 2,
      p_value = NA_real_
    ))
  }
  
  centered <- x - mean(x)
  
  m2 <- mean(centered^2)
  m3 <- mean(centered^3)
  m4 <- mean(centered^4)
  
  if (!is.finite(m2) || m2 <= 0) {
    return(c(
      statistic = NA_real_,
      df = 2,
      p_value = NA_real_
    ))
  }
  
  skewness <- m3 / m2^(3 / 2)
  kurtosis <- m4 / m2^2
  
  JB <- n / 6 * (
    skewness^2 +
      (kurtosis - 3)^2 / 4
  )
  
  c(
    statistic = JB,
    df = 2,
    p_value = pchisq(
      JB,
      df = 2,
      lower.tail = FALSE
    )
  )
}


# Companion matrix roots.
companion_roots <- function(A_list) {
  
  p <- length(A_list)
  k <- nrow(A_list[[1L]])
  
  # VAR(1): roots are eigenvalues of A1.
  if (p == 1L) {
    return(eigen(A_list[[1L]])$values)
  }
  
  # Companion matrix for VAR(p).
  companion <- matrix(
    0,
    nrow = k * p,
    ncol = k * p
  )
  
  # First block row: [A1 A2 ... Ap]
  companion[
    seq_len(k),
    seq_len(k * p)
  ] <- do.call(cbind, A_list)
  
  # Lower block rows: identity matrices.
  companion[
    (k + 1L):(k * p),
    seq_len(k * (p - 1L))
  ] <- diag(k * (p - 1L))
  
  eigen(companion)$values
}


# Stability assessment.
#
# For a cointegrated system with rank r, k-r unit roots
# are expected because the system contains k-r common
# stochastic trends. These roots should not be treated
# as evidence of instability.
#
# The transitory roots are checked separately.

stability_summary <- function(
    model_type,
    parameters,
    k
) {
  
  eig <- companion_roots(
    parameters$A_list
  )
  
  mod <- Mod(eig)
  
  if (model_type == "VAR") {
    
    return(list(
      roots = eig,
      status = if (all(mod < 1)) {
        "stable"
      } else {
        "unstable"
      },
      max_modulus = max(mod),
      expected_unit_roots = 0L,
      transitory_roots = eig
    ))
  }
  
  n_unit_expected <- k - VECM_RANK
  
  # Identify roots closest to +1.
  near_one <- order(
    abs(eig - 1)
  )[seq_len(n_unit_expected)]
  
  transitory <- eig[-near_one]
  
  explosive <- any(
    Mod(transitory) > 1 + 1e-6
  )
  
  list(
    roots = eig,
    status = if (explosive) {
      "explosive transitory root"
    } else {
      "no explosive transitory root"
    },
    max_modulus = max(mod),
    expected_unit_roots = n_unit_expected,
    transitory_roots = transitory
  )
}


# --------------------------- 7. OUTPUT HELPERS ----------------

array_to_long <- function(
    arr,
    model_id,
    irf_type,
    lower = NULL,
    upper = NULL
) {
  
  dn <- dimnames(arr)
  
  grid <- expand.grid(
    horizon = as.integer(dn[[1L]]),
    response = dn[[2L]],
    impulse = dn[[3L]],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  grid$value <- as.vector(arr)
  grid$model <- model_id
  grid$irf_type <- irf_type
  
  if (!is.null(lower)) {
    grid$lower <- as.vector(lower)
  }
  
  if (!is.null(upper)) {
    grid$upper <- as.vector(upper)
  }
  
  grid
}


save_irf_plot <- function(
    df,
    model_id,
    irf_type
) {
  
  p <- ggplot(
    df,
    aes(x = horizon, y = value)
  ) +
    
    geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      colour = "grey55"
    ) +
    
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.18,
      na.rm = TRUE
    ) +
    
    geom_line(
      linewidth = 0.45
    ) +
    
    facet_grid(
      response ~ impulse,
      scales = "free_y"
    ) +
    
    labs(
      title = paste(
        irf_type,
        "impulse responses:",
        model_id
      ),
      subtitle = paste(
        "Rows: response variables;",
        "columns: impulse variables"
      ),
      x = "Horizon (trading days)",
      y = "Response"
    ) +
    
    theme_minimal(base_size = 8) +
    
    theme(
      strip.text = element_text(size = 6),
      plot.title = element_text(size = 11),
      panel.spacing = grid::unit(0.45, "lines")
    )
  
  ggsave(
    file.path(
      PLOTS_DIR,
      paste0(
        model_id,
        "_",
        tolower(irf_type),
        ".png"
      )
    ),
    p,
    width = 14,
    height = max(
      8,
      1.65 * length(unique(df$response))
    ),
    dpi = 300,
    limitsize = FALSE
  )
}


fevd_array_to_long <- function(
    arr,
    model_id,
    fevd_type
) {
  
  df <- array_to_long(
    arr,
    model_id,
    fevd_type
  )
  
  names(df)[names(df) == "impulse"] <- "shock"
  names(df)[names(df) == "irf_type"] <- "fevd_type"
  
  df
}


save_fevd_plot <- function(
    df,
    model_id,
    fevd_type
) {
  
  p <- ggplot(
    df,
    aes(x = horizon, y = value)
  ) +
    
    geom_line(
      linewidth = 0.5
    ) +
    
    facet_grid(
      response ~ shock,
      scales = "fixed"
    ) +
    
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20)
    ) +
    
    labs(
      title = paste(
        fevd_type,
        ":",
        model_id
      ),
      subtitle = paste(
        "Rows: variables whose forecast error variance is decomposed;",
        "columns: shock sources"
      ),
      x = "Forecast horizon (trading days)",
      y = "Contribution, %"
    ) +
    
    theme_minimal(base_size = 8) +
    
    theme(
      strip.text = element_text(size = 6),
      plot.title = element_text(size = 11),
      panel.spacing = grid::unit(0.45, "lines")
    )
  
  ggsave(
    file.path(
      PLOTS_DIR,
      paste0(
        model_id,
        "_",
        tolower(fevd_type),
        ".png"
      )
    ),
    p,
    width = 14,
    height = max(
      8,
      1.65 * length(unique(df$response))
    ),
    dpi = 300,
    limitsize = FALSE
  )
}


# Plot GFEVD estimates with bootstrap confidence intervals.

save_gfevd_bootstrap_plot <- function(
    df,
    model_id,
    confidence_level = GFEVD_BOOT_CI
) {
  
  p <- ggplot(
    df,
    aes(
      x = horizon,
      y = value
    )
  ) +
    
    # Bootstrap confidence interval.
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper
      ),
      fill = "#4C78A8",
      alpha = 0.20,
      na.rm = TRUE
    ) +
    
    # Point estimate.
    geom_line(
      colour = "#1F4E79",
      linewidth = 0.55
    ) +
    
    facet_grid(
      response ~ shock,
      scales = "fixed"
    ) +
    
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20)
    ) +
    
    labs(
      title = paste(
        "GFEVD with bootstrap confidence intervals:",
        model_id
      ),
      
      subtitle = paste0(
        round(confidence_level * 100),
        "% percentile confidence intervals; ",
        "moving block residual bootstrap"
      ),
      
      x = "Forecast horizon (trading days)",
      y = "Contribution, %"
    ) +
    
    theme_minimal(
      base_size = 8
    ) +
    
    theme(
      strip.text = element_text(size = 6),
      plot.title = element_text(size = 11),
      panel.spacing = grid::unit(
        0.45,
        "lines"
      )
    )
  
  ggsave(
    filename = file.path(
      PLOTS_DIR,
      paste0(
        model_id,
        "_gfevd_bootstrap.png"
      )
    ),
    
    plot = p,
    
    width = 14,
    
    height = max(
      8,
      1.65 * length(unique(df$response))
    ),
    
    dpi = 300,
    
    limitsize = FALSE
  )
  
}

# --------------------------- 8. MAIN LOOP ---------------------

set.seed(BOOT_SEED)

all_models <- list()
model_rows <- list()
diag_rows <- list()
root_rows <- list()
irf_rows <- list()
fevd_rows <- list()
gfevd_bootstrap_rows <- list()

for (dataset in ALL_DATASETS) {
  
  model_id <- dataset
  
  # Map period datasets to their corresponding levels datasets.
  levels_id <- sub(
    "^period",
    "levels",
    dataset
  )
  
  is_vecm <- levels_id %in% VECM_DATASETS
  
  data_file_id <- if (is_vecm) {
    levels_id
  } else {
    dataset
  }
  
  message(
    "Estimating ",
    model_id,
    " [",
    if (is_vecm) {
      paste0(
        "VECM(r=", VECM_RANK,
        ", lag=", VECM_LAGGED_DIFFERENCES,
        "; VAR(", VECM_VAR_ORDER, ") in levels)"
      )
    } else {
      "VAR(1)"
    },
    "]"
  )
  
  y <- read_model_data(data_file_id)
  k <- ncol(y)
  
  if (is_vecm) {
    
    fitted_model <- fit_vecm(y)
    
    model_type <- "VECM"
    
    parameters <- extract_vecm_parameters(
      fitted_model
    )
    
  } else {
    
    fitted_model <- fit_var1(y)
    
    model_type <- "VAR"
    
    parameters <- extract_var_parameters(
      fitted_model,
      y
    )
  }
  
  all_models[[model_id]] <- fitted_model
  
  saveRDS(
    fitted_model,
    file.path(
      MODELS_DIR,
      paste0(model_id, ".rds")
    )
  )
  
  capture.output(
    summary(fitted_model),
    file = file.path(
      TABLES_DIR,
      paste0(model_id, "_summary.txt")
    )
  )
  
  
  # ---------------- MODEL INVENTORY ----------------
  
  model_rows[[model_id]] <- data.frame(
    model = model_id,
    data_used = data_file_id,
    type = model_type,
    n_obs = nrow(y),
    n_variables = k,
    variables = paste(
      colnames(y),
      collapse = " | "
    ),
    rank = if (is_vecm) {
      VECM_RANK
    } else {
      NA_integer_
    },
    lagged_differences = if (is_vecm) {
      VECM_LAGGED_DIFFERENCES
    } else {
      NA_integer_
    },
    vecm_level_VAR_order = if (is_vecm) {
      VECM_VAR_ORDER
    } else {
      NA_integer_
    },
    VAR_lag_order = if (is_vecm) {
      VECM_VAR_ORDER
    } else {
      1L
    },
    deterministic_specification = if (is_vecm) {
      "Case II: restricted constant"
    } else {
      "Unrestricted intercept"
    },
    stringsAsFactors = FALSE
  )
  
  
  # ---------------- DIAGNOSTICS ----------------
  
  E <- parameters$residuals
  
  pm <- portmanteau_test(
    fitted_model,
    model_type,
    lags = PORTMANTEAU_LAGS
  )
  
  diag_rows[[
    paste0(model_id, "_portmanteau_asymptotic")
  ]] <- data.frame(
    model = model_id,
    test = paste0(
      "Portmanteau asymptotic (h=",
      pm$lags,
      ", level VAR p=",
      pm$level_VAR_order,
      ")"
    ),
    variable = NA_character_,
    statistic = pm$asymptotic["statistic"],
    df = pm$asymptotic["df"],
    p_value = pm$asymptotic["p_value"],
    stringsAsFactors = FALSE
  )
  
  diag_rows[[
    paste0(model_id, "_portmanteau_adjusted")
  ]] <- data.frame(
    model = model_id,
    test = paste0(
      "Portmanteau adjusted (h=",
      pm$lags,
      ", level VAR p=",
      pm$level_VAR_order,
      ")"
    ),
    variable = NA_character_,
    statistic = pm$adjusted["statistic"],
    df = pm$adjusted["df"],
    p_value = pm$adjusted["p_value"],
    stringsAsFactors = FALSE
  )
  
  for (j in seq_len(k)) {
    
    a <- arch_lm_test(E[, j])
    jb <- jarque_bera_test(E[, j])
    
    diag_rows[[
      paste0(model_id, "_arch_", j)
    ]] <- data.frame(
      model = model_id,
      test = paste0("ARCH-LM(", ARCH_LAGS, ")"),
      variable = colnames(y)[j],
      statistic = a["statistic"],
      df = a["df"],
      p_value = a["p_value"],
      stringsAsFactors = FALSE
    )
    
    diag_rows[[
      paste0(model_id, "_jb_", j)
    ]] <- data.frame(
      model = model_id,
      test = "Jarque-Bera (univariate residual)",
      variable = colnames(y)[j],
      statistic = jb["statistic"],
      df = jb["df"],
      p_value = jb["p_value"],
      stringsAsFactors = FALSE
    )
  }
  
  
  # ---------------- STABILITY ----------------
  
  stab <- stability_summary(
    model_type,
    parameters,
    k
  )
  
  root_rows[[model_id]] <- data.frame(
    model = model_id,
    root = as.character(stab$roots),
    modulus = as.numeric(Mod(stab$roots)),
    expected_unit_roots = stab$expected_unit_roots,
    stability_status = stab$status,
    stringsAsFactors = FALSE
  )
  
  
  # ---------------- IRF AND BOOTSTRAP ----------------
  
  point_irfs <- calc_irfs(
    parameters$A_list,
    parameters$Sigma,
    colnames(y),
    HORIZON
  )
  
  boot <- bootstrap_model(
    model_type,
    y,
    fitted_model,
    parameters,
    runs = BOOT_RUNS,
    seed = BOOT_SEED + match(model_id, ALL_DATASETS)
  )
  
  oirf_df <- array_to_long(
    point_irfs$oirf,
    model_id,
    "OIRF",
    boot$oirf_lower,
    boot$oirf_upper
  )
  
  girf_df <- array_to_long(
    point_irfs$girf,
    model_id,
    "GIRF",
    boot$girf_lower,
    boot$girf_upper
  )
  
  irf_model_df <- rbind(
    oirf_df,
    girf_df
  )
  
  irf_rows[[model_id]] <- irf_model_df
  
  save_irf_plot(
    oirf_df,
    model_id,
    "OIRF"
  )
  
  save_irf_plot(
    girf_df,
    model_id,
    "GIRF"
  )
  
  write.csv(
    irf_model_df,
    file.path(
      TABLES_DIR,
      paste0(model_id, "_irf_values.csv")
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  
  # ---------------- FEVD ----------------
  
  ofevd_arr <- calc_ofevd(
    model = fitted_model,
    model_type = model_type,
    parameters = parameters,
    variable_names = colnames(y),
    horizon = HORIZON
  )
  
  gfevd_arr <- calc_gfevd(
    A_list = parameters$A_list,
    Sigma = parameters$Sigma,
    variable_names = colnames(y),
    horizon = HORIZON
  )
  
  
  # ---------------- GFEVD BOOTSTRAP ----------------
  
  gfevd_boot <- bootstrap_gfevd(
    model_type = model_type,
    y = y,
    point_parameters = parameters,
    runs = GFEVD_BOOT_RUNS,
    confidence_level = GFEVD_BOOT_CI,
    block_length = GFEVD_BLOCK_LENGTH,
    seed = BOOT_SEED +
      10000L +
      match(model_id, ALL_DATASETS)
  )
  
  # Store bootstrap metadata.
  gfevd_bootstrap_rows[[model_id]] <- data.frame(
    model = model_id,
    requested_runs = GFEVD_BOOT_RUNS,
    successful_runs = gfevd_boot$successful_runs,
    failed_runs = gfevd_boot$failed_runs,
    block_length = gfevd_boot$block_length,
    confidence_level = gfevd_boot$confidence_level,
    stringsAsFactors = FALSE
  )
  
  # Convert point estimates into a long-format data frame.
  gfevd_df <- fevd_array_to_long(
    gfevd_arr,
    model_id,
    "GFEVD"
  )
  
  # Attach bootstrap confidence intervals.
  gfevd_df$lower <- as.vector(
    gfevd_boot$lower
  )
  
  gfevd_df$upper <- as.vector(
    gfevd_boot$upper
  )
  
  # Save GFEVD with confidence intervals.
  save_gfevd_bootstrap_plot(
    gfevd_df,
    model_id,
    GFEVD_BOOT_CI
  )
  
  # Keep the original OFEVD output structure.
  ofevd_df <- fevd_array_to_long(
    ofevd_arr,
    model_id,
    "OFEVD"
  )
  
  ofevd_df$lower <- NA_real_
  ofevd_df$upper <- NA_real_
  
  # Combine orthogonalized and generalized FEVD.
  fevd_model_df <- rbind(
    ofevd_df,
    gfevd_df
  )
  
  fevd_rows[[model_id]] <- fevd_model_df
  
  # Preserve the original point-estimate plots.
  save_fevd_plot(
    ofevd_df,
    model_id,
    "OFEVD"
  )
  
  save_fevd_plot(
    gfevd_df,
    model_id,
    "GFEVD"
  )
  
  # Save point estimates and bootstrap confidence intervals.
  write.csv(
    fevd_model_df,
    file.path(
      TABLES_DIR,
      paste0(model_id, "_fevd_values.csv")
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  message(
    "  IRF bootstrap successful: ",
    boot$successful_runs,
    " / ",
    BOOT_RUNS,
    "; GFEVD bootstrap successful: ",
    gfevd_boot$successful_runs,
    " / ",
    GFEVD_BOOT_RUNS,
    "; block length: ",
    gfevd_boot$block_length
  )


# --------------------------- 9. SAVE CONSOLIDATED RESULTS ----

model_table <- do.call(
  rbind,
  model_rows
)

diagnostic_table <- do.call(
  rbind,
  diag_rows
)

root_table <- do.call(
  rbind,
  root_rows
)

irf_table <- do.call(
  rbind,
  irf_rows
)

fevd_table <- do.call(
  rbind,
  fevd_rows
)


write.csv(
  model_table,
  file.path(TABLES_DIR, "model_inventory.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  diagnostic_table,
  file.path(TABLES_DIR, "diagnostics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  root_table,
  file.path(TABLES_DIR, "roots_and_stability.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  irf_table,
  file.path(TABLES_DIR, "all_irf_values.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  fevd_table,
  file.path(TABLES_DIR, "all_fevd_values.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

saveRDS(
  all_models,
  file.path(MODELS_DIR, "all_models.rds")
)

gfevd_bootstrap_table <- do.call(
  rbind,
  gfevd_bootstrap_rows
)

write.csv(
  gfevd_bootstrap_table,
  file.path(
    TABLES_DIR,
    "gfevd_bootstrap_summary.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# --------------------------- 10. REPRODUCIBILITY --------------

writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt"),
  useBytes = TRUE
)


# Save variable ordering used for Cholesky OIRFs.
# The order is important because OIRFs and OFEVD depend on it.

ordering_rows <- lapply(
  ALL_DATASETS,
  function(dataset) {
    
    levels_id <- sub(
      "^period",
      "levels",
      dataset
    )
    
    is_vecm <- levels_id %in% VECM_DATASETS
    
    data_file_id <- if (is_vecm) {
      levels_id
    } else {
      dataset
    }
    
    x <- read_model_data(data_file_id)
    
    data.frame(
      model = dataset,
      position = seq_along(colnames(x)),
      variable = colnames(x),
      stringsAsFactors = FALSE
    )
  }
)

write.csv(
  do.call(rbind, ordering_rows),
  file.path(
    TABLES_DIR,
    "variable_ordering_for_OIRF.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# Save model specification and diagnostic notes.
writeLines(
  c(
    "MODEL SPECIFICATION",
    "",
    "VECM models:",
    paste(VECM_DATASETS, collapse = ", "),
    paste0("Cointegration rank r = ", VECM_RANK),
    paste0(
      "Lagged differences = ",
      VECM_LAGGED_DIFFERENCES
    ),
    paste0(
      "Equivalent VAR order in levels = ",
      VECM_VAR_ORDER
    ),
    "Deterministic specification: Case II, restricted constant.",
    "",
    "VAR models:",
    paste(VAR_DATASETS, collapse = ", "),
    "VAR order = 1.",
    "Deterministic specification: unrestricted intercept.",
    "",
    "IRF:",
    "OIRF uses Cholesky decomposition and depends on variable ordering.",
    "GIRF does not depend on variable ordering.",
    "",
    "FEVD:",
    "OFEVD is based on orthogonalized shocks.",
    "GFEVD is based on generalized shocks and is normalized to 100 percent.",
    "",
    "Bootstrap:",
    paste0("Bootstrap replications requested: ", BOOT_RUNS),
    paste0("Confidence level: ", BOOT_CI),
    "Bootstrap:",
    paste0("IRF bootstrap replications: ", BOOT_RUNS),
    paste0("IRF confidence level: ", BOOT_CI),
    "IRF confidence intervals use a common-vector",
    "Rademacher wild residual bootstrap.",
    
    paste0(
      "GFEVD bootstrap replications: ",
      GFEVD_BOOT_RUNS
    ),
    paste0(
      "GFEVD confidence level: ",
      GFEVD_BOOT_CI
    ),
    "GFEVD confidence intervals use circular",
    "multivariate moving block residual bootstrap.",
    "Residual vectors are resampled jointly.",
    "Block length follows the n^(1/3) rule unless specified.",
    "Bootstrap intervals are percentile-based.",
    "",
    "Stability:",
    "For a cointegrated system with rank r, k-r unit roots are expected.",
    "These roots represent common stochastic trends and are not",
    "classified as explosive transitory roots."
  ),
  file.path(
    TABLES_DIR,
    "model_specification.txt"
  )
)


message(
  "All models completed. Results saved to: ",
  normalizePath(OUTPUT_DIR)
)