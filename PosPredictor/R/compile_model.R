#' Find a Usable Boost Include Path for rstan
#'
#' @return Character string: path to Boost headers, or empty string if not found.
#' @keywords internal
.find_boost_path <- function() {
  # 1. Try BH package first (includes headers on most platforms)
  bh_inc <- system.file("include", package = "BH")
  if (nchar(bh_inc) && file.exists(file.path(bh_inc, "boost"))) {
    return(bh_inc)
  }
  # 2. Try system Boost (common Linux / macOS paths)
  candidates <- c("/usr/include", "/usr/local/include",
                  "/opt/homebrew/include")
  for (p in candidates) {
    if (file.exists(file.path(p, "boost"))) return(p)
  }
  ""  # not found
}

#' Compile the Hierarchical Stan Model and Save as RDS
#'
#' Compiles the Stan model bundled with PosPredictor and saves the resulting
#' \code{stanmodel} object as an RDS file for instant reuse.  Priors are
#' passed as data, so changing prior parameters does NOT require recompilation.
#'
#' @param output_path Character. Path where the compiled model RDS will be
#'   saved.  Defaults to \code{~/.PosPredictor/compiled_model.rds}.
#' @param verbose Logical. Print compilation progress (default TRUE).
#' @return Invisibly returns the path to the saved RDS file.
#' @export
compile_stan_model <- function(output_path = NULL, verbose = TRUE) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required. Install it with:\n",
         "  install.packages('rstan')")
  }

  stan_file <- system.file("stan", "hierarchical_model.stan",
                            package = "PosPredictor")
  if (!nchar(stan_file)) {
    stop("Stan model file not found in PosPredictor installation.")
  }

  if (is.null(output_path)) {
    dir_path <- file.path(path.expand("~"), ".PosPredictor")
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
    output_path <- file.path(dir_path, "compiled_model.rds")
  }

  if (verbose) message("Compiling Stan model (this takes ~60-120 seconds the first time)...")

  boost_lib <- .find_boost_path()
  if (nchar(boost_lib)) {
    compiled <- rstan::stan_model(file = stan_file, boost_lib = boost_lib)
  } else {
    compiled <- rstan::stan_model(file = stan_file)
  }

  saveRDS(compiled, file = output_path)
  if (verbose) message("Compiled model saved to: ", output_path)
  invisible(output_path)
}

#' Load the Pre-compiled Stan Model
#'
#' Loads the pre-compiled Stan model RDS.  If the file does not exist,
#' calls \code{\link{compile_stan_model}} first.
#'
#' @param rds_path Character. Path to the compiled model RDS.
#'   Defaults to \code{~/.PosPredictor/compiled_model.rds}.
#' @param verbose Logical. Verbosity (default TRUE).
#' @return A \code{stanmodel} object.
#' @export
load_stan_model <- function(rds_path = NULL, verbose = TRUE) {
  if (is.null(rds_path)) {
    rds_path <- file.path(path.expand("~"), ".PosPredictor",
                          "compiled_model.rds")
  }
  if (!file.exists(rds_path)) {
    if (verbose) message("Compiled model not found. Compiling now...")
    compile_stan_model(output_path = rds_path, verbose = verbose)
  }
  readRDS(rds_path)
}

#' Run MCMC via Stan to Compute PoS
#'
#' Fits the three-level hierarchical model with Hamiltonian Monte Carlo and
#' returns posterior samples plus the estimated PoS.
#'
#' @param hist_data data.frame from \code{\link{simulate_historical_data}}.
#' @param current_y_os Numeric. Current trial observed log(HR) for OS.
#' @param current_y_pfs Numeric. Current trial observed log(HR) for PFS.
#' @param current_se_os Numeric. Current trial SE for OS.
#' @param current_se_pfs Numeric. Current trial SE for PFS.
#' @param current_within_corr Numeric. Within-trial correlation (default 0.65).
#' @param target_os Numeric. OS success threshold (default \code{log(0.74)}).
#' @param mu_os_prior_mean Numeric. Prior mean for population OS log(HR).
#' @param mu_os_prior_sd Numeric. Prior SD for population OS log(HR).
#' @param mu_pfs_prior_mean Numeric. Prior mean for population PFS log(HR).
#' @param mu_pfs_prior_sd Numeric. Prior SD for population PFS log(HR).
#' @param tau_os_prior_sd Numeric. Half-normal SD prior for tau_os.
#' @param tau_pfs_prior_sd Numeric. Half-normal SD prior for tau_pfs.
#' @param rho_z_prior_sd Numeric. SD for Fisher-z prior on rho.
#' @param iter Integer. Total MCMC iterations per chain (default 2000).
#' @param warmup Integer. Warmup iterations (default 1000).
#' @param chains Integer. Number of chains (default 4).
#' @param adapt_delta Numeric. Target acceptance rate (default 0.97).
#' @param rds_path Character or NULL. Path to pre-compiled model RDS.
#' @param seed Integer. Random seed (default 123).
#' @return A list with: fit (stanfit), pos (scalar PoS), summary (data.frame),
#'   rhat_ok (logical), ess_ok (logical).
#' @export
compute_pos_mcmc <- function(hist_data,
                              current_y_os,
                              current_y_pfs,
                              current_se_os,
                              current_se_pfs,
                              current_within_corr = 0.65,
                              target_os           = log(0.74),
                              mu_os_prior_mean    = -0.30,
                              mu_os_prior_sd      =  0.50,
                              mu_pfs_prior_mean   = -0.40,
                              mu_pfs_prior_sd     =  0.50,
                              tau_os_prior_sd     =  0.25,
                              tau_pfs_prior_sd    =  0.25,
                              rho_z_prior_sd      =  1.50,
                              iter                = 2000,
                              warmup              = 1000,
                              chains              = 4,
                              adapt_delta         = 0.97,
                              rds_path            = NULL,
                              seed                = 123) {

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required.")
  }

  K <- nrow(hist_data)

  # Build within-trial covariance list for historical trials
  W_hist <- lapply(seq_len(K), function(k) {
    matrix(
      c(hist_data$se_os[k]^2,
        hist_data$within_corr[k] * hist_data$se_os[k] * hist_data$se_pfs[k],
        hist_data$within_corr[k] * hist_data$se_os[k] * hist_data$se_pfs[k],
        hist_data$se_pfs[k]^2),
      nrow = 2
    )
  })

  W_current <- matrix(
    c(current_se_os^2,
      current_within_corr * current_se_os * current_se_pfs,
      current_within_corr * current_se_os * current_se_pfs,
      current_se_pfs^2),
    nrow = 2
  )

  # Observed y: Stan array[K] vector[2] expects K x 2 matrix (row = trial, col = endpoint)
  y_hist_mat <- matrix(
    c(hist_data$log_hr_os, hist_data$log_hr_pfs),
    nrow = K, ncol = 2
  )
  # Stan array[K] matrix[2,2] expects a 3D array with dims [K, 2, 2]
  W_hist_arr <- array(0, dim = c(K, 2, 2))
  for (k in seq_len(K)) {
    W_hist_arr[k, , ] <- W_hist[[k]]
  }

  stan_data <- list(
    K                 = K,
    y_hist            = y_hist_mat,    # K x 2
    W_hist            = W_hist_arr,    # K x 2 x 2
    y_current         = c(current_y_os, current_y_pfs),
    W_current         = W_current,
    target_os         = target_os,
    mu_os_prior_mean  = mu_os_prior_mean,
    mu_os_prior_sd    = mu_os_prior_sd,
    mu_pfs_prior_mean = mu_pfs_prior_mean,
    mu_pfs_prior_sd   = mu_pfs_prior_sd,
    tau_os_prior_sd   = tau_os_prior_sd,
    tau_pfs_prior_sd  = tau_pfs_prior_sd,
    rho_z_prior_sd    = rho_z_prior_sd
  )

  model <- load_stan_model(rds_path = rds_path, verbose = TRUE)

  fit <- rstan::sampling(
    model,
    data    = stan_data,
    iter    = iter,
    warmup  = warmup,
    chains  = chains,
    seed    = seed,
    control = list(adapt_delta = adapt_delta)
  )

  # Posterior PoS: mean of generated quantity pos_os
  pos_samples <- rstan::extract(fit, "pos_os")$pos_os
  pos_val     <- mean(pos_samples)

  # Convergence diagnostics
  fit_summary <- as.data.frame(rstan::summary(fit)$summary)
  fit_summary$parameter <- rownames(fit_summary)

  key_params <- fit_summary[fit_summary$parameter %in%
                              c("mu[1]", "mu[2]", "tau_os", "tau_pfs",
                                "rho", "theta_current[1]", "pos_os"), ]

  rhat_ok <- all(key_params[["Rhat"]] < 1.01, na.rm = TRUE)
  ess_ok  <- all(key_params[["n_eff"]] > 400,  na.rm = TRUE)

  list(
    fit     = fit,
    pos     = pos_val,
    summary = key_params,
    rhat_ok = rhat_ok,
    ess_ok  = ess_ok
  )
}
