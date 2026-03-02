#' Simulate Historical Trial Data
#'
#' Generates synthetic trial-level log(HR) data for OS and PFS following the
#' three-level hierarchical model: population parameters -> trial-specific
#' effects -> observed data.
#'
#' @param K Integer. Number of historical trials (default 27).
#' @param mu_os Numeric. True population mean log(HR) for OS (default -0.25).
#' @param mu_pfs Numeric. True population mean log(HR) for PFS (default -0.35).
#' @param tau_os Numeric. Between-trial SD for OS (default 0.15).
#' @param tau_pfs Numeric. Between-trial SD for PFS (default 0.18).
#' @param rho Numeric. Between-trial correlation (default 0.75).
#' @param seed Integer. Random seed for reproducibility (default 42).
#' @return A data.frame with columns: trial, log_hr_os, log_hr_pfs,
#'   se_os, se_pfs, within_corr, n_os_events, n_pfs_events,
#'   true_theta_os, true_theta_pfs.
#' @export
simulate_historical_data <- function(K = 27,
                                      mu_os  = -0.25,
                                      mu_pfs = -0.35,
                                      tau_os  = 0.15,
                                      tau_pfs = 0.18,
                                      rho     = 0.75,
                                      seed    = 42) {
  set.seed(seed)

  # Between-trial covariance matrix (Level 3 -> Level 2)
  Sigma_between <- matrix(
    c(tau_os^2,             rho * tau_os * tau_pfs,
      rho * tau_os * tau_pfs, tau_pfs^2),
    nrow = 2
  )

  # Level 2: Draw true trial-specific effects theta_k ~ MVN(mu, Sigma)
  mu_vec  <- c(mu_os, mu_pfs)
  theta_k <- MASS::mvrnorm(K, mu = mu_vec, Sigma = Sigma_between)

  # Level 1: Within-trial parameters (event counts drive SE)
  n_os_events  <- sample(180:300, K, replace = TRUE)
  n_pfs_events <- sample(250:380, K, replace = TRUE)

  # SE ≈ 2 / sqrt(D) for log(HR) from balanced trials
  se_os  <- 2 / sqrt(n_os_events)  + runif(K, 0, 0.02)
  se_pfs <- 2 / sqrt(n_pfs_events) + runif(K, 0, 0.02)

  # Within-trial correlation (typically 0.55-0.75)
  within_corr <- runif(K, min = 0.55, max = 0.75)

  # Level 1: Observed data y_k ~ MVN(theta_k, W_k)
  log_hr_os  <- numeric(K)
  log_hr_pfs <- numeric(K)
  for (k in seq_len(K)) {
    W_k <- matrix(
      c(se_os[k]^2,
        within_corr[k] * se_os[k] * se_pfs[k],
        within_corr[k] * se_os[k] * se_pfs[k],
        se_pfs[k]^2),
      nrow = 2
    )
    obs_k <- MASS::mvrnorm(1, mu = theta_k[k, ], Sigma = W_k)
    log_hr_os[k]  <- obs_k[1]
    log_hr_pfs[k] <- obs_k[2]
  }

  data.frame(
    trial          = paste0("Trial_", seq_len(K)),
    log_hr_os      = log_hr_os,
    log_hr_pfs     = log_hr_pfs,
    se_os          = se_os,
    se_pfs         = se_pfs,
    within_corr    = within_corr,
    n_os_events    = n_os_events,
    n_pfs_events   = n_pfs_events,
    true_theta_os  = theta_k[, 1],
    true_theta_pfs = theta_k[, 2],
    stringsAsFactors = FALSE
  )
}
