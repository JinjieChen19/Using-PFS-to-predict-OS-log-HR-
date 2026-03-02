#' Compute PoS via Closed-Form Conjugate Multivariate-Normal Model
#'
#' Uses a conjugate bivariate-normal Bayesian update to compute the posterior
#' probability of success (PoS) for the current trial's OS endpoint.
#'
#' Model summary:
#' \itemize{
#'   \item Prior: \eqn{\mu \sim MVN(\mu_0, \Sigma_0)}
#'   \item Marginal for trial k: \eqn{y_k \sim MVN(\mu, \Sigma_{between} + W_k)}
#'   \item After observing historical data: posterior \eqn{\mu | data_{hist}}
#'   \item Predictive prior for current trial: \eqn{\theta_{cur} \sim MVN(\hat\mu, \Sigma_{between} + \hat V_\mu)}
#'   \item Update with current data: \eqn{y_{cur} \sim MVN(\theta_{cur}, W_{cur})}
#'   \item PoS = \eqn{P(\theta_{cur,OS} < target_{OS})}
#' }
#'
#' @param hist_data data.frame from \code{\link{simulate_historical_data}}.
#' @param current_y_os Numeric. Current trial observed log(HR) for OS.
#' @param current_y_pfs Numeric. Current trial observed log(HR) for PFS.
#' @param current_se_os Numeric. Current trial SE for OS log(HR).
#' @param current_se_pfs Numeric. Current trial SE for PFS log(HR).
#' @param current_within_corr Numeric. Within-trial correlation for current trial.
#' @param target_os Numeric. Success threshold for OS log(HR) (e.g. log(0.74)).
#' @param mu_prior Numeric vector length 2. Prior means for (OS, PFS) log(HR).
#' @param sigma_prior Numeric vector length 2. Prior SDs for (OS, PFS) log(HR).
#' @param tau Numeric vector length 2. Between-trial SDs for (OS, PFS).
#' @param rho_between Numeric. Between-trial correlation.
#' @return A named list with: pos, mu_posterior, sigma_posterior,
#'   mu_prior_hist, sigma_prior_hist, credible_interval_os.
#' @export
compute_pos_closed_form <- function(hist_data,
                                     current_y_os,
                                     current_y_pfs,
                                     current_se_os,
                                     current_se_pfs,
                                     current_within_corr = 0.65,
                                     target_os   = log(0.74),
                                     mu_prior    = c(-0.30, -0.40),
                                     sigma_prior = c(0.50, 0.50),
                                     tau         = c(0.15, 0.18),
                                     rho_between = 0.75) {

  # Between-trial covariance
  Sigma_between <- matrix(
    c(tau[1]^2,               rho_between * tau[1] * tau[2],
      rho_between * tau[1] * tau[2], tau[2]^2),
    nrow = 2
  )

  # --- Step 1: Update population prior with historical data ---
  # Prior precision
  Sigma_prior_mat <- diag(sigma_prior^2)
  Lambda_post <- solve(Sigma_prior_mat)
  eta_post    <- Lambda_post %*% mu_prior

  K <- nrow(hist_data)
  for (k in seq_len(K)) {
    W_k <- matrix(
      c(hist_data$se_os[k]^2,
        hist_data$within_corr[k] * hist_data$se_os[k] * hist_data$se_pfs[k],
        hist_data$within_corr[k] * hist_data$se_os[k] * hist_data$se_pfs[k],
        hist_data$se_pfs[k]^2),
      nrow = 2
    )
    # Marginal covariance for trial k: Sigma_between + W_k
    V_k     <- Sigma_between + W_k
    V_k_inv <- solve(V_k)
    Lambda_post <- Lambda_post + V_k_inv
    eta_post    <- eta_post + V_k_inv %*% c(hist_data$log_hr_os[k],
                                             hist_data$log_hr_pfs[k])
  }

  Sigma_post_hist <- solve(Lambda_post)
  mu_post_hist    <- Sigma_post_hist %*% eta_post

  # --- Step 2: Predictive prior for current trial's true effect ---
  # theta_current | hist_data ~ MVN(mu_post_hist, Sigma_between + Sigma_post_hist)
  Sigma_pred_prior <- Sigma_between + Sigma_post_hist

  # --- Step 3: Update with current trial data ---
  W_current <- matrix(
    c(current_se_os^2,
      current_within_corr * current_se_os * current_se_pfs,
      current_within_corr * current_se_os * current_se_pfs,
      current_se_pfs^2),
    nrow = 2
  )

  Lambda_final <- solve(Sigma_pred_prior) + solve(W_current)
  eta_final    <- solve(Sigma_pred_prior) %*% mu_post_hist +
                  solve(W_current) %*% c(current_y_os, current_y_pfs)

  Sigma_final <- solve(Lambda_final)
  mu_final    <- Sigma_final %*% eta_final

  # --- Step 4: PoS = P(theta_OS < target_os) ---
  pos_os <- stats::pnorm(target_os,
                         mean = mu_final[1],
                         sd   = sqrt(Sigma_final[1, 1]))

  # 95% credible interval for theta_OS
  ci_os <- stats::qnorm(c(0.025, 0.975),
                        mean = mu_final[1],
                        sd   = sqrt(Sigma_final[1, 1]))

  list(
    pos                  = pos_os,
    mu_posterior         = as.numeric(mu_final),
    sigma_posterior      = Sigma_final,
    mu_prior_hist        = as.numeric(mu_post_hist),
    sigma_prior_hist     = Sigma_post_hist,
    credible_interval_os = ci_os
  )
}
