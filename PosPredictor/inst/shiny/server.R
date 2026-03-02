library(shiny)
library(ggplot2)
library(MASS)
library(DT)

# Source package functions: use installed package namespace if available,
# otherwise source R/ files relative to the shiny/ directory itself.
pkg_available <- requireNamespace("PosPredictor", quietly = TRUE)
if (pkg_available) {
  simulate_historical_data <- PosPredictor::simulate_historical_data
  compute_pos_closed_form  <- PosPredictor::compute_pos_closed_form
  compile_stan_model       <- PosPredictor::compile_stan_model
  load_stan_model          <- PosPredictor::load_stan_model
  compute_pos_mcmc         <- PosPredictor::compute_pos_mcmc
} else {
  # Running from source: shiny/ is inside inst/shiny/, so R/ is two levels up
  pkg_r_dir <- normalizePath(
    file.path(dirname(sys.frame(1)$ofile), "..", "..", "R"),
    mustWork = FALSE
  )
  if (dir.exists(pkg_r_dir)) {
    for (f in list.files(pkg_r_dir, pattern = "\\.R$", full.names = TRUE)) {
      source(f, local = FALSE)
    }
  }
}

# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # ---- Reactive data store --------------------------------------------------
  hist_data <- reactiveVal(NULL)

  # Prior specification table (Tab 1)
  output$prior_table <- renderTable({
    data.frame(
      Parameter = c("mu_os", "mu_pfs", "tau_os", "tau_pfs", "rho"),
      Meaning   = c("Pop. mean log(HR) OS", "Pop. mean log(HR) PFS",
                    "Between-trial SD for OS", "Between-trial SD for PFS",
                    "Between-trial correlation"),
      Prior     = c("Gaussian", "Gaussian", "half-Normal", "half-Normal",
                    "Fisher-z transform"),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

  # =========================================================================
  # Tab 2: Data Input
  # =========================================================================
  observeEvent(input$simulate_btn, {
    req(input$K, input$mu_os_true, input$mu_pfs_true,
        input$tau_os_true, input$tau_pfs_true,
        input$rho_true, input$sim_seed)
    df <- simulate_historical_data(
      K       = input$K,
      mu_os   = input$mu_os_true,
      mu_pfs  = input$mu_pfs_true,
      tau_os  = input$tau_os_true,
      tau_pfs = input$tau_pfs_true,
      rho     = input$rho_true,
      seed    = input$sim_seed
    )
    hist_data(df)
  })

  # Auto-simulate on startup
  isolate({
    hist_data(simulate_historical_data())
  })

  output$scatter_plot <- renderPlot({
    df <- hist_data()
    req(df)

    # Current trial point
    cur <- data.frame(
      log_hr_pfs   = input$cur_y_pfs,
      log_hr_os    = input$cur_y_os,
      n_os_events  = 0,
      trial        = "Current Trial"
    )

    # Weighted regression
    wt <- df$n_os_events
    if (sum(wt) > 0) {
      fit_lm <- lm(log_hr_os ~ log_hr_pfs, data = df, weights = wt)
      slope_val <- round(coef(fit_lm)[2], 3)
      r2_val    <- round(summary(fit_lm)$r.squared, 3)
      lm_label  <- paste0("Slope = ", slope_val, ",  R\u00b2 = ", r2_val)
    } else {
      lm_label <- ""
    }

    p <- ggplot(df, aes(x = log_hr_pfs, y = log_hr_os, size = n_os_events)) +
      geom_point(aes(color = "Historical Trials"), alpha = 0.75) +
      geom_smooth(method = "lm", formula = y ~ x,
                  aes(weight = n_os_events),
                  se = TRUE, color = "firebrick", linewidth = 0.8) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "grey50", linewidth = 0.6) +
      geom_point(data = cur,
                 aes(x = log_hr_pfs, y = log_hr_os, color = "Current Trial"),
                 size = 6, shape = 3, stroke = 2) +
      geom_hline(yintercept = input$target_os,
                 linetype = "dotted", color = "purple", linewidth = 0.7) +
      scale_color_manual(
        name   = "",
        values = c("Historical Trials" = "steelblue",
                   "Current Trial"     = "darkorange")
      ) +
      scale_size_continuous(name = "OS events", range = c(2, 8)) +
      labs(
        title    = "log(HR) for PFS vs log(HR) for OS across Trials",
        subtitle = lm_label,
        x        = "log(HR) — PFS",
        y        = "log(HR) — OS"
      ) +
      theme_bw(base_size = 13) +
      theme(legend.position = "bottom")
    print(p)
  })

  output$hist_table <- DT::renderDataTable({
    df <- hist_data()
    req(df)
    out <- df[, c("trial", "log_hr_os", "log_hr_pfs",
                  "se_os", "se_pfs", "within_corr",
                  "n_os_events", "n_pfs_events")]
    out[, 2:7] <- round(out[, 2:7], 4)
    out
  }, options = list(pageLength = 10))

  # =========================================================================
  # Tab 3: Closed-Form PoS
  # =========================================================================
  cf_result <- reactiveVal(NULL)

  observeEvent(input$run_cf_btn, {
    df <- hist_data()
    req(df, input$cur_y_os, input$cur_y_pfs,
        input$cur_se_os, input$cur_se_pfs)

    res <- tryCatch(
      compute_pos_closed_form(
        hist_data           = df,
        current_y_os        = input$cur_y_os,
        current_y_pfs       = input$cur_y_pfs,
        current_se_os       = input$cur_se_os,
        current_se_pfs      = input$cur_se_pfs,
        current_within_corr = input$cur_within_corr,
        target_os           = input$target_os,
        mu_prior            = c(input$cf_mu_os_prior_mean,
                                input$cf_mu_pfs_prior_mean),
        sigma_prior         = c(input$cf_sigma_prior_os,
                                input$cf_sigma_prior_pfs),
        tau                 = c(input$cf_tau_os, input$cf_tau_pfs),
        rho_between         = input$cf_rho_between
      ),
      error = function(e) {
        showNotification(paste("Error:", conditionMessage(e)), type = "error")
        NULL
      }
    )
    cf_result(res)
  })

  output$cf_pos_box <- renderUI({
    res <- cf_result()
    if (is.null(res)) {
      return(p("Press 'Compute Closed-Form PoS' to see results.",
               style = "color:grey;"))
    }
    pos_pct <- round(res$pos * 100, 1)
    ci      <- round(res$credible_interval_os, 3)
    color   <- if (pos_pct >= 70) "#27ae60" else if (pos_pct >= 50) "#f39c12" else "#e74c3c"

    div(
      style = paste0("background:", color, "; color:white; padding:20px;",
                     " border-radius:8px; margin-bottom:10px;"),
      h3(paste0("Closed-Form PoS = ", pos_pct, "%"),
         style = "margin:0; font-size:28px;"),
      p(paste0("Posterior mean θ_OS = ",
               round(res$mu_posterior[1], 4),
               " (95% CI: ", ci[1], " to ", ci[2], ")"),
        style = "margin:5px 0 0 0;"),
      p(paste0("Success threshold: target log(HR) = ", round(input$target_os, 3),
               " (HR ≤ ", round(exp(input$target_os), 3), ")"),
        style = "margin:5px 0 0 0;")
    )
  })

  output$cf_posterior_os_plot <- renderPlot({
    res <- cf_result()
    req(res)
    mu_os  <- res$mu_posterior[1]
    sd_os  <- sqrt(res$sigma_posterior[1, 1])
    x_seq  <- seq(mu_os - 4 * sd_os, mu_os + 4 * sd_os, length.out = 400)
    df_plt <- data.frame(x = x_seq, y = dnorm(x_seq, mu_os, sd_os))

    ggplot(df_plt, aes(x, y)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_area(data = subset(df_plt, x < input$target_os),
                fill = "steelblue", alpha = 0.35) +
      geom_vline(xintercept = input$target_os,
                 color = "red", linetype = "dashed") +
      geom_vline(xintercept = mu_os, color = "navy", linetype = "dotted") +
      labs(title = "Posterior Distribution: θ_current, OS",
           subtitle = paste0("Shaded area = PoS = ",
                             round(pnorm(input$target_os, mu_os, sd_os) * 100, 1), "%"),
           x = "θ_OS (log HR)", y = "Density") +
      theme_bw(base_size = 12)
  })

  output$cf_posterior_pfs_plot <- renderPlot({
    res <- cf_result()
    req(res)
    mu_pfs  <- res$mu_posterior[2]
    sd_pfs  <- sqrt(res$sigma_posterior[2, 2])
    x_seq   <- seq(mu_pfs - 4 * sd_pfs, mu_pfs + 4 * sd_pfs, length.out = 400)
    df_plt  <- data.frame(x = x_seq, y = dnorm(x_seq, mu_pfs, sd_pfs))

    ggplot(df_plt, aes(x, y)) +
      geom_line(color = "darkorange", linewidth = 1.2) +
      geom_vline(xintercept = mu_pfs, color = "darkorange4", linetype = "dotted") +
      labs(title = "Posterior Distribution: θ_current, PFS",
           x = "θ_PFS (log HR)", y = "Density") +
      theme_bw(base_size = 12)
  })

  output$cf_meta_summary <- renderTable({
    res <- cf_result()
    req(res)
    data.frame(
      Parameter = c("μ_OS (meta-analytic)",  "μ_PFS (meta-analytic)",
                    "θ_OS (current, posterior)", "θ_PFS (current, posterior)"),
      Estimate  = round(c(res$mu_prior_hist[1], res$mu_prior_hist[2],
                          res$mu_posterior[1],  res$mu_posterior[2]), 4),
      SD        = round(c(sqrt(res$sigma_prior_hist[1, 1]),
                          sqrt(res$sigma_prior_hist[2, 2]),
                          sqrt(res$sigma_posterior[1, 1]),
                          sqrt(res$sigma_posterior[2, 2])), 4),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  # =========================================================================
  # Tab 4: MCMC PoS (Stan)
  # =========================================================================
  mcmc_result <- reactiveVal(NULL)

  # Compile Stan model
  observeEvent(input$compile_btn, {
    rds <- if (nchar(trimws(input$rds_path)) > 0) input$rds_path else NULL
    withProgress(message = "Compiling Stan model...", value = 0.5, {
      tryCatch(
        compile_stan_model(output_path = rds, verbose = FALSE),
        error = function(e) {
          showNotification(paste("Compilation error:", conditionMessage(e)),
                           type = "error", duration = 10)
        }
      )
    })
    showNotification("Stan model compiled and saved!", type = "message")
  })

  # Run MCMC
  observeEvent(input$run_mcmc_btn, {
    df <- hist_data()
    req(df)
    rds_path_val <- if (nchar(trimws(input$rds_path)) > 0) input$rds_path else NULL

    withProgress(message = "Running MCMC (this may take a few minutes)...",
                 value = 0.3, {
      res <- tryCatch(
        compute_pos_mcmc(
          hist_data           = df,
          current_y_os        = input$cur_y_os,
          current_y_pfs       = input$cur_y_pfs,
          current_se_os       = input$cur_se_os,
          current_se_pfs      = input$cur_se_pfs,
          current_within_corr = input$cur_within_corr,
          target_os           = input$target_os,
          mu_os_prior_mean    = input$mc_mu_os_mean,
          mu_os_prior_sd      = input$mc_mu_os_sd,
          mu_pfs_prior_mean   = input$mc_mu_pfs_mean,
          mu_pfs_prior_sd     = input$mc_mu_pfs_sd,
          tau_os_prior_sd     = input$mc_tau_os_sd,
          tau_pfs_prior_sd    = input$mc_tau_pfs_sd,
          rho_z_prior_sd      = input$mc_rho_z_sd,
          iter                = input$mc_iter,
          warmup              = input$mc_warmup,
          chains              = input$mc_chains,
          adapt_delta         = input$mc_adapt,
          rds_path            = rds_path_val,
          seed                = input$mc_seed
        ),
        error = function(e) {
          showNotification(paste("MCMC error:", conditionMessage(e)),
                           type = "error", duration = 20)
          NULL
        }
      )
      setProgress(1)
    })
    mcmc_result(res)
  })

  output$mcmc_pos_box <- renderUI({
    res <- mcmc_result()
    if (is.null(res)) {
      return(p("Press 'Run MCMC' to compute the Bayesian PoS.",
               style = "color:grey;"))
    }
    pos_pct <- round(res$pos * 100, 1)
    color   <- if (pos_pct >= 70) "#27ae60" else if (pos_pct >= 50) "#f39c12" else "#e74c3c"
    div(
      style = paste0("background:", color, "; color:white; padding:20px;",
                     " border-radius:8px; margin-bottom:10px;"),
      h3(paste0("MCMC PoS = ", pos_pct, "%"), style = "margin:0; font-size:28px;"),
      p(paste0("Based on ", input$mc_chains, " chains × ",
               input$mc_iter - input$mc_warmup, " post-warmup draws"),
        style = "margin:5px 0 0 0;")
    )
  })

  output$mcmc_convergence_ui <- renderUI({
    res <- mcmc_result()
    req(res)
    rhat_ok_txt <- if (res$rhat_ok) "✓ Rhat < 1.01" else "✗ Rhat ≥ 1.01"
    ess_ok_txt  <- if (res$ess_ok)  "✓ ESS > 400"   else "✗ ESS ≤ 400"
    rhat_col    <- if (res$rhat_ok) "green" else "red"
    ess_col     <- if (res$ess_ok)  "green" else "red"
    fluidRow(
      column(3,
        div(style = paste0("border:1px solid ", rhat_col,
                           "; padding:10px; border-radius:6px;"),
            h5(rhat_ok_txt, style = paste0("color:", rhat_col, "; margin:0;")))
      ),
      column(3,
        div(style = paste0("border:1px solid ", ess_col,
                           "; padding:10px; border-radius:6px;"),
            h5(ess_ok_txt,  style = paste0("color:", ess_col,  "; margin:0;")))
      )
    )
  })

  output$mcmc_trace_plot <- renderPlot({
    res <- mcmc_result()
    req(res)
    pars <- c("mu[1]", "mu[2]", "tau_os", "tau_pfs", "rho")
    bayesplot::mcmc_trace(rstan::extract(res$fit, permuted = FALSE), pars = pars)
  })

  output$mcmc_density_plot <- renderPlot({
    res <- mcmc_result()
    req(res)
    pars <- c("mu[1]", "mu[2]", "tau_os", "tau_pfs", "rho",
              "theta_current[1]", "pos_os")
    bayesplot::mcmc_dens_overlay(
      rstan::extract(res$fit, permuted = FALSE), pars = pars
    )
  })

  output$mcmc_pairs_plot <- renderPlot({
    res <- mcmc_result()
    req(res)
    bayesplot::mcmc_pairs(
      res$fit,
      pars        = c("mu[1]", "mu[2]", "tau_os", "tau_pfs", "rho"),
      off_diag_args = list(size = 0.5, alpha = 0.3)
    )
  })

  output$mcmc_diag_table <- renderTable({
    res <- mcmc_result()
    req(res)
    df <- res$summary[, c("parameter", "mean", "sd", "2.5%", "50%",
                          "97.5%", "n_eff", "Rhat")]
    df[, -1] <- round(df[, -1], 4)
    df
  }, striped = TRUE, bordered = TRUE)
}
