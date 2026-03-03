library(shiny)
library(ggplot2)
library(MASS)
library(DT)

# --- 1. Dynamic Package Loading ---
# Check if the package is installed to use its namespace, 
# otherwise source files from the development directory.
pkg_available <- requireNamespace("PosPredictor", quietly = TRUE)

if (pkg_available) {
  simulate_historical_data <- PosPredictor::simulate_historical_data
  compute_pos_closed_form  <- PosPredictor::compute_pos_closed_form
  compile_stan_model       <- PosPredictor::compile_stan_model
  load_stan_model          <- PosPredictor::load_stan_model
  compute_pos_mcmc         <- PosPredictor::compute_pos_mcmc
} else {
  # Compatibility mode for running from source code
  pkg_r_dir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
  if (dir.exists(pkg_r_dir)) {
    for (f in list.files(pkg_r_dir, pattern = "\\.R$", full.names = TRUE)) {
      source(f, local = FALSE)
    }
  }
}

# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # ---- Reactive Data Storage ----
  hist_data <- reactiveVal(NULL)

  # Tab 1: Display the prior parameters
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
  # Tab 2: Data Input & Simulation
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

  # Initialize data on app startup
  isolate({
    hist_data(simulate_historical_data())
  })

  # Scatter plot for PFS vs OS log(HR)
  output$scatter_plot <- renderPlot({
    df <- hist_data()
    req(df)

    cur <- data.frame(
      log_hr_pfs   = input$cur_y_pfs,
      log_hr_os    = input$cur_y_os,
      n_os_events  = 0,
      trial        = "Current Trial"
    )

    # Calculate weighted linear regression line
    wt <- df$n_os_events
    if (sum(wt) > 0) {
      fit_lm <- lm(log_hr_os ~ log_hr_pfs, data = df, weights = wt)
      slope_val <- round(coef(fit_lm)[2], 3)
      r2_val    <- round(summary(fit_lm)$r.squared, 3)
      lm_label  <- paste0("Slope = ", slope_val, ",  R\u00b2 = ", r2_val)
    } else {
      lm_label <- ""
    }

    # Plotting using ggplot2
    ggplot(df, aes(x = log_hr_pfs, y = log_hr_os, size = n_os_events)) +
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
  })

  # =========================================================================
  # Tab 4: MCMC PoS (Stan) - Final logic for Stan model
  # =========================================================================
  mcmc_result <- reactiveVal(NULL)

  # Trigger Stan compilation
  observeEvent(input$compile_btn, {
    # Using NULL rds_path triggers the internal smart path discovery
    withProgress(message = "Compiling Stan model...", value = 0.5, {
      tryCatch(
        compile_stan_model(verbose = FALSE), 
        error = function(e) {
          showNotification(paste("Compilation error:", conditionMessage(e)),
                           type = "error", duration = 10)
        }
      )
    })
    showNotification("Stan model compiled successfully!", type = "message")
  })

  # Run MCMC sampling
  observeEvent(input$run_mcmc_btn, {
    df <- hist_data()
    req(df)

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
          rds_path            = NULL, # Forces usage of the intelligent default path
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
}
