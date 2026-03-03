library(shiny)
library(ggplot2)
library(DT)

# ===========================================================================
# PosPredictor User Interface
# ===========================================================================

ui <- navbarPage(
  title = "PosPredictor: Using PFS to Predict OS",
  theme = NULL,
  id    = "main_nav",
  
  # Ensure MathJax is loaded for LaTeX formulas
  header = tags$head(withMathJax()),

  # ===========================================================================
  # Tab 1: Model Description
  # ===========================================================================
  tabPanel(
    "Model Description",
    fluidPage(
      h2("Hierarchical Random-Effects Meta-Analysis Model"),
      p("This application computes the posterior Probability of Success (PoS)",
        "for an ongoing OS readout, borrowing strength from historical trials",
        "through a three-level hierarchical (random-effects) meta-analytic model."),

      hr(),
      h3("Model Structure"),

      h4("Level 3 — Population Hyperparameters"),
      p("The population mean log(HR) vector follows a bivariate normal prior:"),
      p("$$\\boldsymbol{\\mu} = (\\mu_{OS}, \\mu_{PFS})^\\top, \\quad 
           \\mu_{OS} \\sim N(m_{OS}, s_{OS}^2), \\quad 
           \\mu_{PFS} \\sim N(m_{PFS}, s_{PFS}^2)$$"),
      p("The between-trial covariance matrix is:"),
      p("$$\\boldsymbol{\\Sigma} = \\begin{pmatrix} 
             \\tau_{OS}^2 & \\rho\\,\\tau_{OS}\\tau_{PFS} \\\\ 
             \\rho\\,\\tau_{OS}\\tau_{PFS} & \\tau_{PFS}^2 
           \\end{pmatrix}$$"),
      p("Priors: \\(\\tau_{OS}, \\tau_{PFS} \\sim \\text{half-Normal}\\);",
        " \\(\\rho\\) is modelled via Fisher-z transform: \\(\\tanh^{-1}(\\rho) \\sim N(0, \\sigma_z^2)\\)."),

      h4("Level 2 — Trial-Specific Effects"),
      p("For each trial \\(k = 1, \\ldots, K\\) (historical) and the current trial:"),
      p("$$\\boldsymbol{\\theta}_k = (\\theta_{k,OS}, \\theta_{k,PFS})^\\top 
           \\sim \\text{MVN}(\\boldsymbol{\\mu}, \\boldsymbol{\\Sigma})$$"),
      p("These represent the ", strong("true"), " treatment effects in each trial."),
      
      h4("Level 1 — Observed Data (Likelihood)"),
      p("Observed log(HR) pairs \\(\\mathbf{y}_k = (y_{k,OS}, y_{k,PFS})^\\top\\) are",
        " noisy measurements of the true effects:"),
      p("$$\\mathbf{y}_k \\sim \\text{MVN}(\\boldsymbol{\\theta}_k, W_k)$$"),
      
      h4("Probability of Success"),
      p("$$\\text{PoS}_{OS} = \\Pr(\\theta_{current,OS} < \\text{target}_{OS} \\mid \\text{data})$$"),

      hr(),
      h3("Prior Specifications"),
      tableOutput("prior_table"),

      hr(),
      h3("Computational Approaches"),
      tags$ul(
        tags$li(strong("Closed-Form (Tab 3):"), 
                " Conjugate bivariate-normal update for instant results."),
        tags$li(strong("Full Bayesian MCMC (Tab 4):"), 
                " Hamiltonian Monte Carlo via Stan. The model is automatically compiled",
                " and stored in a standard system directory to avoid permission issues.")
      )
    )
  ),

  # ===========================================================================
  # Tab 2: Data Input
  # ===========================================================================
  tabPanel(
    "Data Input",
    fluidPage(
      h2("Trial Data Input"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Historical Data Simulation"),
          sliderInput("K", "Number of historical trials", 5, 40, 27, step = 1),
          sliderInput("mu_os_true", "True μ_OS", -0.8, 0.2, -0.25, step = 0.01),
          sliderInput("mu_pfs_true", "True μ_PFS", -0.8, 0.2, -0.35, step = 0.01),
          sliderInput("tau_os_true", "True τ_OS", 0.01, 0.5, 0.15, step = 0.01),
          sliderInput("tau_pfs_true", "True τ_PFS", 0.01, 0.5, 0.18, step = 0.01),
          sliderInput("rho_true", "True ρ (between-trial)", -0.95, 0.95, 0.75, step = 0.05),
          numericInput("sim_seed", "Random seed", 42, min = 1),
          actionButton("simulate_btn", "Simulate Data", class = "btn-primary btn-block"),
          
          hr(),
          h4("Current Trial Interim"),
          numericInput("cur_y_os",  "Current log(HR) OS",  -0.30, step = 0.01),
          numericInput("cur_y_pfs", "Current log(HR) PFS", -0.40, step = 0.01),
          numericInput("cur_se_os",  "Current SE (OS)",   0.12, min = 0.01, step = 0.01),
          numericInput("cur_se_pfs", "Current SE (PFS)",  0.10, min = 0.01, step = 0.01),
          sliderInput("cur_within_corr", "Within-trial corr.", 0.2, 0.95, 0.65, step = 0.05),
          numericInput("target_os", "Success Threshold (log HR)", log(0.74), step = 0.01)
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            tabPanel("Scatter Plot", br(), plotOutput("scatter_plot", height = "500px")),
            tabPanel("Data Table", br(), DT::dataTableOutput("hist_table"))
          )
        )
      )
    )
  ),

  # ===========================================================================
  # Tab 3: Closed-Form PoS
  # ===========================================================================
  tabPanel(
    "Closed-Form PoS",
    fluidPage(
      h2("Conjugate MVN Model Results"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Prior Hyperparameters"),
          numericInput("cf_mu_os_prior_mean",  "Prior mean μ_OS",  -0.30, step = 0.01),
          numericInput("cf_mu_pfs_prior_mean", "Prior mean μ_PFS", -0.40, step = 0.01),
          numericInput("cf_sigma_prior_os",    "Prior SD σ_OS",     0.50, min = 0.01, step = 0.05),
          numericInput("cf_sigma_prior_pfs",   "Prior SD σ_PFS",    0.50, min = 0.01, step = 0.05),
          hr(),
          h4("Heterogeneity Settings"),
          numericInput("cf_tau_os",   "τ_OS (Between SD)",  0.15, min = 0.01, step = 0.01),
          numericInput("cf_tau_pfs",  "τ_PFS (Between SD)", 0.18, min = 0.01, step = 0.01),
          sliderInput("cf_rho_between", "ρ (Between Corr.)", -0.95, 0.95, 0.75, step = 0.05),
          br(),
          actionButton("run_cf_btn", "Compute PoS", class = "btn-success btn-block")
        ),
        mainPanel(
          width = 9,
          uiOutput("cf_pos_box"),
          fluidRow(
            column(6, plotOutput("cf_posterior_os_plot",  height = "350px")),
            column(6, plotOutput("cf_posterior_pfs_plot", height = "350px"))
          ),
          br(),
          h4("Meta-Analytic Summary"),
          tableOutput("cf_meta_summary")
        )
      )
    )
  ),

  # ===========================================================================
  # Tab 4: MCMC PoS (Stan)
  # ===========================================================================
  tabPanel(
    "MCMC PoS (Stan)",
    fluidPage(
      h2("Full Bayesian MCMC via Stan"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Bayesian Priors"),
          numericInput("mc_mu_os_mean",  "Prior mean μ_OS",  -0.30, step = 0.01),
          numericInput("mc_mu_os_sd",    "Prior SD μ_OS",     0.50, min = 0.01, step = 0.05),
          numericInput("mc_mu_pfs_mean", "Prior mean μ_PFS", -0.40, step = 0.01),
          numericInput("mc_mu_pfs_sd",   "Prior SD μ_PFS",    0.50, min = 0.01, step = 0.05),
          numericInput("mc_tau_os_sd",   "τ_OS Prior SD",  0.25, min = 0.01, step = 0.05),
          numericInput("mc_tau_pfs_sd",  "τ_PFS Prior SD", 0.25, min = 0.01, step = 0.05),
          numericInput("mc_rho_z_sd",    "Fisher-z SD", 1.50, min = 0.10, step = 0.10),
          
          hr(),
          h4("MCMC Engine"),
          numericInput("mc_iter",    "Total Iterations",  2000, min = 500,  step = 500),
          numericInput("mc_chains",  "Chains", 4, min = 1, step = 1),
          numericInput("mc_seed",    "MCMC Seed", 123, min = 1),
          
          hr(),
          h4("Model Management"),
          helpText("Model files are automatically managed in the app data directory."),
          actionButton("compile_btn", "Compile Model", class = "btn-warning btn-block"),
          br(),
          actionButton("run_mcmc_btn", "Run MCMC Analysis", class = "btn-danger btn-block btn-lg")
        ),
        mainPanel(
          width = 9,
          uiOutput("mcmc_pos_box"),
          uiOutput("mcmc_convergence_ui"),
          br(),
          tabsetPanel(
            tabPanel("Trace Plots", br(), plotOutput("mcmc_trace_plot", height = "500px")),
            tabPanel("Densities", br(), plotOutput("mcmc_density_plot", height = "500px")),
            tabPanel("Pairs", br(), plotOutput("mcmc_pairs_plot", height = "500px")),
            tabPanel("Diagnostics", br(), tableOutput("mcmc_diag_table"))
          )
        )
      )
    )
  )
)
