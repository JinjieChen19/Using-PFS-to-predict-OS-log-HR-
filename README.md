# Using PFS to Predict OS — log(HR) Surrogate Analysis

This repository contains two complementary tools for the meta-analytic
surrogate-endpoint analysis of **Progression-Free Survival (PFS) as a
surrogate for Overall Survival (OS)** using trial-level log(Hazard Ratio) data.

---

## Repository Structure

```
.
├── analysis.py               # Stand-alone Python weighted-regression script
├── requirements.txt          # Python dependencies
└── PosPredictor/             # R package — hierarchical Bayesian PoS app
    ├── DESCRIPTION
    ├── NAMESPACE
    ├── LICENSE
    ├── R/
    │   ├── simulate_data.R   # Simulate historical trial data
    │   ├── closed_form.R     # Closed-form conjugate MVN PoS
    │   ├── compile_model.R   # Compile / load Stan model + MCMC PoS
    │   └── run_app.R         # Launch the Shiny application
    └── inst/
        ├── shiny/
        │   ├── ui.R          # Four-tab Shiny UI
        │   └── server.R      # Shiny server logic
        └── stan/
            └── hierarchical_model.stan  # Three-level Stan model
```

---

## 1 · Python Script (`analysis.py`)

A self-contained script that fits a **weighted linear regression** of
`log(OS HR)` on `log(PFS HR)` across 20 published oncology trials and
produces a publication-ready scatter plot with a 95 % confidence band.

### Requirements

```bash
pip install -r requirements.txt
```

### Run

```bash
python analysis.py
```

Output: regression summary printed to the console and a plot saved to
`pfs_vs_os_log_hr.png`.

---

## 2 · R Package (`PosPredictor`)

An R package that ships a **four-tab Shiny application** for computing the
posterior **Probability of Success (PoS)** of an OS endpoint, borrowing
strength from historical trials through a three-level hierarchical
random-effects meta-analytic model.

### Installation

```r
# Install dependencies first
install.packages(c("shiny", "ggplot2", "MASS", "mvtnorm",
                   "DT", "rstan", "bayesplot"))

# Install PosPredictor from the local source
install.packages("PosPredictor", repos = NULL, type = "source")
```

### Launch the App

```r
library(PosPredictor)
run_app()
```

### Key Exported Functions

| Function | Description |
|---|---|
| `simulate_historical_data()` | Generate synthetic trial-level log(HR) data |
| `compute_pos_closed_form()` | Closed-form conjugate bivariate-normal PoS |
| `compile_stan_model()` | Compile the hierarchical Stan model to an RDS |
| `load_stan_model()` | Load the pre-compiled Stan model |
| `compute_pos_mcmc()` | Full Bayesian MCMC PoS via Stan |
| `run_app()` | Launch the interactive four-tab Shiny application |

### Statistical Model

The package implements a **three-level hierarchical model**:

* **Level 3 — Population hyperparameters**  
  μ = (μ_OS, μ_PFS) ~ MVN prior; τ_OS, τ_PFS ~ half-Normal; ρ via Fisher-z.

* **Level 2 — Trial-specific true effects**  
  θ_k ~ MVN(μ, Σ), non-centered parameterisation for efficient MCMC.

* **Level 1 — Observed data (likelihood)**  
  y_k ~ MVN(θ_k, W_k), where W_k is the known within-trial covariance.

**PoS** = P(θ_current,OS < target_OS | data).

Two computational backends are provided:
* **Closed-form** — conjugate sequential Bayesian update (instant).
* **Full Bayesian MCMC** — Hamiltonian Monte Carlo via Stan (HMC, 4 chains,
  adapt_delta = 0.97); priors are passed as data so no recompilation is needed
  when hyperparameters change.

---

## License

MIT — see [`PosPredictor/LICENSE`](PosPredictor/LICENSE).
