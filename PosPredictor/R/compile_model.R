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
  # 2. Try RcppParallel (common on macOS/Linux for rstan)
  rcpp_inc <- system.file("include", package = "RcppParallel")
  if (nchar(rcpp_inc) && file.exists(file.path(rcpp_inc, "boost"))) {
    return(rcpp_inc)
  }
  # 3. Try common system Boost paths
  candidates <- c("/usr/include", "/usr/local/include", "/opt/homebrew/include")
  for (p in candidates) {
    if (file.exists(file.path(p, "boost"))) return(p)
  }
  ""  # not found
}

#' Get Default Path for Compiled Model
#' 
#' @return Character string of the RDS file path.
#' @keywords internal
.get_default_rds_path <- function() {
  # Use R's standard cross-platform user data directory
  # This avoids permission issues with hidden folders like ~/.PosPredictor
  data_dir <- tools::R_user_dir("PosPredictor", which = "data")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(file.path(data_dir, "compiled_model.rds"))
}

#' Compile the Hierarchical Stan Model and Save as RDS
#'
#' @param output_path Character. Path where the compiled model RDS will be saved.
#' @param verbose Logical. Print compilation progress (default TRUE).
#' @return Invisibly returns the path to the saved RDS file.
#' @export
compile_stan_model <- function(output_path = NULL, verbose = TRUE) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required. Install it with: install.packages('rstan')")
  }

  stan_file <- system.file("stan", "hierarchical_model.stan", package = "PosPredictor")
  if (!nchar(stan_file)) {
    stop("Stan model file not found in PosPredictor installation.")
  }

  if (is.null(output_path)) {
    output_path <- .get_default_rds_path()
  }

  # Ensure the parent directory exists before saving
  dir_path <- dirname(output_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
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
#' @param rds_path Character. Path to the compiled model RDS.
#' @param verbose Logical. Verbosity (default TRUE).
#' @return A \code{stanmodel} object.
#' @export
load_stan_model <- function(rds_path = NULL, verbose = TRUE) {
  if (is.null(rds_path)) {
    rds_path <- .get_default_rds_path()
  }
  
  if (!file.exists(rds_path)) {
    if (verbose) message("Compiled model not found. Compiling now...")
    compile_stan_model(output_path = rds_path, verbose = verbose)
  }
  
  # Final check to ensure file was actually created and is readable
  if (!file.exists(rds_path)) {
    stop("Failed to create or access the compiled model at: ", rds_path)
  }
  
  readRDS(rds_path)
}
