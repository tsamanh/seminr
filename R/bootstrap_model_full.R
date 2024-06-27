#' seminr bootstrap_model_full Function
#' 
#' returning all estimates of each bootstrap
#' 
#' 
#'
#' @param seminr_model A fully estimated model with associated data, measurement model and structural model
#'
#' @param nboot A parameter specifying the number of bootstrap iterations to perform, default
#'        value is 500. If 0 then no bootstrapping is performed.
#'
#' @param cores A parameter specifying the maximum number of cores to use in the parallelization.
#'
#' @param seed A parameter to specify the seed for reproducibility of results. Default is NULL.
#'
#' @param ... A list of parameters passed on to the estimation method.
#'
#' @export

bootstrap_model_full <- function(seminr_model, nboot = 500, cores = NULL, seed = NULL, ...) {
  
  # Create a cluster for parallel computing
  suppressWarnings(ifelse(is.null(cores), cl <- parallel::makeCluster(parallel::detectCores(), setup_strategy = "sequential"), cl <- parallel::makeCluster(cores, setup_strategy = "sequential")))
  
  # prepare parameters for cluster export (model parameters)
  data <- seminr_model$rawdata
  measurement_model <- seminr_model$measurement_model
  structural_model <- seminr_model$smMatrix
  inner_weights <- seminr_model$inner_weights
  missing_value <- seminr_model$settings$missing_value
  maxIt <- seminr_model$settings$maxIt
  stopCriterion <- seminr_model$settings$stopCriterion
  missing <- seminr_model$settings$missing
  cross_loadings <- cor(data, seminr_model$construct_scores)
  
  # Check for and create random seed if NULL
  if (is.null(seed)) {seed <- sample.int(100000, size = 1)}
  
  
  # Function to perform a single bootstrap resample and run seminr::estimate_pls
  bootstrap_sample <- function(i, data, measurement_model, structural_model, inner_weights, original_model, original_cross_loadings) {
    set.seed(seed + i)
    resampled_data <- data[sample(1:nrow(data), replace = TRUE), ]
    boot_model <- tryCatch({
      seminr::estimate_pls(data = resampled_data,
                           measurement_model = measurement_model,
                           structural_model = structural_model,
                           inner_weights = inner_weights)
      
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(boot_model)) {
      return(list(data = resampled_data, boot_path_coef = NA, boot_outer_loadings = NA, boot_outer_weights = NA, boot_rSquared = NA))
    }
    
    path_coef <- boot_model$path_coef
    outer_loadings <- boot_model$outer_loadings
    outer_weights <- boot_model$outer_weights
    rSquared <- boot_model$rSquared

    
    return(list(boot_data = resampled_data,
                boot_path_coef = path_coef,
                boot_outer_loadings = outer_loadings,
                boot_outer_weights = outer_weights,
                boot_rSquared = rSquared))
  }
  
  # Export necessary variables and functions to the cluster
  parallel::clusterExport(cl, varlist = c("data",
                                          "measurement_model",
                                          "structural_model",
                                          "inner_weights",
                                          "bootstrap_sample"))
  
  # Perform bootstrap in parallel
  bootstrap_results <- parLapply(cl, 1:nboot, bootstrap_sample, data, measurement_model, structural_model, inner_weights)
  
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Separate the raw data and model estimates
  boot_raw_data <- lapply(bootstrap_results, function(res) res$boot_data)
  boot_path_coef <- lapply(bootstrap_results, function(res) res$boot_path_coef)
  boot_outer_loadings <- lapply(bootstrap_results, function(res) res$boot_outer_loadings)
  boot_outer_weights <- lapply(bootstrap_results, function(res) res$boot_outer_weights)
  boot_rSquared <- lapply(bootstrap_results, function(res) res$boot_rSquared)

  return(list(boot_raw_data = boot_raw_data,
              boot_path_coef = boot_path_coef,
              boot_outer_loadings = boot_outer_loadings,
              boot_outer_weights = boot_outer_weights,
              boot_rSquared = boot_rSquared))
}

library(parallel)
library(seminr)

boot_test <- bootstrap_model_full(pls_model_clc1, nboot = 10)
pls_model_clc1$inner_weights
