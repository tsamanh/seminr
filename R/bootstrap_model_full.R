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
    
    
    boot_HTMT <- HTMT(boot_model)
    boot_total_effects <- total_effects(boot_model$path_coef)
    path_coef <- boot_model$path_coef
    outer_loadings <- boot_model$outer_loadings
    outer_weights <- boot_model$outer_weights
    rSquared <- boot_model$rSquared
    scores <- boot_model$construct_scores
    boot_cross_loadings <- cor(resampled_data, scores)
    
    # Sign change: benchmark on the original sample for the sign
    
    sc_path_coef <- boot_sign_change(path_coef, original_model$path_coef)
    sc_outer_loadings <- boot_sign_change(outer_loadings, original_model$outer_loadings)
    sc_outer_weights <- boot_sign_change(outer_weights, original_model$outer_weights)
    sc_cross_loadings <- boot_sign_change(boot_cross_loadings, original_cross_loadings)
    sc_construct_scores <- scores
    
    # Store the bootstrap index that is flipping, if it doesn't flip, return 0
    if (any(sc_path_coef != path_coef)) {
      sc_index <- i
    } else {
      sc_index <- 0
    }
    
    #Construct-level sign change: compare sum of the sums vs difference of the sums
    cl_change <- boot_construct_level_change(boot_model, original_model)
    cl_path_coef <- cl_change$cl_boot_path_coef
    cl_outer_loadings <- cl_change$cl_boot_outer_loadings
    cl_outer_weights <- cl_change$cl_boot_outer_weights
    cl_construct_scores <- cl_change$cl_boot_construct_scores
    cl_cross_loadings <- cl_change$cl_boot_cross_loadings
    cl_index <- cl_change$cl_change

    # Create an object to be returned
    boot_object <- list(boot_data = resampled_data,
                        boot_path_coef = path_coef,
                        boot_outer_loadings = outer_loadings,
                        boot_outer_weights = outer_weights,
                        boot_rSquared = rSquared,
                        boot_construct_scores = scores,
                        boot_total_effects = boot_total_effects,
                        boot_HTMT = boot_HTMT,
                        boot_cross_loadings = cross_loadings,
                        sc_path_coef = sc_path_coef,
                        sc_outer_loadings = sc_outer_loadings,
                        sc_outer_weights = sc_outer_weights,
                        sc_cross_loadings = sc_cross_loadings,
                        sc_construct_scores = sc_construct_scores,
                        sc_index = sc_index,
                        cl_path_coef = cl_path_coef,
                        cl_outer_loadings = cl_outer_loadings,
                        cl_outer_weights = cl_outer_weights,
                        cl_cross_loadings = cl_cross_loadings,
                        cl_construct_scores = cl_construct_scores,
                        cl_index = cl_index)
    
    return(boot_object)
  }
  
  # Export necessary variables and functions to the cluster
  parallel::clusterExport(cl, varlist = c("data",
                                          "measurement_model",
                                          "structural_model",
                                          "inner_weights",
                                          "bootstrap_sample",
                                          "convert_to_table_output",
                                          "HTMT",
                                          "total_effects",
                                          "seminr_model",
                                          "cross_loadings",
                                          "boot_sign_change",
                                          "boot_construct_level_change"), envir = environment())
  
  # Perform bootstrap in parallel
  bootstrap_results <- parallel::parLapply(cl, 1:nboot, bootstrap_sample, data, measurement_model, structural_model, inner_weights, seminr_model, cross_loadings)
  
  # Stop the cluster
  parallel::stopCluster(cl)
  
  # Separate the raw data and model estimates
  seminr_model$boot_raw_data <- lapply(bootstrap_results, function(res) res$boot_data)
  seminr_model$boot_path_coef <- lapply(bootstrap_results, function(res) res$boot_path_coef)
  seminr_model$boot_outer_loadings <- lapply(bootstrap_results, function(res) res$boot_outer_loadings)
  seminr_model$boot_outer_weights <- lapply(bootstrap_results, function(res) res$boot_outer_weights)
  seminr_model$boot_rSquared <- lapply(bootstrap_results, function(res) res$boot_rSquared)
  seminr_model$boot_construct_scores <- lapply(bootstrap_results, function(res) res$boot_construct_scores)
  seminr_model$boot_HTMT <- lapply(bootstrap_results, function(res) res$boot_HTMT)
  seminr_model$boot_total_effects <- lapply(bootstrap_results, function(res) res$boot_HTMT)
  seminr_model$boot_cross_loadings <- lapply(bootstrap_results, function(res) res$boot_cross_loadings)
  seminr_model$sc_path_coef <- lapply(bootstrap_results, function(res) res$sc_path_coef)
  seminr_model$sc_outer_loadings <- lapply(bootstrap_results, function(res) res$sc_outer_loadings)
  seminr_model$sc_outer_weights <- lapply(bootstrap_results, function(res) res$sc_outer_weights)
  seminr_model$sc_cross_loadings <- lapply(bootstrap_results, function(res) res$sc_cross_loadings)
  seminr_model$sc_construct_scores <- lapply(bootstrap_results, function(res) res$sc_construct_scores)
  seminr_model$sc_index <- lapply(bootstrap_results, function(res) res$sc_index)
  seminr_model$cl_path_coef <- lapply(bootstrap_results, function(res) res$cl_path_coef)
  seminr_model$cl_outer_loadings <- lapply(bootstrap_results, function(res) res$cl_outer_loadings)
  seminr_model$cl_outer_weights <- lapply(bootstrap_results, function(res) res$cl_outer_weights)
  seminr_model$cl_cross_loadings <- lapply(bootstrap_results, function(res) res$cl_cross_loadings)
  seminr_model$cl_construct_scores <- lapply(bootstrap_results, function(res) res$cl_construct_scores)
  seminr_model$cl_index <- lapply(bootstrap_results, function(res) res$cl_index)
  
  
  class(seminr_model) <- c("boot_seminr_model_full")
  message("SEMinR Model successfully bootstrapped")
  
  return(seminr_model)
}
