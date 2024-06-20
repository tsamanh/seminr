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
total_effects <- function(path_coef) {
  output <- path_coef
  paths <- path_coef
  while (sum(paths) > 0) {
    paths <- paths %*% path_coef
    output <- output + paths
  }
  return(output)
}

convert_to_table_output <- function(matrix) {
  class(matrix) <- append(class(matrix), "table_output")
  return(matrix)
}


HTMT <- function(seminr_model) {
  if (is.null(seminr_model$hoc)) {
    constructs <- intersect(unique(seminr_model$smMatrix),unique(seminr_model$mmMatrix[,1 ]))
  } else {
    constructs <- intersect(unique(c(seminr_model$smMatrix, seminr_model$first_stage_model$smMatrix)),unique(seminr_model$mmMatrix[,1 ]))
  }
  
  HTMT <- matrix(, nrow=length(constructs), ncol=length(constructs),
                 dimnames = list(constructs,constructs))
  for (constructi in constructs[1:(length(constructs)-1)]) {
    for (constructj in constructs[(which(constructs == constructi)+1):length(constructs)]) {
      manifesti <- seminr_model$mmMatrix[seminr_model$mmMatrix[, 1] == constructi, "measurement"]
      manifestj <- seminr_model$mmMatrix[seminr_model$mmMatrix[, 1] == constructj, "measurement"]
      item_correlation_matrix <- abs(stats::cor(seminr_model$data[, manifesti],seminr_model$data[, manifestj]))
      HTHM <- mean(item_correlation_matrix)
      if(length(manifesti)>1 ) {
        cor_matrix <- abs(stats::cor(seminr_model$data[, manifesti], seminr_model$data[, manifesti]))
        diag(cor_matrix) <- 0
        MTHM <- (2/(length(manifesti)*(length(manifesti)-1)))*(sum(cor_matrix[!lower.tri(cor_matrix)]))
      } else {
        MTHM <- 1
      }
      if(length(manifestj)>1) {
        cor_matrix2 <- abs(stats::cor(seminr_model$data[, manifestj], seminr_model$data[, manifestj]))
        diag(cor_matrix2) <- 0
        MTHM <- sqrt(MTHM * (2/(length(manifestj)*(length(manifestj)-1)))*(sum(cor_matrix2[!lower.tri(cor_matrix2)])))
      } else {
        MTHM <- sqrt(1 * MTHM)
      }
      HTMT[constructi, constructj] <- HTHM / MTHM
    }
  }
  convert_to_table_output(HTMT)
}

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
  
  # Check for and create random seed if NULL
  if (is.null(seed)) {seed <- sample.int(100000, size = 1)}
  
  # Function to perform a single bootstrap resample and run seminr::estimate_pls
  bootstrap_sample <- function(i, data, measurement_model, structural_model, inner_weights) {
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
