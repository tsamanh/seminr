#' function to flatten the nested list

flatten_list <- function(boot_list) {
  flat <- array(unlist(boot_list), dim = c(nrow(boot_list[[1]]),
                                           ncol(boot_list[[1]]),
                                           length(boot_list)))
  dimnames(flat) <- list(rownames(boot_list[[1]]), colnames(boot_list[[1]]), NULL)
  return(flat)
}

#' Summary function for the bootstrap full model
#' @export


summary.boot_seminr_model_full <- function (boot_object, alpha = 0.05, ...) {
  stopifnot(inherits(boot_object, "boot_seminr_model_full"))
  
  # Convert the bootstrapped list into 3D array
  flat_path <- flatten_list(boot_object$boot_path_coef)
  flat_weights <- flatten_list(boot_object$boot_outer_weights)
  flat_loadings <- flatten_list(boot_object$boot_outer_loadings)
  flat_HTMT <- flatten_list(boot_object$boot_HTMT)
  flat_cross_loadings <- flatten_list(boot_object$boot_cross_loadings)
  
  # Convert the sc boostrapped list into 3D array
  flat_sc_path <- flatten_list(boot_object$sc_path_coef)
  flat_sc_weights <- flatten_list(boot_object$sc_outer_weights)
  flat_sc_loadings <- flatten_list(boot_object$sc_outer_loadings)
  flat_sc_cross_loadings <- flatten_list(boot_object$sc_cross_loadings)
  
  # Convert the cl bootstrapped list into 3D array
  flat_cl_path <- flatten_list(boot_object$cl_path_coef)
  flat_cl_weights <- flatten_list(boot_object$cl_outer_weights)
  flat_cl_loadings <- flatten_list(boot_object$cl_outer_loadings)
  flat_cl_cross_loadings <- flatten_list(boot_object$cl_cross_loadings)
  
  # Convert the di bootstrapped list into 3D array
  flat_di_path <- flatten_list(boot_object$di_path_coef)
  flat_di_weights <- flatten_list(boot_object$di_outer_weights)
  flat_di_loadings <- flatten_list(boot_object$di_outer_loadings)
  flat_di_cross_loadings <- flatten_list(boot_object$di_cross_loadings)
  
  
  # Get the Mean, SD, and Confidence interval of the arrays
  path_summary <- parse_boot_array(boot_object$path_coef, flat_path, alpha = alpha)
  weights_summary <- parse_boot_array(boot_object$outer_weights, flat_weights, alpha = alpha)
  loadings_summary <- parse_boot_array(boot_object$outer_loadings, flat_loadings, alpha = alpha)
  HTMT_summary <- parse_boot_array(HTMT(boot_object), flat_HTMT, alpha = alpha)
  cross_loadings_summary <- parse_boot_array(cor(boot_object$rawdata, boot_object$construct_scores), flat_cross_loadings, alpha = alpha)
  
  sc_path_summary <- parse_boot_array(boot_object$path_coef, flat_sc_path, alpha = alpha)
  sc_weights_summary <- parse_boot_array(boot_object$outer_weights, flat_sc_weights, alpha = alpha) 
  sc_loadings_summary <- parse_boot_array(boot_object$outer_loadings, flat_sc_loadings, alpha = alpha)
  sc_cross_loadings_summary <- parse_boot_array(cor(boot_object$rawdata, boot_object$construct_scores), flat_sc_cross_loadings, alpha = alpha)
  
  cl_path_summary <- parse_boot_array(boot_object$path_coef, flat_cl_path, alpha = alpha)
  cl_weights_summary <- parse_boot_array(boot_object$outer_weights, flat_cl_weights, alpha = alpha)
  cl_loadings_summary <- parse_boot_array(boot_object$outer_loadings, flat_cl_loadings, alpha = alpha)
  cl_cross_loadings_summary <- parse_boot_array(cor(boot_object$rawdata, boot_object$construct_scores), flat_cl_cross_loadings, alpha = alpha)
  
  di_path_summary <- parse_boot_array(boot_object$path_coef, flat_di_path, alpha = alpha)
  di_weights_summary <- parse_boot_array(boot_object$outer_weights, flat_di_weights, alpha = alpha)
  di_loadings_summary <- parse_boot_array(boot_object$outer_loadings, flat_di_loadings, alpha = alpha)
  di_cross_loadings_summary <- parse_boot_array(cor(boot_object$rawdata, boot_object$construct_scores), flat_di_cross_loadings, alpha = alpha)
  
  boot_summary_object <- list(nboot = boot_object$boots,
                              boot_path_summary = path_summary,
                              boot_weights_summary = weights_summary,
                              boot_loadings_summary = loadings_summary,
                              boot_HTMT_summary = HTMT_summary,
                              boot_cross_loadings_summary = cross_loadings_summary,
                              sc_path_summary = sc_path_summary,
                              sc_weights_summary = sc_weights_summary,
                              sc_loadings_summary = sc_loadings_summary,
                              sc_cross_loadings_summary = sc_cross_loadings_summary,
                              cl_path_summary = cl_path_summary,
                              cl_weights_summary = cl_weights_summary,
                              cl_loadings_summary = cl_loadings_summary,
                              cl_cross_loadings_summary = cl_cross_loadings_summary,
                              di_path_summary = di_path_summary,
                              di_weights_summary = di_weights_summary,
                              di_loadings_summary = di_loadings_summary,
                              di_cross_loadings_summary = di_cross_loadings_summary)
  
  return(boot_summary_object)
}

