#' seminr bootstrap sign_change function
#' 
#' For sign correction during bootstrapping. Keep the sign of the original algorithm for each bootstrap
#' 
 
boot_sign_change <- function (boot_table, original_table) {
  if (any(sign(boot_table) != sign(original_table))) {
    index_sc <- which(sign(boot_table) != sign(original_table))
    sc_boot_table <- boot_table
    sc_boot_table[index_sc] <- sc_boot_table[index_sc] * (-1)
    return(sc_boot_table)
  } else {
    return(boot_table)
  }
}


boot_construct_level_change <- function(boot_model, original_model) {
  
  # Calculate the sum of boot weights and sum of original weights, then calculate the sum of sums and the difference of the sums
  sum_boot_weights <- apply(boot_model$outer_weights, 2, sum)
  sum_original_weights <- apply(original_model$outer_weights, 2, sum)
  sum_of_sum <- abs(sum_boot_weights + sum_original_weights)
  diff_of_sum <- abs(sum_boot_weights - sum_original_weights)
  flag <- sum_of_sum < diff_of_sum
  
  # Create the construct-level sign change objects
  cl_boot_outer_weights <- boot_model$outer_weights
  cl_boot_outer_loadings <- boot_model$outer_loadings
  cl_boot_construct_scores <- boot_model$construct_scores
  cl_boot_path_coef <- boot_model$path_coef
  
  # Flip the sign of the weights, loadings, construct scores if its sum of sums < difference of sums
  for (i in 1:length(flag))  {
    if (flag[i]) {
      cl_boot_outer_weights[,i] <- cl_boot_outer_weights[,i] * (-1)
      cl_boot_outer_loadings[,i] <- cl_boot__outer_loadings[,i] * (-1)
      cl_boot_construct_scores[,i] <- cl_boot_construct_scores[,i] * (-1)
    }
  }
  
  # Flip the sign of the path if one of the variable in an IV-DV pair flip sign
  for (i in 1:nrow(boot_model$structural_model)) {
    IV <- boot_model$structural_model[i,]["source"]
    DV <- boot_model$structural_model[i,]["target"]
    
    if ((flag[[IV]] - flag[[DV]])) {
      cl_boot_path_coef[IV,DV] <- cl_boot_path_coef[IV,DV] * (-1)
    }
  }
  
  # Calculate the new cross loadings
  cl_boot_cross_loadings <- cor(boot_model$rawdata, cl_boot_construct_scores)
  
  # change = TRUE if there is any flip
  cl_change <- any(flag)

  # Prepare the object for export
  cl_boot_object <- list(cl_boot_outer_weights = cl_boot_outer_weights,
                         cl_boot_outer_loadings = cl_boot_outer_loadings,
                         cl_boot_construct_scores = cl_boot_construct_scores,
                         cl_boot_cross_loadings = cl_boot_cross_loadings,
                         cl_boot_path_coef = cl_boot_path_coef,
                         cl_change = cl_change)
  
  return(cl_boot_object)
}
