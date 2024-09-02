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
      cl_boot_outer_loadings[,i] <- cl_boot_outer_loadings[,i] * (-1)
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
  
  # mark the constructs that flip
  cl_flip_constructs <- names(which(flag))

  # Prepare the object for export
  cl_boot_object <- list(cl_boot_outer_weights = cl_boot_outer_weights,
                         cl_boot_outer_loadings = cl_boot_outer_loadings,
                         cl_boot_construct_scores = cl_boot_construct_scores,
                         cl_boot_cross_loadings = cl_boot_cross_loadings,
                         cl_boot_path_coef = cl_boot_path_coef,
                         cl_change = cl_change,
                         cl_flip_constructs = cl_flip_constructs)
  
  return(cl_boot_object)
}


boot_dominant_indicator <- function(boot_model, original_model) {
  
  # Find the highest loading item of each construct
  construct_list <- original_model$constructs
  highest_loadings_list <- c()

  for (construct in construct_list) {
    loadings <- original_model$outer_loadings
    highest_loading <- which.max(abs(loadings[,construct]))
    highest_loadings_list <- c(highest_loadings_list, highest_loading)
  }
  
  names(highest_loadings_list) <- construct_list
  
  highest_loadings_list
  # Create the dominant indicator sign change objects
  di_outer_weights <- boot_model$outer_weights
  di_outer_loadings <- boot_model$outer_loadings
  di_construct_scores <- boot_model$construct_scores
  di_path_coef <- boot_model$path_coef
  
  # Create a flag vector to check if a construct score flip
  flag <- rep(FALSE, length(construct_list))
  names(flag) <- construct_list
  
  boot_model$outer_loadings[,construct] <- boot_model$outer_loadings[,construct] * (-1)
  
  for (construct in construct_list) {
    highest_loadings_var_index <- highest_loadings_list[construct]
    
    # Compare the sign of the loadings and flip if the sign is different from the original
    original_loadings <- original_model$outer_loadings[,construct]
    boot_loadings <- boot_model$outer_loadings[,construct]
    
    if (sign(original_loadings[highest_loadings_var_index]) != sign(boot_loadings[highest_loadings_var_index])) {
      di_outer_loadings[,construct] <- di_outer_loadings[,construct] * (-1)
      di_outer_weights[,construct] <- di_outer_weights[,construct] * (-1)
      di_construct_scores[,construct] <- di_construct_scores[,construct] * (-1)
      flag[construct] <- TRUE
    }
  }   

    # Flip the path coef if one of the IV-DV flag is TRUE
  for (i in 1:nrow(boot_model$structural_model)) {
    IV <- boot_model$structural_model[i,]["source"]
    DV <- boot_model$structural_model[i,]["target"]

    if ((flag[IV] - flag[DV])) {
      di_path_coef[IV,DV] <- di_path_coef[IV,DV] * (-1)
    }
  }
  
  # Recalculate the crossloadings
  di_cross_loadings <- cor(boot_model$rawdata, di_construct_scores)
  
  # Mark this bootstrap as flipped
  di_change <- any(flag)
  
  # Mark the constructs that flip
  di_flip_constructs <- names(which(flag))
  
  # Create an object to be returned
  di_object <- list(di_outer_weights = di_outer_weights,
                 di_outer_loadings = di_outer_loadings,
                 di_construct_scores = di_construct_scores,
                 di_cross_loadings = di_cross_loadings,
                 di_path_coef = di_path_coef,
                 di_change = di_change,
                 di_flip_constructs = di_flip_constructs)

  return(di_object)
}

boot_construct_score_correction <- function(boot_model, original_model) {
  boot_item <- boot_model$data
  boot_construct_scores <- boot_model$construct_scores
  csc_outer_weights <- boot_model$outer_weights
  csc_outer_loadings <- boot_model$outer_loadings
  
  # Initiate weights matrix to store the original weights
  original_weights <- matrix(0, nrow = ncol(boot_item), ncol = length(original_model$constructs))
  colnames(original_weights) <- colnames(original_model$outer_weights)
  rownames(original_weights) <- colnames(boot_item)
  
  # Add the original weights
  for (i in 1:ncol(original_weights)) {
    index <- which(original_model$outer_weights[,i] != 0)
    original_weights[index,i] <- original_model$outer_weights[index,i]
  }
  
  # Calculate new scores using original weights
  new_scores <- scale(boot_item) %*% original_weights
  
  # Get the sum and differences
  sum_scores <- new_scores + boot_construct_scores
  diff_scores <- new_scores - boot_construct_scores
  
  sum_sum <- apply(abs(sum_scores), 2, sum)
  sum_diff <- apply(abs(diff_scores), 2, sum)
  
  flip_constructs <- which(sum_sum < sum_diff)
  
  csc_boot_scores <- boot_construct_scores
  for (i in 1:ncol(csc_boot_scores)) {
    if (i %in% flip_constructs) {
      csc_boot_scores[,i] <- csc_boot_scores[,i] * -1
      csc_outer_weights[,i] <- csc_outer_weights[,i] * (-1)
      csc_outer_loadings[,i] <- csc_outer_loadings[,i] * (-1)
    }
  }

  csc_path_coef <- boot_model$path_coef
  csc_rSquare <- boot_model$rSquared
  
  dependent_var <- unique(boot_model$smMatrix[,"target"])
  for (i in 1:length(dependent_var)) {
    
    # get the name of independent variables
    independent_var_index <- which(boot_model$smMatrix[,"target"] == dependent_var[i])
    ind_var <- boot_model$smMatrix[,"source"][independent_var_index]
    
    # construct the formula for regression
    formula <- paste(dependent_var[i], "~", paste(ind_var, collapse = " + "))
    
    # run regression
    reg <- lm(formula, data = as.data.frame(csc_boot_scores))
    sum_reg <- summary(reg)
    
    # Extract the path
    csc_path_coef[ind_var, dependent_var[i]] <- reg$coefficients[-1]
    
    # Extract the Rsquare and adjusted Rsquare
    csc_rSquare[,dependent_var[i]] <- c(sum_reg$r.squared, sum_reg$adj.r.squared)
  }
  

  # Extract the cross loadings
  csc_cross_loadings <- cor(boot_item, csc_boot_scores)
  
  # Flag if there is a flip
  if (length(flip_constructs) != 0) {
    flag <- TRUE
  } else {
    flag <- FALSE
  }
  
  flag_constructs <- names(flip_constructs)
  
  # Return the object
  csc_object <- list(csc_construct_scores = csc_boot_scores,
                     csc_outer_loadings = csc_outer_loadings,
                     csc_outer_weights = csc_outer_weights,
                     csc_path_coef = csc_path_coef,
                     csc_rSquare = csc_rSquare,
                     csc_flag = flag,
                     csc_flip_constructs = flag_constructs)  
  
  return(csc_object)
}
