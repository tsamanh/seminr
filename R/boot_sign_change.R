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


