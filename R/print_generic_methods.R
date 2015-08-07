#' Genereric print method for S3 object fluidSCraw
#' 
#' This function prints an S3 fluidSCraw object
#' @param x fluidSCraw object
#' @keywords print fluidSCraw
#' @return prints first 6 lines of dataframe and additional information
#' @export
#' @examples
#' print.fluidSCraw(x)
print.fluidSCraw <- function(x) {
  
  print(head(x$data))
  cat("\n")
  cat("Number of unique assays/genes: ",length(x$unique_assays), "\n")
  cat("Number of unique samples/cells: ", length(x$unique_samples))
  
  
}


#' Genereric print method for S3 object fluidSCproc
#' 
#' This function prints an S3 fluidSCproc object
#' @param x fluidSCproc object
#' @keywords print fluidSproc
#' @return prints first 6 lines of dataframe and additional information
#' @export
#' @examples
#' print.fluidSCproc(x)
print.fluidSCproc <- function(x) {
  
  print(head(x$data))
  cat("\n")
  cat("Number of unique assays/genes: ",length(x$unique_assays), "\n")
  cat("Number of unique samples/cells: ", length(x$unique_samples), "\n")
  cat("Limit of detection used: ", x$proc_info$LoD)
  
}

