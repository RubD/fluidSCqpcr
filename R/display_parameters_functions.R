#' Function to calculate parameters of distribution
#' 
#' This function calculates distribution parameters like range, quantile, median, ...
#' @param rawFluidCt data slot from S3 fluidSCraw or fluidSCproc object
#' @param reference_genes reference genes in your assays
#' @param groupColumn create facets based on column with grouping 
#' @return outputs text in console with information
#' @export
#' @examples
#' display_parameters_function()

display_parameters_function <- function(rawFluidCt, reference_genes = NULL, groupColumn = NULL) {
  
  if(!is.null(groupColumn)) cat("For group: ", unique(rawFluidCt[ ,groupColumn]))
  cat("\n")
  all_assays_median <- median(rawFluidCt[rawFluidCt$Value < 999 & !is.na(rawFluidCt$Value), "Value"])
  cat("median Ct value for all assays = ",all_assays_median)
  cat("\n")
  
  if(!is.null(reference_genes)) {
    
    median_vec <- NULL
    
    for(gene in reference_genes) {
      
      assay_median <- median(rawFluidCt[rawFluidCt$Value < 999 & !is.na(rawFluidCt$Value) & rawFluidCt$Assays == gene, "Value"])
      cat("median Ct value for ", gene," = ",assay_median)
      cat("\n")
      
      median_vec <- c(median_vec, assay_median)
    }
    
    median_vec <- c(median_vec, all_assays_median)
    
    geom_mean <- prod(median_vec)^(1/length(median_vec))
    aritm_mean <- sum(median_vec)/length(median_vec)
    
    range_res <- range(rawFluidCt[rawFluidCt$Value < 999 & !is.na(rawFluidCt$Value), "Value"])
    quantile_res <- quantile(rawFluidCt[rawFluidCt$Value < 999 & !is.na(rawFluidCt$Value), "Value"])
    
    cat("Geometric mean for all assays and reference genes (",reference_genes,") =", geom_mean)
    cat("\n")
    cat("Arithmic mean for all assays and reference genes (",reference_genes,") =", aritm_mean)
    cat("\n")
    cat("Range for all Ct values = ",range_res)
    cat("\n")
    cat("Quantile results for all Ct values = ",quantile_res)
    cat("\n")
  }
  
  else {
    
    range_res <- range(rawFluidCt[rawFluidCt$Value < 999 & !is.na(rawFluidCt$Value), "Value"])
    quantile_res <- quantile(rawFluidCt[rawFluidCt$Value < 999 & !is.na(rawFluidCt$Value), "Value"])
    
    cat("Range for all Ct values = ",range_res)
    cat("\n")
    cat("Quantile results for all Ct values = ",quantile_res)
    cat("\n")
    
    
  }
  
  cat("\n")
  cat("\n")
  
}

#' Function to output parameters of distribution
#' 
#' This function outputs distribution parameters like range, quantile, median, ... calls display_parameters_function()
#' @param rawFluidCt data slot from S3 fluidSCraw or fluidSCproc object
#' @param reference_genes reference genes in your assays
#' @param groupColumn create facets based on column with grouping 
#' @return outputs text in console with information
#' @export
#' @examples
#' display_parameters_distribution()


# function to use display_paramters_function on 1 or more conditions
display_parameters_distribution <- function(fluidSCraw, reference_genes = NULL, groupColumn = NULL) {
  
  
  rawFluidCt <- fluidSCraw$data
  
  if(!is.null(groupColumn)) {
    
    tempsplit <- split(rawFluidCt, list(rawFluidCt[ ,groupColumn]))
    
    lapply(tempsplit, FUN = function(x) display_parameters_function(x, reference_genes = reference_genes, groupColumn = groupColumn))
    
    
  } 
  else {
    
    display_parameters_function(rawFluidCt, reference_genes = reference_genes, groupColumn = groupColumn)
    
  }
  
}
