#' Constructor function for fluidSCraw S3 object from dataframe 
#' 
#' This function creates a fluidSCraw S3 object.
#' @param x dataframe with at least 4 columns: Samples, Assays, Value and Call
#' @keywords construct fluidSCraw
#' @return S3 object fluidSCraw for pre-processing steps
#' @export
#' @examples
#' fluidSCraw()

fluidSCraw <- function(x) {
  
  if(class(x) != "data.frame" | !all(c("Assays","Samples","Value","Call") %in% colnames(x))) stop("input needs to be dataframe with at least 4 columns: Samples, Assays, Value and Call")
  
  uniqSamples <- as.character(unique(x[,"Samples"]))
  uniqAssays <- as.character(unique(x[,"Assays"]))
  
  structure(list(data = x, unique_samples = uniqSamples , unique_assays = uniqAssays), class ="fluidSCraw")
  
}
