#' Constructor function for fluidSCproc S3 object from dataframe 
#' 
#' This function creates a fluidSCproc S3 object.
#' @param x dataframe with at least 5 columns: Samples, Assays, Value, Call and log2Ex
#' @param LoD need to provide Limit of Detection (LoD) that you used to calculate log2 expression values (log2Ex)
#' @keywords construct fluidSCproc
#' @return S3 object fluidSCproc for analysis and visualization
#' @export
#' @examples
#' fluidSCproc()
#' 
fluidSCproc <- function(x, LoD = NULL) {
  
  if(class(x) != "data.frame" | !all(c("Assays","Samples","Value","Call","log2Ex") %in% colnames(x))) stop("input needs to be dataframe with at least 5 columns: Samples, Assays, Value, Call and log2Ex")
  if(is.null(LoD)) stop("you need provide information about LoD used")
  
  uniqSamples <- as.character(unique(x[,"Samples"]))
  uniqAssays <- as.character(unique(x[,"Assays"]))
  procInfo <- list(LoD = LoD)
  
  structure(list(data = x, unique_samples = uniqSamples , unique_assays = uniqAssays, proc_info = procInfo), class ="fluidSCproc")
  
}