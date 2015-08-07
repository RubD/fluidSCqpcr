
#' Function to identify bad samples
#' 
#' This function determines the quality of the samples based on an indication of missing data. 
#' The more an assay is expressed in all samples (expression probability), the heavier its weight.
#' @param fluidSCraw fluidSCraw S3 object
#' @param sizefactor factor to multiply median missing data value, defaults to 3 see details
#' @return returns character vector with names of bad Samples
#' @export
#' @references Livak et al. "Methods for qPCR gene expression profiling applied to 1440 lymphoblastoid single cells"
#' @details sizefactor 1st step) calculate expression probability per assay
#' 2nd step) failed calls get the value of this expression probability and are summed per sample
#' 3th step) if sum of sample i > median for all samples * sizefactor => consider it a bad sample
#' @examples
#' identifyBadAssays()

identifyBadSamples <- function(fluidSCraw, sizefactor = 3) {
  
  
  rawFluidCt <- fluidSCraw$data
  
  # counttable for pass and fail calls for each primer
  AssayCalls <- lapply(split(rawFluidCt$Call, as.character(rawFluidCt$Assays)), FUN = function(x) table(x))
  # count expr probability for each assay
  exprProb <- lapply(AssayCalls, FUN = function(x) x[2] / (x[1] + x[2]))
  rawFluidCt <- rawFluidCt[order(rawFluidCt$Samples, rawFluidCt$Assays), ]
  rawFluidCt$exprProb <- as.vector(unlist(exprProb))
  
  # give failed calls a weight (from the exprProb)
  rawFluidCt$failWeight <- ifelse(rawFluidCt$Call == "Fail", rawFluidCt$exprProb, 0)
  # sum the failweights per sample
  sampleFail <- lapply(split(rawFluidCt$failWeight, rawFluidCt$Samples), FUN = function(x) sum(x))
  
  # determine bad samples based on failweight relative to median missing data
  medianMissingData <- median(unlist(sampleFail))
  culSamples <- names(sampleFail[sampleFail > sizefactor * medianMissingData])
  
  return(culSamples)
  
}
