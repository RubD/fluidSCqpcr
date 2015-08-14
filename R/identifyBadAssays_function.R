
#' Function to identifie bad assays
#' 
#' This function determines if assays are of good quality based on the ratio of good versus bad quality calls
#' @param fluidSCraw fluidSCraw S3 object
#' @param min_ratio_goodVSbad ratio of number of good calls over bad calls, ratio = #good/#bad. Defaults to 1, i.e. minimum half of the calls must be good
#' @param min_good_datapoints minimum number of good quality datapoints needed
#' @return returns list with overview of valued and result, i.e. names of assays that do not meet the requirements
#' @export
#' @examples
#' identifyBadAssays()

identifyBadAssays <- function(fluidSCraw, min_ratio_goodVSbad = 1, min_good_datapoints = 0) {
  
  ## DESCRIPTION ##
  # identifies bad primers based on:
  #    - expressed, i.e. not 999
  #    _ ratio of good calls / bad calls below treshold
  
  ## USAGE ##
  # identifyBadPrimers(fluidSCraw, min_ratio_goodVSbad)
  
  ## ARGUMENTS ##
  # fluidSCraw = fluidSCraw Object
  # min_ratio_goodVSbad = ratio of passed vs failed calls per primer
  
  rawFluidCt <- fluidSCraw$data
  
  temp_results <- rawFluidCt[rawFluidCt$Value < 999,]
  badPrimersList <- lapply(split(temp_results$Call, temp_results$Assays), FUN = function(x) table(x))
  badPrimers <- lapply(badPrimersList, FUN = function(x) ifelse( (x[2] / x[1]) < min_ratio_goodVSbad, "bad", "good"))
  lowPrimers <- lapply(badPrimersList, FUN = function(x) ifelse(x[2] < min_good_datapoints, "low", "ok"))
  
  badGenes <- names(badPrimers)[badPrimers == "bad"]
  lowGenes <- names(lowPrimers)[lowPrimers == "low"]
  
  return(list(overview = badPrimersList, result_bad = badGenes, result_low = lowGenes))
  
}
