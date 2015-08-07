
#' Function to call quality of values based on list with melting temparatures for all assays
#' 
#' This function takes as input a list with melt temperatures for all assays and calculates whether data points are within range 
#' @param fluidSCraw fluidSCraw S3 object
#' @param meltTempRange range around melt temperature peak 
#' @param named_listMeltTemp list with melting temperatures and names that correspond to the assays
#' @keywords melt temperature fluidSCraw
#' @return fluidSCraw S3 object were new Call column has been added or replaced
#' @export
#' @examples
#' decide_Call_meltTemp()

decide_Call_meltTemp <- function(fluidSCraw, meltTempRange = c(2,2), named_listMeltTemp = NULL) {
  
  ## DESCRIPTION ##
  # Recreates the quality calls based on self provided melt temperatures for all assays
  
  ## USAGE ##
  # decide_Call_meltTemp(fluidSCraw, meltTempRange = c(2,2), named_listMeltTemp = self_provided_list)
  
  ## ARGUMENTS ##
  # fluidSCraw = fluidSCraw Object
  # meltTempRange = lower and upper limit range around 'true' melt temperature
  # names_listMeltTemp = list with temperatures and names for the assays
  
  
  rawFluidCt <- fluidSCraw$data
  
  if(!all(unique(rawFluidCt$Assays) %in% names(named_listMeltTemp))) stop("the named list with melt temperatures must contain all the unique assays from the fluidSCraw object")
  
  newCall <- NULL
  for(i in 1:nrow(rawFluidCt)) {
    
    
    gene <- rawFluidCt[i ,"Assays"]
    melttemp <- rawFluidCt[i, "Melt_temp"]
    trueMelttemp <- as.vector(named_listMeltTemp[names(named_listMeltTemp) == gene])
    
    # determine whether melt temperature is around 'true' melt temperature
    lowerLimit <- trueMelttemp - meltTempRange[1]
    upperLimit <- trueMelttemp + meltTempRange[2]
    x <- ifelse(melttemp > lowerLimit & melttemp < upperLimit, "Pass","Fail")
    newCall <- c(newCall, x)
    
  }
  
  rawFluidCt$Call <- newCall
  
  return(fluidSCraw(rawFluidCt))
}
