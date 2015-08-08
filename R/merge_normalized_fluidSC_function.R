
#' Function to merge normalized fluidSCproc objects
#'
#' This function merges multiple fluidSCproc S3 objects and creates a new cell-to-cell median normalization log2 expression value for combined samples
#' @param fluidSCproc_list list of fluidSCproc S3 objects, need to have the same column names
#' @return returns a fluidSCproc S3 object
#' @export
#' @details 
#' 
#' Normalization is based on cell-to-cell median. Offset value between median log2 expression in sample vs average of all median values.
#' End result is that all cells will have the same median log2 expression.
#' @references Livak et al. "Methods for qPCR gene expression profiling applied to 1440 lymphoblastoid single cells"
#' @examples
#' merge_normalized_fluidSC()




merge_normalized_fluidSC <- function(fluidSCproc_list) {
  
  # merge different fluidSCproc S3 objects with rbind
  datalist <- lapply(fluidSCproc_list, FUN = function(x) x$data)
  mergedFluidSCproc <- do.call("rbind", datalist)
  
  ## NORMALIZE based on cell-to-cell median ##
  CellMedian <- aggregate(mergedFluidSCproc[!mergedFluidSCproc$log2Ex == 0,]$log2Ex, by = list(mergedFluidSCproc[!mergedFluidSCproc$log2Ex == 0,]$Samples), FUN = function(x) median(x))
  AllCellsAverage <- mean(CellMedian$x)
  offsetValue <- AllCellsAverage - CellMedian$x
  CellMedian$offsetValue <- offsetValue
  
  log2ExNormMerge <- NULL
  for(i in 1:nrow(mergedFluidSCproc)) {
    
    row <- mergedFluidSCproc[i, ]
    
    if(row[,"log2Ex"] == 0 ) x <- 0
    else {
      cell <- as.character(row[,"Samples"])
      offset <- CellMedian[CellMedian$Group.1 == cell,]$offsetValue
      x <- row[,"log2Ex"] + offset
    }
    
    log2ExNormMerge <- c(log2ExNormMerge, x)
    
  }
  
  mergedFluidSCproc$log2ExNormMerge <- log2ExNormMerge
  mergedFluidSCproc$lin2ExNormMerge <- 2^mergedFluidSCproc$log2ExNormMerge
  
  return(fluidSCproc(mergedFluidSCproc, LoD = "mixed - merged"))
}
