
#' Function to normalize fluidSCraw object
#'
#' This function processes and normalizes a fluidSCraw S3 object. It also calculates a log2 expression value based on a Limit of Detection (LoD) score
#' @param fluidSCraw fluidSCraw S3 object
#' @param LoD Limit of Detection value, defaults to 26
#' @return returns a fluidSCproc S3 object
#' @export
#' @details 
#' - if value = 999 and call = Fail, replace value with LoD (not detected)
#' - if value > LoD and call = Pass, replace value with LoD (not detected)
#' - if value < LoD and call = Fail, replace value with mean expression of assay in all cells (missing data)
#' 
#' Normalization is based on cell-to-cell median. Offset value between median log2 expression in sample vs average of all median values.
#' End result is that all cells will have the same median log2 expression.
#' @references Livak et al. "Methods for qPCR gene expression profiling applied to 1440 lymphoblastoid single cells"
#' @examples
#' normalize_fluidSC()

normalize_fluidSC <- function(fluidSCraw, LoD = 26) {
  
  
  rawFluidCt <- fluidSCraw$data
  
  ## ADD log2exp scores based on LoD ##
  # mean value for gene / sample group
  temp_table <- rawFluidCt[rawFluidCt$Value < 999 & rawFluidCt$Call == "Pass", ]
  meanGeneExpression <- aggregate(temp_table$Value, by = list(temp_table$Assays), FUN = function(x) mean(x))
  
  
  newVal <- NULL
  for(i in 1:nrow(rawFluidCt)) {
    
    row <- rawFluidCt[i, ]
    
    if(row[,"Value"] < LoD & row[,"Call"] == "Pass") x <- row[,"Value"]
    else if(row[,"Value"] == 999 & row[,"Call"] == "Fail") x <- LoD
    else if(row[,"Value"] > LoD & row[,"Call"] == "Pass") x <- LoD
    else if(row[,"Value"] < LoD & row[,"Call"] == "Fail") {
      gene <- as.character(row[,"Assays"])
      x <- meanGeneExpression[meanGeneExpression$Group.1 == gene,]$x
    }
    
    newVal <- c(newVal, x)
    
  }
  
  rawFluidCt$newValue <- newVal
  rawFluidCt$log2Ex <- LoD - newVal
  rawFluidCt$lin2Ex <- 2^rawFluidCt$log2Ex
  
  
  
  ## NORMALIZE based on cell-to-cell median ##
  CellMedian <- aggregate(rawFluidCt[!rawFluidCt$log2Ex == 0,]$log2Ex, by = list(rawFluidCt[!rawFluidCt$log2Ex == 0,]$Samples), FUN = function(x) median(x))
  AllCellsAverage <- mean(CellMedian$x)
  offsetValue <- AllCellsAverage - CellMedian$x
  CellMedian$offsetValue <- offsetValue
  
  
  log2ExNorm <- NULL
  for(i in 1:nrow(rawFluidCt)) {
    
    row <- rawFluidCt[i, ]
    
    if(row[,"log2Ex"] == 0 ) x <- 0
    else {
      cell <- as.character(row[,"Samples"])
      offset <- CellMedian[CellMedian$Group.1 == cell,]$offsetValue
      x <- row[,"log2Ex"] + offset
    }
    
    log2ExNorm <- c(log2ExNorm, x)
    
  }
  
  rawFluidCt$log2ExNorm <- log2ExNorm
  rawFluidCt$lin2ExNorm <- 2^rawFluidCt$log2ExNorm
  
  
  
  return(fluidSCproc(rawFluidCt, LoD = LoD))
  
}
