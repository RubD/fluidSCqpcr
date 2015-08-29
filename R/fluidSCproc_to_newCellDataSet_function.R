
#' Function to create a newCellDataSet object from fluidSCproc object to use in the monocle package
#'
#' This function creates newCellDatqSet object. This is the starting object for the monocle package.
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2ExNorm"
#' @param sampleData character vector with column names that contain information about the samples
#' @param featureData character vector with column names that contain information about the features (assays)
#' @param NAvalues how to process NA values, remove or replace with not_detected_value (default = 0)
#' @param not_detected_value value to replace NA values with
#' @return returns a t-SNE plot
#' @export
#' @details NA
#' @references http://monocle-bio.sourceforge.net/tutorial.html
#' @examples
#' fluidSCproc_to_newCellDataSet()

fluidSCproc_to_newCellDataSet <- function(fluidSCproc,based_on_values = "log2ExNormMerge",
                                          sampleData = NULL, featureData = NULL,
                                          NAvalues = c("remove_assays","replace_with_not_detected"), not_detected_value = 0, 
                                           ...) {
  
  ## load libraries
  library(monocle)
  
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  
  ## for merged fluidSCproc objects, what to do with NA values?
  NAvalues <- match.arg(NAvalues)
  
  # OR remove genes with NA values
  if(NAvalues == "remove_assays") {
    NAgenes <- colnames(normMatrix)[apply(normMatrix, 2, FUN = function(x) any(is.na(x)))]
    if(identical(NAgenes, character(0))) {
      normMatrix <- normMatrix
    }
    else {
      normMatrix <- normMatrix[,-(which(colnames(normMatrix) %in% NAgenes))]
      normFluidCt <- normFluidCt[!normFluidCt$Assays %in% NAgenes,]
    }
  } 
  # replace NA values with 0
  else if(NAvalues == "replace_with_not_detected") {
    normMatrix[is.na(normMatrix)] <- not_detected_value
  } else stop("Choose one of the 2 options to take care of NA values")
  
  
  # 3. data sets to make expressionSet (newCellDataSet)
  expression_matrix <- t(as.matrix(normMatrix))
  
  
  # phenotype sample data
  pdataFrame <- unique(normFluidCt[,c("Samples",sampleData), drop = F])
  rownames(pdataFrame) <- pdataFrame$Samples
  pdataFrame <- pdataFrame[match(colnames(expression_matrix), pdataFrame$Samples), ,drop = F]
  pdataFrame <- new("AnnotatedDataFrame", data = pdataFrame)
  
  # feature assay data
  fdataFrame <- unique(normFluidCt[,c("Assays", featureData), drop = F])
  rownames(fdataFrame) <- fdataFrame$Assays
  fdataFrame <- fdataFrame[match(rownames(expression_matrix), fdataFrame$Assays), ,drop = F]
  fdataFrame <- new("AnnotatedDataFrame", data = fdataFrame)
  
  #return(list(expression_matrix=expression_matrix, pdataFrame=pdataFrame, fdataFrame=fdataFrame ))
  return(newCellDataSet(cellData = expression_matrix, phenoData = pdataFrame, featureData = fdataFrame, ...))
  
}

