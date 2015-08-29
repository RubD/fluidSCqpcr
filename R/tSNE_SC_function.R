
#' Function to create a t-SNE map from single cell data
#'
#' This function creates a heatmap?
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param perplexity measure for information about the Shannon entropy, see reference
#' @param clusterColumn name of column with cluster information
#' @param NAvalues how to process NA values, remove or replace with not_detected_value (default = 0)
#' @param not_detected_value value to replace NA values with
#' @return returns a t-SNE plot
#' @export
#' @details NA
#' @references http://lvdmaaten.github.io/tsne/
#' @examples
#' tSNE_SC()

tSNE_SC <- function(fluidSCproc,  based_on_values = "log2Ex", perplexity = 5, clusterColumn = NULL,
                    NAvalues = c("remove_assays","replace_with_not_detected"),not_detected_value = 0) {
  
  # load libraries
  library(ggplot2)
  library(Rtsne)
  
  # check input
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  stopifnot(class(fluidSCproc) == "fluidSCproc")
  normFluidCt <- fluidSCproc$data
  
  # without cluster column to indicate colors
  if(is.null(clusterColumn)) {
    
    normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values)
    
    
    ## for merged fluidSCproc objects, what to do with NA values?
    NAvalues <- match.arg(NAvalues)
    
    # OR remove genes with NA values
    if(NAvalues == "remove_assays") {
      NAgenes <- colnames(normMatrix)[apply(normMatrix, 2, FUN = function(x) any(is.na(x)))]
      if(identical(NAgenes, character(0))) {
        normMatrix <- normMatrix
      }
      else normMatrix <- normMatrix[,-(which(colnames(normMatrix) %in% NAgenes))]
    } 
    # OR replace NA values with 0
    else if(NAvalues == "replace_with_not_detected") {
      normMatrix[is.na(normMatrix)] <- not_detected_value
    } else stop("Choose one of the 2 options to take care of NA values")
    
    
    ## remove genes with zero variance!
    keepgenes <- !apply(normMatrix, MARGIN = 2, function(x) var(x,na.rm = T)) == 0
    keepgenes <- names(keepgenes)[keepgenes]
    normMatrix <- normMatrix[, colnames(normMatrix) %in% c("Samples",keepgenes)]
    
    ## apply tSNE method
    tsne_out <- Rtsne(as.matrix(normMatrix[,2:ncol(normMatrix)]), perplexity = perplexity)
    dfr <- as.data.frame(tsne_out$Y)
    
    pl <- ggplot(dfr, aes(V1, V2))
    pl <- pl + geom_point(size = 4) + labs(list(x = "Dimension 1", y = "Dimension 2"))
    pl <- pl + theme_bw() + theme(axis.title = element_text(size = 22))
    return(pl)
  } 
  
  # with cluster column to indicate colors
  else {
    
    myformula <- paste0("Samples + ",clusterColumn, " ~ Assays")
    normMatrix <- dcast(normFluidCt, formula = eval(parse(text = myformula)), value.var = based_on_values)
    
    
    ## for merged fluidSCproc objects, what to do with NA values?
    NAvalues <- match.arg(NAvalues)
    
    # OR remove genes with NA values
    if(NAvalues == "remove_assays") {
      NAgenes <- colnames(normMatrix)[apply(normMatrix, 2, FUN = function(x) any(is.na(x)))]
      if(identical(NAgenes, character(0))) {
        normMatrix <- normMatrix
      }
      else normMatrix <- normMatrix[,-(which(colnames(normMatrix) %in% NAgenes))]
    } 
    # OR replace NA values with 0
    else if(NAvalues == "replace_with_not_detected") {
      normMatrix[is.na(normMatrix)] <- not_detected_value
    } else stop("Choose one of the 2 options to take care of NA values")
    
    
    # remove genes with zero variance!
    keepgenes <- !apply(normMatrix, MARGIN = 2, function(x) var(x,na.rm = T)) == 0
    keepgenes <- names(keepgenes)[keepgenes]
    normMatrix <- normMatrix[, colnames(normMatrix) %in% c("Samples",clusterColumn,keepgenes)]
    
    # apply tSNE method
    tsne_out <- Rtsne(as.matrix(normMatrix[,3:ncol(normMatrix)]), perplexity = perplexity)
    dfr <- as.data.frame(tsne_out$Y)
    dfr <- cbind(dfr, cluster = factor(normMatrix[,clusterColumn]))
    
    pl <- ggplot(dfr, aes(V1, V2, color = cluster))
    pl <- pl + geom_point(size = 4) + labs(list(x = "Dimension 1", y = "Dimension 2"))
    pl <- pl + theme_bw() + theme(axis.title = element_text(size = 22))
    return(pl)
    
  }
  
  
}
