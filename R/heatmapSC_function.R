
#' Function to create a heatmap from single cell data
#'
#' This function creates a heatmap?
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param distance_method choose method to calculute distance matrix: pearson, spearman, euclidean, maxium, manhatten, binary
#' @param clust_method choose method to cluster distance matrix: average, complete, ward.D, ward.D2, single, median
#' @param NAvalues how to process NA values, remove or replace with not_detected_value (default = 0)
#' @param not_detected_value value to replace NA values with
#' @param hmcolor give your own colors or colorpalette
#' @return returns a heatmap plot
#' @export
#' @details NA
#' @examples
#' heatmapSC()


heatmapSC <- function(fluidSCproc,  based_on_values = "log2Ex", distance_method = c("pearson","spearman","euclidean","maximum","manhatten","binary"),
                      clust_method = c("average","complete","ward.D","ward.D2","single","median"),NAvalues = c("remove_assays","replace_with_not_detected"),not_detected_value = 0,
                      hmcolor =  NULL, ...) {
  
  
  library(gplots)
  library(reshape2)
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  clust_method <- match.arg(clust_method)
  distance_method <- match.arg(distance_method)
  if(is.null(hmcolor)) hmcolor <- colorRampPalette( c("#313695", "#4575b4", "#74add1",  "#ffffbf", "#f46d43", "#d73027","#a50026"))(48)
  
  c("#313695", "#4575b4", "#74add1",  "#ffffbf", "#f46d43", "#d73027","#a50026")
  
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  
  ## for merged fluidSCproc objects, what to do with NA values?
  NAvalues <- match.arg(NAvalues)
  
  # remove genes with NA values
  if(NAvalues == "remove_assays") {
    NAgenes <- colnames(normMatrix)[apply(normMatrix, 2, FUN = function(x) any(is.na(x)))]
    if(identical(NAgenes, character(0))) {
      normMatrix <- normMatrix
    }
    else normMatrix <- normMatrix[,-(which(colnames(normMatrix) %in% NAgenes))]
  } 
  # replace NA values with 0
  else if(NAvalues == "replace_with_not_detected") {
    normMatrix[is.na(normMatrix)] <- not_detected_value
  } else stop("Choose one of the 2 options to take care of NA values")
  
  
  
  
  
  # remove genes with variation (not informative)
  var0Genes <- colnames(normMatrix)[apply(normMatrix, 2, var) == 0]
  normMatrix <- normMatrix[,!colnames(normMatrix) %in% var0Genes]
  
  
  if(distance_method %in% c("pearson","spearman")) {
    
    my_heatmap <- heatmap.2(as.matrix(t(normMatrix)), hclustfun = function(x) hclust(x, method = clust_method),
                            distfun = function(x) as.dist(1 - cor(t(x),method = distance_method)),
                            trace = "none", scale = "row",col = hmcolor, ... )
    
    
  } else if(distance_method %in% c("euclidean","maximum","manhatten","binary")) {
    
    my_heatmap <- heatmap.2(as.matrix(t(normMatrix)),hclustfun = function(x) hclust(x, method = clust_method),
                            distfun = function(x) dist(x, method = distance_method),
                            trace = "none", scale = "row",col = hmcolor, ... )
    
  } else stop("distance method is not known")
  
}
