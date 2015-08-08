
#' Function to create a heatmap from single cell data
#'
#' This function creates a heatmap?
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param distance_method choose method to calculute distance matrix: pearson, spearman, euclidean, maxium, manhatten, binary
#' @param clust_method choose method to cluster distance matrix: average, complete, ward.D, ward.D2, single, median
#' @param hmcolor give your own colors or colorpalette
#' @return returns a heatmap plot
#' @export
#' @details NA
#' @examples
#' heatmapSC()


heatmapSC <- function(fluidSCproc,  based_on_values = "log2Ex", distance_method = c("pearson","spearman","euclidean","maximum","manhatten","binary"),
                      clust_method = c("average","complete","ward.D","ward.D2","single","median"), hmcolor =  NULL, ...) {
  
  
  library(gplots)
  
  clust_method <- match.arg(clust_method)
  distance_method <- match.arg(distance_method)
  if(is.null(hmcolor)) hmcolor <- colorRampPalette( c("#313695", "#4575b4", "#74add1",  "#ffffbf", "#f46d43", "#d73027","#a50026"))(48)
  
  c("#313695", "#4575b4", "#74add1",  "#ffffbf", "#f46d43", "#d73027","#a50026")
  
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  # remove genes with variation (not informative)
  var0Genes <- colnames(normMatrix)[apply(normMatrix, 2, var) == 0]
  normMatrix <- normMatrix[,!colnames(normMatrix) %in% var0Genes]
  
  
  if(distance_method %in% c("pearson","spearman")) {
    
    my_heatmap <- heatmap.2(as.matrix(t(normMatrix)), distfun = function(x) as.dist(1 - cor(t(x),method = distance_method)), trace = "none", scale = "row",col = hmcolor, ... )
    
    
  } else if(distance_method %in% c("euclidean","maximum","manhatten","binary")) {
    
    my_heatmap <- heatmap.2(as.matrix(t(normMatrix)), distfun = function(x) dist(x, method = distance_method), trace = "none", scale = "row",col = hmcolor, ... )
    
  } else stop("distance method is not known")
  
}
