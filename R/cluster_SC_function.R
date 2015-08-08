

#' Function to cluster single cell data
#'
#' This function merges several clustering methods and appends the result to the fluidSCproc object
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param clustMethod choose method for clustering samples
#' @param nrClust number of clusters you want to detect with Kmeans
#' @param distMethod choose method to create distance matrix
#' @param corrMethod choose correlation score algorithm to use if clustMethod = Correlation
#' @param hclustMethod choose method to cluster distance matrix
#' @param selected_assays make selection of assays if wanted
#' @param cluster_column_name name of column that will be created with the clustering results
#' @return returns a fluidSCproc S3 object appended with the cluster results
#' @export
#' @details NA 
#' @examples
#' cluster_SC()


cluster_SC <- function(fluidSCproc,  based_on_values = "log2Ex", clustMethod = c("Kmeans","Hierarchical", "Correlation"), nrClust = 2,
                       distMethod = "euclidean", corrMethod = "pearson", hclustMethod = "average", selected_assays = NULL, cluster_column_name = "clust") {
  
  usedLoD <- fluidSCproc$proc_info$LoD
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  
  scores <- normMatrix
  if(!is.null(selected_assays)) scores <- scores[, selected_assays] # create subset from assays
  
  
  ## CLUSTERING ##
  if (clustMethod == "Kmeans") {
    K <- kmeans(scores, centers = nrClust, nstart = 25, iter.max = 1000, algorithm = "Hartigan-Wong")
    I <- K$cluster
  }
  
  else if (clustMethod == "Hierarchical") {
    H <- hclust(dist(scores, method = distMethod), method = hclustMethod)
    I <- cutree(H, k = nrClust)
  }
  
  else if (clustMethod == "Correlation") {
    distfun <- function(x) as.dist(1 - cor(t(x), method = corrMethod))
    mycordist <- distfun(scores)
    hclustfun <- function(x) hclust(x, method = hclustMethod)
    Hcor <- hclustfun(mycordist)
    I <- cutree(Hcor, k = nrClust)
  }
  
  
  ## merge cluster data with original data ##
  prepDfr <- as.data.frame(I); prepDfr$Samples <- rownames(prepDfr); colnames(prepDfr)[[1]] <- cluster_column_name
  if(!is.null(selected_assays)) normFluidCt <- normFluidCt[normFluidCt$Assays %in% selected_assays, ]
  mergeDfr <- merge(normFluidCt, prepDfr, by = "Samples" )
  
  
  return(fluidSCproc(mergeDfr, usedLoD))
  
  
}
