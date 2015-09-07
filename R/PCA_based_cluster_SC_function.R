#' Function to cluster single cell data
#'
#' This function merges several clustering methods and appends the result to the fluidSCproc object
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param scaleData boolean to scale data before clustering, default = F
#' @param clustMethod choose method for clustering samples
#' @param firstPC/secondPC principal components to use
#' @param nrClust number of clusters you want to detect with Kmeans
#' @param hclustMethod choose method to cluster distance matrix
#' @param cluster_column_name name of column that will be created with the clustering results
#' @param return_sign_loadMatrix boolean to return loadings multiplied by their sign for all PCs, default to FALSE
#' @param NAvalues how to handle NA values, remove or replace with not_detected_value
#' @param not_detected_value value to replace NA values
#' @param ... additional parameters for distance methods 
#' @return returns a fluidSCproc S3 object appended with the cluster results
#' @export
#' @details NA 
#' @examples
#' PCA_based_cluster_SC()

PCA_based_cluster_SC <- function(fluidSCproc,  based_on_values = "log2ExNorm", scaleData = T, PCAscaled = F,
                                 clustMethod = c("Kmeans","Hierarchical","Correlation"),
                                 firstPC = "PC1", secondPC = "PC2", nrClust = 2,
                                 hclustMethod = c("ward.D","ward.D2","single","complete","average","centroid"),
                                 cluster_column_name = "PCAclust", return_sign_loadMatrix = F,
                                 NAvalues = c("remove_assays","replace_with_not_detected"), not_detected_value = 0, ...) {
  
  
  ### checks ###
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  stopifnot(class(fluidSCproc) == "fluidSCproc")
  
  
  usedLoD <- fluidSCproc$proc_info$LoD
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  
  ### for merged fluidSCproc objects, what to do with NA values? ###
  NAvalues <- match.arg(NAvalues)
  
  # OR remove genes with NA values
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
  
  
  # first remove genes with zero variance!
  keepgenes <- !apply(normMatrix, MARGIN = 2, function(x) var(x,na.rm = T)) == 0
  keepgenes <- names(keepgenes)[keepgenes]
  normMatrix <- normMatrix[, colnames(normMatrix) %in% c(keepgenes)]
  
  ### scale data if wanted ###
  if(scaleData) {
    normMatrix <- scale(normMatrix)
  }
  
  
  # do Principal Component Analysis
  PCAdata <- prcomp(normMatrix, scale. = PCAscaled)
  
  # geneloadings represent contribution of variables (gene expressions) to individual principal components
  geneloadings <- PCAdata$rotation
  signMatrix <- ifelse(sign(geneloadings) == 1, 1, -1) # negative vs positive contribution
  loadingMatrix <- apply(geneloadings, MARGIN = 2, FUN = function(x) abs(x) / sum(abs(x)) * 100) # relative contribution
  sign_loadingMatrix <- signMatrix * loadingMatrix 
  
  
  # cluster data
  scores <- PCAdata$x
  clustMethod <- match.arg(clustMethod)
  
  if (clustMethod == "Kmeans") {
    K <- kmeans(scores[, c(firstPC, secondPC)], centers = nrClust, ...)
    I <- K$cluster
  }
  
  else if (clustMethod == "Hierarchical") {
    H <- hclust(dist(scores[, c(firstPC, secondPC)], ...), 
                method = hclustMethod)
    I <- cutree(H, k = nrClust)
  }
  
  else if (clustMethod == "Correlation") {
    distfun <- function(x) as.dist(1 - cor(t(x), ...))
    mycordist <- distfun(scores[, c(firstPC, secondPC)])
    hclustfun <- function(x) hclust(x, method = hclustMethod)
    Hcor <- hclustfun(mycordist)
    I <- cutree(Hcor, k = nrClust)
  }
  
  ## merge cluster data with original data ##
  prepDfr <- as.data.frame(I); prepDfr$Samples <- rownames(prepDfr); colnames(prepDfr)[[1]] <- cluster_column_name
  mergeDfr <- merge(normFluidCt, prepDfr, by = "Samples" )
  
  ifelse(return_sign_loadMatrix,  return(sign_loadingMatrix), return(fluidSCproc(mergeDfr, usedLoD)) ) 
  
}
