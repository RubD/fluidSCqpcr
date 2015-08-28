
#' Function to sparse cluster single cell data using the sparcl package
#'
#' This function is a wrapper for the sparcl package functions HierarchicalSparseCluster and KMeansSparseCluster and appends the result to the fluidSCproc object
#' or gives a vector with gene weights for clustering of the data
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param sparseClust choose method for clustering from sparcl package
#' @param wbound parameter for sparse clustering that determines feature weights
#' @param cluster_column_name name of column that will be created with the clustering results
#' @param nrClust number of clusters you want to detect (K for sparse kmeans clustering)
#' @param return_geneWeights return weights of genes to clustering instead of fluidSCproc object
#' @param selected_assays make selection of assays if wanted
#' @param NAvalues how to handle NA values, remove or replace with not_detected_value
#' @param not_detected_value value to replace NA values 
#' @param ... additional parameters for HierarchicalSparseCluster or KmeansSparseCluster from sparcle package
#' @return returns a fluidSCproc S3 object appended with the cluster results OR a vector with gene weights if return_geneWeights = TRUE
#' @export
#' @details NA 
#' @examples
#' cluster_sparse_SC()
cluster_sparse_SC <- function(fluidSCproc,  based_on_values = "log2ExNorm", sparseClust = c("Hierarchical","Kmeans"), wbound = 4, 
                              cluster_column_name = "clust", nrClust = 2, return_geneWeights = F,selected_assays = NULL,
                              NAvalues = c("remove_assays","replace_with_not_detected"), not_detected_value = 0, ...) {
  
  # load libraries
  library(sparcl)
  
  # checks
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  stopifnot(class(fluidSCproc) == "fluidSCproc")
  
  usedLoD <- fluidSCproc$proc_info$LoD
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  cluster_column_name <- cluster_column_name
  sparseClust = match.arg(sparseClust)
  
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
  
  
  scores <- normMatrix
  if(!is.null(selected_assays)) scores <- scores[, selected_assays] # create subset from assays
  
  
  
  if(sparseClust == "Hierarchical") {
    
    # sparse clustering samples hierarchical
    testClust <- HierarchicalSparseCluster(as.matrix(normMatrix), wbound = wbound, ...)
    testCut <- cutree(testClust$hc, k = nrClust)
    sampleNames <- rownames(as.matrix(normMatrix))
    clusDfr <- as.data.frame(testCut); clusDfr$Samples <- sampleNames; colnames(clusDfr)[[1]] <- cluster_column_name
    
    # weight of genes in clustering
    geneWeights <- as.vector(testClust$ws)
    names(geneWeights) <- colnames(normMatrix)
    geneWeights <- geneWeights[order(geneWeights,decreasing = T)]
    geneWdfr <- as.data.frame(geneWeights); geneWdfr$Assays <- rownames(geneWdfr)
    colnames(geneWdfr)[[1]] <- paste0(colnames(geneWdfr)[[1]],"_", cluster_column_name)
    
    # merge sample clustering and gene weight information
    normFluidCt <- merge(normFluidCt, clusDfr, by = "Samples")
    normFluidCt <- merge(normFluidCt, geneWdfr, by = "Assays")
    
    
    ifelse(return_geneWeights, return(geneWeights), return(fluidSCproc(normFluidCt, usedLoD)))
    
    
  } else if(sparseClust == "Kmeans") {
    
    # sparse clustering Kmeans
    KClust <- KMeansSparseCluster(as.matrix(normMatrix), K = nrClust, wbounds = wbound, ...)
    KmClust <- KClust[[1]]$Cs
    sampleNames <- rownames(as.matrix(normMatrix))
    clusDfr <- as.data.frame(KmClust); clusDfr$Samples <- sampleNames; colnames(clusDfr)[[1]] <- cluster_column_name
    
    # weight of genes in clustering
    geneWeights <- KClust[[1]]$ws
    geneWeights <- geneWeights[order(geneWeights,decreasing = T)]
    geneWdfr <- as.data.frame(geneWeights); geneWdfr$Assays <- rownames(geneWdfr)
    colnames(geneWdfr)[[1]] <- paste0(colnames(geneWdfr)[[1]],"_", cluster_column_name)
    
    # merge sample clustering and gene weight information
    normFluidCt <- merge(normFluidCt, clusDfr, by = "Samples")
    normFluidCt <- merge(normFluidCt, geneWdfr, by = "Assays")
    
    
    ifelse(return_geneWeights, return(geneWeights), return(fluidSCproc(normFluidCt, usedLoD)))
    
  }
  
  
  
}



#' Function to determine optimal wbound parameter for cluster_sparse_SC function
#'
#' This function is a wrapper for the sparcl package functions HierarchicalSparseCluster.permute and KMeansSparseCluster.permute
#' It will given an indication about the wbound parameter
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param sparseClust choose method for clustering from sparcl package
#' @param nperms number of permutations
#' @param wbounds vector of values for sparse clustering that determines feature weights
#' @param K number of cluster if chosen Kmeans as clustering method
#' @param selected_assays make selection of assays if wanted
#' @param NAvalues how to handle NA values, remove or replace with not_detected_value
#' @param not_detected_value value to replace NA values 
#' @param ... additional parameters for HierarchicalSparseCluster.permute or KmeansSparseCluster.permute from sparcle package
#' @return returns info about the optimal wbound parameter
#' @export
#' @details NA 
#' @examples
#' cluster_sparse_SC_permute()
cluster_sparse_SC_permute <- function(fluidSCproc, based_on_values = "log2ExNorm",sparseClust = c("Hierarchical","Kmeans"), nperms = 10, wbounds = c(1.5,2:6), K = NULL,
                                      selected_assays = NULL, NAvalues = c("remove_assays","replace_with_not_detected"), not_detected_value = 0, ...) {
  
  # load libraries
  library(sparcl)
  
  # checks
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  stopifnot(class(fluidSCproc) == "fluidSCproc")
  
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
    else normMatrix <- normMatrix[,-(which(colnames(normMatrix) %in% NAgenes))]
  } 
  # OR replace NA values with 0
  else if(NAvalues == "replace_with_not_detected") {
    normMatrix[is.na(normMatrix)] <- not_detected_value
  } else stop("Choose one of the 2 options to take care of NA values")
  
  
  if(!is.null(selected_assays)) normMatrix <- normMatrix[, selected_assays] # create subset from assays
  
  
  if(sparseClust == "Hierarchical") {
    
    perm.out <- HierarchicalSparseCluster.permute(as.matrix(normMatrix), nperms = nperms, wbounds = wbounds, ...)
    return(perm.out)
  }
  
  else if(sparseClust == "Kmeans") {
    
    if(is.null(K)) stop("When Kmeans permute selected, choose K for number of clusters, e.g. check with Kmeans explained variance")
    
    kmperm.out <- KMeansSparseCluster.permute(x = as.matrix(normMatrix), nperms = nperms, wbounds = wbounds, K = K)
    return(kmperm.out)
    
  }
  
  
}



