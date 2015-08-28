
#' Function to visualize explained variance by different Kmeans cluster. Creates scree (elbow) plot.
#'
#' This function will perform Kmeans clustering on your fluidSCproc object for a range of defined clusters (1 to test_nrClust) and 
#' the output will be a scree-plot, where you see the additional variance that is explained with each cluster.
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param test_nrClust max number of clusters to test
#' @param selected_assays make selection of assays if wanted
#' @param NAvalues how to handle NA values, remove or replace with not_detected_value
#' @param not_detected_value value to replace NA values 
#' @param ... additional parameters for HierarchicalSparseCluster.permute or KmeansSparseCluster.permute from sparcle package
#' @return returns info about the optimal wbound parameter
#' @export
#' @details NA 
#' @examples
#' Kmeans_cluster_explained_variance()
Kmeans_cluster_explained_variance <- function(fluidSCproc, based_on_values = "log2Ex", test_nrClust = 10,
                                              selected_assays = NULL, NAvalues = c("remove_assays","replace_with_not_detected"), not_detected_value = 0,...) {
  
  
  # load libraries
  library(ggplot2)
  
  
  # perform checks
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
  
  scores <- normMatrix
  if(!is.null(selected_assays)) scores <- scores[, selected_assays] # create subset from assays
  
  
  ##  calculate variance explained by K clusters scores
  vec <- vector(mode = "numeric", length = test_nrClust)
  
  for(i in 1:test_nrClust) {
    
    K <- kmeans(scores, centers = i, ...)
    varExpl <- (K$betweenss / K$totss) * 100
    vec[i] <- round(varExpl,1)
    
  }
  
  myfr <- data.frame(nrClust = 1:test_nrClust, varExplained = vec)
  
  # calculate change per step
  changevec <- rep(NA, nrow(myfr) -1)
  
  for(i in nrow(myfr):2) {
    
    changevec[i] <- myfr$varExplained[i] - myfr$varExplained[i-1]
  }
  
  myfr$change <- changevec
  myfr$varExplained[2] - myfr$varExplained[1] 
  
  # create ggplot graph
  pl <- ggplot(myfr, aes(x = nrClust, y = varExplained))
  pl <- pl + geom_line() + geom_point(size = 4 ) + scale_x_continuous(breaks = 1:test_nrClust) + theme_bw()
  pl <- pl + geom_bar(data = myfr, aes(x = nrClust, y = change), stat = "identity", alpha = 0.4)
  pl <- pl + labs(list(x = "# of clusters", y = "% variance explained (between_ss / total_ss * 100)"))
  
  return(pl)
  
}

#' Function to determine best number of clusters in kmeans clustering
#'
#' This function will calculate gap-statistics for range of K clusters and it will also output results of a bootstrap function on range of kmeans clustering results.
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param test_nrClust max number of clusters to test
#' @param clusGap_B number of monte carlo ("bootstrap") samples for gap statistic
#' @param clusterboot_B Number of resampling runs for each scheme for bootstrap method, see clusterboot
#' @param selected_assays make selection of assays if wanted
#' @param NAvalues how to handle NA values, remove or replace with not_detected_value
#' @param not_detected_value value to replace NA values 
#' @return returns and prints estimations for optimal number of clusters based on 2 different "bootstrap" algorithms (functions)
#' @export
#' @details NA 
#' @examples
#' Kmeans_cluster_statistics()
Kmeans_cluster_statistics <- function(fluidSCproc, based_on_values = "log2Ex", test_nrClust = 10, clusGap_B = 500, clusterboot_B = 20,
                                      selected_assays = NULL, NAvalues = c("remove_assays","replace_with_not_detected"), not_detected_value = 0) {
  
  
  # load libraries
  library(cluster)
  library(fpc)
  
  # perform checks
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
  
  scores <- normMatrix
  if(!is.null(selected_assays)) scores <- scores[, selected_assays] # create subset from assays
  
  
  ## 1. calculate Gap distance
  myGap <- clusGap(x = as.matrix(normMatrix), FUNcluster = kmeans, K.max = test_nrClust, B = clusGap_B)
  
  print("GAP DISTANCE:")
  print(myGap, method = "firstSEmax")
  print(myGap, method = "firstmax")
  print(myGap, method = "Tibs2001SEmax")
  print(myGap, method = "globalSEmax")
  print(myGap, method = "globalmax")
  
  print("")
  print("")
  
  ## 2. bootstrap method
  km_boot <- clusterboot(data = as.matrix(normMatrix), B = clusterboot_B, clustermethod = kmeansCBI, krange = test_nrClust)
  print("BOOTSTRAP METHOD:")
  print(km_boot)
  
  return(list(gap_distance = myGap, km_boot = km_boot ))
}

