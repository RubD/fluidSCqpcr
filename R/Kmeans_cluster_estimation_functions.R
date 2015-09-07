#' Helper function for Kmeans_cluster_explained_variance()
#'
#' Function to calculate the variance that is explained by performing kmeans clustering
#' @param scores matrix with scores
#' @param nrClust number of clusters to test
#' @return returns percentage of variance explained
variance_explained_Kmeans <- function(scores, nrClust) {
  
  ss <- function(x) sum(scale(x, scale = F)^2) # function to calculate sum of squares
  
  K <- kmeans(scores, centers = nrClust)
  
  # total sum of squares = sum of squares within clusters plus sum of squares between clusters
  # totss = tot.withinss + betweenss
  totss <- ss(scores) # sum of squares for all data (= 1 cluster)
  
  betweenss <- ss(K$centers[K$cluster, ]) # sum of squares between clusters
  withinss <- sapply(split(as.data.frame(scores), K$cluster), ss) # sum of squares in different clusters
  tot.withinss <- sum(withinss) # sum of sum of squares in different clusters
  
  varExpl <- (betweenss / totss) *100
  
  return(varExpl)
}


#' Helper function for Kmeans_cluster_explained_variance()
#'
#' Function to calculate the variance that is explained by performing sparse kmeans clustering
#' @param scores matrix with scores
#' @param nrClust number of clusters to test
#' @param wbound weight for kmeans clustering bound to features (genes)
#' @return returns percentage of variance explained
variance_explained_sparseKmeans <- function(scores, nrClust, wbound) {
  
  ss <- function(x) sum(scale(x, scale = F)^2) # function to calculate sum of squares
  
  sparseK <- KMeansSparseCluster(scores, K = nrClust, wbounds = wbound)
  # total sum of squares = sum of squares within clusters plus sum of squares between clusters
  # totss = tot.withinss + betweenss
  totss <- ss(scores) # sum of squares for all data (= 1 cluster)
  
  withinss <- sapply(split(as.data.frame(scores), sparseK[[1]]$Cs), ss) # sum of squares in different clusters
  tot.withinss <- sum(withinss) # sum of sum of squares in different clusters
  
  betweenss <- totss - tot.withinss
  
  varExpl <- (betweenss / totss) *100
  
  return(varExpl)
}



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
Kmeans_cluster_explained_variance <- function(fluidSCproc, based_on_values = "log2Ex", method = c("Kmeans", "sparseKmeans"),
                                              scaleData = F, test_nrClust = 10, wbound = 2,
                                              selected_assays = NULL,
                                              NAvalues = c("remove_assays","replace_with_not_detected"), not_detected_value = 0, ...) {
  
  # references
  # http://stackoverflow.com/questions/8637460/k-means-return-value-in-r
  # http://stats.stackexchange.com/questions/48520/interpreting-result-of-k-means-clustering-in-r
  
  # load libraries
  library(ggplot2)
  
  
  # perform checks
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  stopifnot(class(fluidSCproc) == "fluidSCproc")
  
  method <- match.arg(method)
  
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
  
  # scale data if wanted
  if(scaleData) {
    
    # first remove genes with zero variance!
    keepgenes <- !apply(normMatrix, MARGIN = 2, function(x) var(x,na.rm = T)) == 0
    keepgenes <- names(keepgenes)[keepgenes]
    normMatrix <- normMatrix[, colnames(normMatrix) %in% c(keepgenes)]
    
    # second scale data
    normMatrix <- scale(normMatrix)
  }
  
  scores <- normMatrix
  if(!is.null(selected_assays)) scores <- scores[, selected_assays] # create subset from assays
  
  
  ##  calculate variance explained by K clusters scores
  vec <- vector(mode = "numeric", length = test_nrClust)
  
  if(method == "Kmeans") {
    
    for(i in 1:test_nrClust) {
      
      #K <- kmeans(scores, centers = i, ...)
      #varExpl <- (K$betweenss / K$totss) * 100
      varExpl <- variance_explained_Kmeans(scores = scores, nrClust = i)
      vec[i] <- round(varExpl,1)
      
    }
    
  } else if(method == "sparseKmeans") {
    
    print("OK1")
    for(i in 2:test_nrClust) {
      
      varExpl <- variance_explained_sparseKmeans(scores = scores, nrClust = i, wbound = wbound)
      vec[i] <- round(as.numeric(varExpl),1)
      
    }
    
  } else stop("choose method 'Kmeans' or 'sparseKmeans'")
  
  
  
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



#' Helper function for Kmeans_cluster_statistics()
#'
#' @param inputMatriw matrix with scores
#' @param sparseK number of clusters for sparse Kmeans clustering
#' @param wbound weight for kmeans clustering bound to features (genes)
#' @return returns list that can be used with clusGap function
clusGap_sparseKmeans <- function(inputMatrix, sparseK, wbound = 2) {
  
  result <- KMeansSparseCluster(inputMatrix, K = sparseK, wbounds = wbound)
  return(list(cluster = result[[1]]$Cs)) 
  
}


#' Helper function for Kmeans_cluster_statistics()
#'
#' @param inputMatriw matrix with scores
#' @param sparseK number of clusters for sparse Kmeans clustering
#' @param wbound weight for kmeans clustering bound to features (genes)
#' @return returns list that can be used with clusterboot function
sparseKmeansCBI <- function(inputMatrix, sparseK, wbound = 2 ) {
  
  result <- KMeansSparseCluster(inputMatrix, K = sparseK, wbounds = wbound)
  
  partition <- result[[1]]$Cs
  cl <- list()
  nc <- sparseK
  for(i in 1:nc) cl[[i]] <- partition == i
  out <- list(result = result, nc = nc, clusterlist = cl, partition = partition,
              clustermethod = "sparseKmeans")
  return(out)
  
  
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
Kmeans_cluster_statistics <- function(fluidSCproc, based_on_values = "log2Ex", scaleData = F, test_nrClust = 10,
                                      method = c("Kmeans", "sparseKmeans"), wbound = 2,
                                      selected_assays = NULL, NAvalues = c("remove_assays","replace_with_not_detected"),
                                      not_detected_value = 0,
                                      clusGap_B = 100, clusterboot_B = 100, ...) {
  
  
  # load libraries
  library(cluster)
  library(fpc)
  
  method <- match.arg(method)
  
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
  
  # scale data if wanted
  if(scaleData) {
    
    # first remove genes with zero variance!
    keepgenes <- !apply(normMatrix, MARGIN = 2, function(x) var(x,na.rm = T)) == 0
    keepgenes <- names(keepgenes)[keepgenes]
    normMatrix <- normMatrix[, colnames(normMatrix) %in% c(keepgenes)]
    
    # second scale data
    normMatrix <- scale(normMatrix)
  }
  
  if(!is.null(selected_assays)) normMatrix <- normMatrix[, selected_assays] # create subset from assays
  
  
  ## 1. calculate Gap distance
  if(method == "Kmeans") {
    myGap <- clusGap(x = as.matrix(normMatrix), FUNcluster = kmeans, K.max = test_nrClust, B = clusGap_B)
  }
  else if(method == "sparseKmeans") {
    myGap <- clusGap(x = as.matrix(normMatrix), FUNcluster = clusGap_sparseKmeans, wbound = wbound, K.max = test_nrClust, B = clusGap_B)
  }
  
  
  print("clusGAP DISTANCE:")
  print(myGap, method = "firstSEmax")
  print(myGap, method = "firstmax")
  print(myGap, method = "Tibs2001SEmax")
  print(myGap, method = "globalSEmax")
  print(myGap, method = "globalmax")
  
  print("")
  print("")
  
  
  ## 2. bootstrap method
  if(method == "Kmeans") {
    km_boot <- clusterboot(data = as.matrix(normMatrix), B = clusterboot_B, clustermethod = kmeansCBI, krange = test_nrClust)
    
  }
  else if(method == "sparseKmeans") {
    km_boot <- clusterboot(data = as.matrix(normMatrix), B = clusterboot_B, clustermethod = sparseKmeansCBI, sparseK = test_nrClust, wbound = wbound)
    
  }
  
  
  print("BOOTSTRAP METHOD:")
  print(km_boot)
  
  return(list(gap_distance = myGap, km_boot = km_boot ))
}