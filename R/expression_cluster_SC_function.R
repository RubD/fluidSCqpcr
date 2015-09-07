#' Function to to create barplots from different cluster groups
#'
#' This function will calculate deviation scores for the different cluster groups from the population average and display them as bar plots.
#' In addition it will also cluster the assays.
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param cluster_column_name name of the column with cluster information
#' @param centrality_function function to use for calculating a centrality score, e.g. mean, median, ...
#' @param log_base_adjust transform deviation scores
#' @param cluster_genes algorithm to cluster assays
#' @param nrClust number of assay clusters
#' @param distMethod method to calculate distance matrix
#' @param clustMethod method to cluster distance matrix
#' @param geneWeights named assay vector with weights, will be used to highlight more important genes
#' @param quantile_breaks divide geneWeights in quantile groups
#' @param quantile_labels give labels to groups created by quantile_breaks
#' @param alpha_range vector with range for alpha-value
#' @param selected_assays select subset of assays
#' @return returns a barplot with facets for different clusters and important assays highlighted based on assay weights
#' @export
#' @details NA 
#' @examples
#' expression_cluster_SC()
expression_cluster_SC <- function(fluidSCproc, based_on_values = "log2ExNorm", cluster_column_name = "clust", centrality_function = mean, log_base_adjust = 10,
                                  cluster_genes = c("Hierarchical", "Correlation"),nrClust = 6, distMethod = "euclidean", clustMethod = "average",
                                  geneWeights = NULL, quantile_breaks = c(0.25, 0.5, 0.75, 0.95, 1), quantile_labels = c("very low","low","medium","high","very high"),
                                  alpha_range = c(0.5,1), selected_assays = NULL ) {
  
  # load libraries
  library(ggplot2)
  
  # perform checks
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  stopifnot(class(fluidSCproc) == "fluidSCproc")
  
  
  cluster_genes <- match.arg(cluster_genes)
  normFluidCt <- fluidSCproc$data
  
  ## TEST: if geneWeights provided, only continue with subet of data
  if(!is.null(geneWeights)) normFluidCt <- normFluidCt[normFluidCt$Assays %in% names(geneWeights), ]
  
  # split dataframe according to cluster column
  splitlist <- split(normFluidCt, f = normFluidCt[,cluster_column_name])
  
  # calculate overall mean per assay (gene)
  overall_mean <- aggregate(normFluidCt[,based_on_values], by = list(normFluidCt$Assays), FUN = function(x) centrality_function(x) )
  
  
  # calculate centrality scores for all assays in the different cluster groups
  mergelist <- list()
  for(split in 1:length(splitlist)) {
    
    nr_cells <- nrow(splitlist[[split]]) / length(unique(as.character(splitlist[[split]]$Assays))) # complex way to find number of cells
    
    temp <- splitlist[[split]]
    aggr_temp <- aggregate(temp[, based_on_values], by = list(temp$Assays), FUN = function(x) centrality_function(x))
    aggr_temp$x <- aggr_temp$x - overall_mean$x
    groupname <- paste0("group",split," ",nr_cells," cells")
    aggr_temp$group <- groupname
    mergelist[[split]] <- aggr_temp
  }
  
  
  mergetest <- do.call("rbind", mergelist)
  mergetest$y = ifelse(mergetest$x == 0, 0, ifelse(mergetest$x > 0, log((mergetest$x + 1),log_base_adjust), -log(((-1 * mergetest$x) +1),log_base_adjust)))
  mergecast <- dcast(mergetest, formula = Group.1 ~ group, value.var = "y"); rownames(mergecast) <- mergecast$Group.1; mergecast <- mergecast[, -1]
  
  
  if(cluster_genes == "Hierarchical") {
    
    # for hierarchical clustering of genes
    mydist <- dist(as.matrix(mergecast), method = distMethod)
    myclust <- hclust(mydist, method = clustMethod)
    labels_in_order <- myclust$labels[myclust$order]
    
    mergetest$Group.1 <- factor(mergetest$Group.1, levels = labels_in_order)
    
    # create gene groups
    genegroups <- cutree(myclust, k = nrClust)
    genegroups <- as.data.frame(genegroups)
    genegroups$genes <- rownames(genegroups)
    
    
    modmergetest <- merge(x = mergetest, y = genegroups, by.x = "Group.1", by.y =  "genes")
    modmergetest$genegroups <- factor(modmergetest$genegroups)
    
    if(!is.null(geneWeights)) {
      
      geneWdfr <- as.data.frame(geneWeights); geneWdfr$genes <- rownames(geneWdfr)
      
      mod2mergetest <- merge(modmergetest, geneWdfr, by.x = "Group.1", by.y  = "genes" )
      
      saveQuantiles <- as.vector(quantile(mod2mergetest$geneWeights, probs = quantile_breaks))
      
      mod2mergetest$gene_weight <- cut(mod2mergetest$geneWeights,
                                       breaks = c(saveQuantiles), 
                                       include.lowest = T, right = T, labels = quantile_labels)
      
      
      if(!is.null(selected_assays)) mod2mergetest <- mod2mergetest[mod2mergetest$Group.1 %in% selected_assays, ]
      
      
      pl <- ggplot(mod2mergetest, aes(x = Group.1, y = y, fill = genegroups, alpha = gene_weight))
      pl <- pl + geom_bar(stat = "identity") + facet_wrap(~ group, ncol = 1 )
      pl <- pl  + theme_bw() + theme(axis.text.x = element_text(angle = 45)) 
      pl <- pl + scale_fill_discrete() + scale_alpha_discrete(range = alpha_range, guide = guide_legend(reverse = T))
      pl
      
      
    } else {
      
      if(!is.null(selected_assays)) modmergetest <- modmergetest[modmergetest$Group.1 %in% selected_assays, ]
      
      pl <- ggplot(modmergetest, aes(x = Group.1, y = y, fill = genegroups))
      pl <- pl + geom_bar(stat = "identity") + facet_wrap(~ group, ncol = 1 )
      pl <- pl  + theme_bw() + theme(axis.text.x = element_text(angle = 45)) 
      pl <- pl + scale_fill_discrete()
      pl
      
    }
    
    
    
  }
  
  else if(cluster_genes == "Correlation") {
    
    # for corr clust
    varvec <- !apply(mergecast, 1, FUN = function(x) var(x)) == 0
    keepgenes <- names(varvec)[varvec]
    
    modmergecast <- mergecast[rownames(mergecast) %in% keepgenes, ]
    mydist <- as.dist(1 - cor(as.matrix(t(modmergecast)), method = distMethod))
    myclust <- hclust(mydist, method = clustMethod)
    
    labels_in_order <- myclust$labels[myclust$order]
    
    modmergetest <- mergetest[mergetest$Group.1 %in% keepgenes, ]
    
    genegroups <- cutree(myclust, k = nrClust)
    genegroups <- as.data.frame(genegroups)
    genegroups$genes <- rownames(genegroups)
    
    modmergetest <- merge(x = modmergetest, y = genegroups,by.x = "Group.1", by.y =  "genes")
    
    modmergetest$Group.1 <- factor(modmergetest$Group.1, levels = labels_in_order)
    modmergetest$genegroups <- factor(modmergetest$genegroups)
    
    
    if(!is.null(geneWeights)) {
      
      geneWdfr <- as.data.frame(geneWeights); geneWdfr$genes <- rownames(geneWdfr)
      
      
      
      mod2mergetest <- merge(modmergetest, geneWdfr, by.x = "Group.1", by.y  = "genes" )
      
      print(head(mod2mergetest))
      
      
      saveQuantiles <- as.vector(quantile(mod2mergetest$geneWeights, probs = quantile_breaks))
      
      print(saveQuantiles)
      
      mod2mergetest$gene_weight <- cut(mod2mergetest$geneWeights,
                                       breaks = c(saveQuantiles), 
                                       include.lowest = T, right = T, labels = quantile_labels)
      
      
      if(!is.null(selected_assays)) mod2mergetest <- mod2mergetest[mod2mergetest$Group.1 %in% selected_assays, ]
      
      
      pl <- ggplot(mod2mergetest, aes(x = Group.1, y = y, fill = genegroups, alpha = gene_weight))
      pl <- pl + geom_bar(stat = "identity") + facet_wrap(~ group, ncol = 1 )
      pl <- pl  + theme_bw() + theme(axis.text.x = element_text(angle = 45)) 
      pl <- pl + scale_fill_discrete() + scale_alpha_discrete(range = alpha_range, guide = guide_legend(reverse = T))
      pl
      
      
    } else {
      
      if(!is.null(selected_assays)) modmergetest <- modmergetest[modmergetest$Group.1 %in% selected_assays, ]
      
      pl <- ggplot(modmergetest, aes(x = Group.1, y = y, fill = genegroups))
      pl <- pl + geom_bar(stat = "identity") + facet_wrap(~ group, ncol = 1 )
      pl <- pl  + theme_bw() + theme(axis.text.x = element_text(angle = 45)) 
      pl <- pl + scale_fill_discrete()
      pl
      
    }  
    
  }
  
  
  return(pl)
  
}
