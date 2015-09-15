

#' Function to booleanize gene expression values in subpopulations based on clustering
#'
#' Create boolean values in subpopulations based on non-parametric Wilcoxon test or with Kmeans clustering (active and inactive cluster)
#' @param fluidSCproc fluidSCproc S3 object~
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param clusterColumn name of column that contains cluster data
#' @param boolean_method create boolean values based on Wilcoxon test or Kmeans clustering of each assay
#' @param Wilcox_pval p-value treshold to classify genes as active or inactive based on one-sided Wilcoxon test ("greater")
#' @param print_heatmap print heatmap from boolean dataframe endresult
#' @return returns a dataframe with boolean values for all subpopulations
#' @export
#' @details NA 
#' @examples
#' booleanize_clusters()

booleanize_clusters <- function(fluidSCproc,  based_on_values = "log2ExNorm", clusterColumn, 
                                boolean_method = c("Wilcox_test","Kmeans_2_clusters"), Wilcox_pval = 0.4,
                                print_heatmap = T) {
  
  
  # check input
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  stopifnot(class(fluidSCproc) == "fluidSCproc")
  normFluidCt <- fluidSCproc$data
  
  boolean_method <- match.arg(boolean_method)
  
  if(boolean_method == "Wilcox_test") {
    
    # decide active or inactive relative to overall population based on Wilcoxon test: 
    # gene expression in subpopulation vs gene expression in all population (one side greater test)
    
    genes <- unique(normFluidCt$Assays[order(normFluidCt$Assays)])
    
    
    newlist <- list()
    number_groups <- length(unique(normFluidCt[ ,clusterColumn]))
    for(i in genes) {
      
      values_background <- normFluidCt[normFluidCt$Assays == i, ][ , based_on_values] # gene expression in all single cells
      
      for(j in 1:number_groups) {
        
        group_values <- normFluidCt[normFluidCt[ ,clusterColumn] == j & normFluidCt$Assays == i, ][ , based_on_values] # gene expressin in 1 subpopulation
        
        # one-sided wilcoxon test
        group_wilcox <- wilcox.test(group_values, y = values_background, alternative = "greater")
        newlist[[i]][[j]] <- group_wilcox$p.value
      }  
      
    }
    
    pvalDfr <- t(as.data.frame(newlist)); colnames(pvalDfr) <- paste0("group",1:number_groups)
    
    booleanDfr <- ifelse(pvalDfr < Wilcox_pval, 1, 0) # create booleanDfr based treshold of wilcox p-values
    
    if(print_heatmap) {
      
      library(gplots)
      number_rows <- nrow(booleanDfr)
      number_cols <- ncol(booleanDfr)
      heatmap.2(booleanDfr, scale = "none", trace = "none", col = c("white","black"), key = F, cexCol = 1,
                sepwidth = c(0.01,0.01), rowsep = 1:number_rows, colsep = 1:number_cols, cexRow = 0.5, lwid = c(1, 4), lhei = c(0.5,4))
    }
    
    
  } 
  
  
  
  
  else if(boolean_method == "Kmeans_2_clusters") {
    
    # decide active or inactive relative to overall population based on Kmeans grouping in 2 clusters 
    
    genes <- unique(normFluidCt$Assays[order(normFluidCt$Assays)])
    
    newlist <- list()
    number_groups <- length(unique(normFluidCt[,clusterColumn]))
    for(i in genes) {
      
      gene_values <- normFluidCt[normFluidCt$Assays == i, c("Samples", based_on_values, clusterColumn)]
      
      if(length(unique(gene_values[,based_on_values])) < 2 ) {
        
        # if all values are the same, i.e. usually all zero or not expressed, set gene to inactive
        
        check <- 0
        
      } else {
        
        Ktest <- kmeans(gene_values[,based_on_values], 2) # create 2 clusters for each gene expression values with Kmeans method
        
        # annotate the active and inactive cluster based on maximum expression level
        activityDfr <- Ktest$centers
        activityDfr$activity <- ifelse(activityDfr[,1] == max(activityDfr), "active","inactive")
        
        # change kmeans clusters 1 and 2 to active or inactive (depends which one has maximum expression)
        gene_values$bool <- Ktest$cluster
        gene_values$activity <- as.vector(activityDfr$activity[gene_values$bool]) 
        
        
        for(j in 1:number_groups) {
          
          # create table for number of active and inactive cells per group for 1 particular gene
          modegroup <- table(gene_values$activity[gene_values[ ,clusterColumn] == j])
          
          # if subpopulation or group has same number of inactive and active cells, set to active
          if(max(modegroup) == min(modegroup)) {
            check <- 1
          } else {
            
            # if most of the cells are active in subgroup then that particular gene is active in that subgroup
            check <- ifelse( names(modegroup[modegroup == max(modegroup)]) == "active", 1, 0)
            
          }
          newlist[[i]][[j]] <- check
        }  
        
      }
      
      
    }
    
    
    booleanDfr <- t(as.data.frame(newlist)); colnames(booleanDfr) <- paste0("group",1:number_groups)
    
    if(print_heatmap) {
      
      library(gplots)
      number_rows <- nrow(booleanDfr)
      number_cols <- ncol(booleanDfr)
      heatmap.2(booleanDfr, scale = "none", trace = "none", col = c("white","black"), key = F, cexCol = 1,
                sepwidth = c(0.01,0.01), rowsep = 1:number_rows, colsep = 1:number_cols, cexRow = 0.5, lwid = c(1, 4), lhei = c(0.5,4))
    }
    
  }
  
  return(booleanDfr)
  
}
