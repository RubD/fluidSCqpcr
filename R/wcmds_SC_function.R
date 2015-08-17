


#' Function to create Multidimensional scaling plot
#'
#' This function creates a MDS plot
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param mds_on do dimentional scaling on Samples or Assays
#' @param genes_selection make a selection for the genes to consider
#' @param cluster_column_name colors single genes based on this column
#' @param show_labels boolean to show gene / sample labels, defaults to TRUE
#' @return returns a mds plot
#' @export
#' @details NA
#' @examples
#' wcmds_SC()

wcmds_SC <- function(fluidSCproc, dimensions = 2,  based_on_values = "log2Ex", mds_on = c("Samples","Assays"),
                     genes_selection = NULL, cluster_column_name = NULL, show_labels = T) {
  
  library("vegan")
  # possible to give weights to samples (or genes)
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  mds_on <- match.arg(mds_on)
  
  normFluidCt <- fluidSCproc$data
  if(!is.null(genes_selection)) normFluidCt <- normFluidCt[normFluidCt$Assays %in% genes_selection,]
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  ## multidimensional scaling on samples
  
  if(mds_on == "Samples") {
    
    d<- dist(normMatrix)
    mds<- wcmdscale(d, k=dimensions)
    
    output <- as.data.frame(as.matrix(mds)); output$Samples <- rownames(output)
    if(!is.null(cluster_column_name)) {
      
      output <- merge(x = output, y= unique(normFluidCt[,c("Samples",cluster_column_name)]), by = "Samples")
      output[,cluster_column_name] <- factor(output[,cluster_column_name])
      
    } 
    
    
    #mOutput <- melt(output)
    
    if(!is.null(cluster_column_name)) {
      
      mdsPlot <- ggplot(output, aes_string("V1", "V2", color = cluster_column_name, label = "Samples"))
      
    } else mdsPlot <- ggplot(output, aes(V1, V2, label = Samples))
    mdsPlot <- mdsPlot + geom_point()
    if(show_labels) mdsPlot <- mdsPlot + geom_text()
    mdsPlot <- mdsPlot + theme_bw()
    return(mdsPlot)
    
  }
  
  ## multidimensional scaling on genes
  if(mds_on == "Assays") {
    
    d<- dist(t(normMatrix))
    mds<- wcmdscale(d, k=dimensions)
    
    output <- as.data.frame(as.matrix(mds)); output$Samples <- rownames(output)
    
    #mOutput <- melt(output)
    
    mdsPlot <- ggplot(output, aes(V1, V2, label = Samples))
    mdsPlot <- mdsPlot + geom_point()
    if(show_labels)mdsPlot <- mdsPlot + geom_text()
    mdsPlot <- mdsPlot + theme_bw()
    return(mdsPlot)
    
  }
  
  
  
}

