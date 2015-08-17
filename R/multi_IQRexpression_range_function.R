

#' Function to show expression range of genes based on quantile probabilities for multiple conditions
#'
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param high_quantile high cutoff for gene expression range, based on quantile probabilities [0 - 1]
#' @param low_quantile low cutoff for gene expression range, based on quantile probabilities [0 - 1]
#' @param sored_selected_assays select subset of assays, will be displayed in same order
#' @param groupCol column name different sample groups, will be displayed in same plot
#' @param print_graphs boolean to print graphs to output
#' @return returns a list with data and plot
#' @export
#' @details NA
#' @examples
#' show_IQRexpression_range()


multi_IQRexpression_range <- function(fluidSCproc, based_on_values = "log2Ex", high_quantile = 0.75, low_quantile = 0.25,
                                      sorted_selected_assays = NULL, groupCol = NULL, print_graphs = T) {
  
  
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  normFluidCt <- fluidSCproc$data
  if(!is.null(sorted_selected_assays)) normFluidCt <- normFluidCt[normFluidCt$Assays %in% sorted_selected_assays, ]
  
  
  # function to create appropriate dataframe with ...
  create_melted_rangeDfr <- function(dfr, splitCol = NULL) {
    
    if(!is.null(splitCol)) groupNr <- unique(dfr[,splitCol])
    normMatrix <- dcast(dfr, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
    #normMatrix <- normMatrix[ ,sorted_selected_assays]
    
    # convert 0's to NA, so that you can exclude them in the quantile function
    normMatrix[normMatrix == 0] <- NA
    
    # calculate percentage that assay was not detected (ND)
    assay_ND <- apply(normMatrix, 2, FUN = function(x) sum(is.na(x)) / length(x) * 100)
    
    # define range based on 3th quartile - 1st quartile
    highIQR <- apply(normMatrix, 2, FUN = function(x) quantile(x, high_quantile, na.rm = T))
    lowIQR <- apply(normMatrix, 2, FUN = function(x) quantile(x, low_quantile, na.rm = T))
    
    
    temdfr <- as.data.frame(cbind(highIQR, lowIQR, assay_ND ))
    temdfr$rawSampleNames <- rownames(temdfr)
    temdfr$sampleNames <- paste0(rownames(temdfr)," (", format(temdfr[,"assay_ND"],digits = 1, nsmall = 1),"%)")
    temdfr$notDetected <- paste0(" (", format(temdfr[,"assay_ND"],digits = 1, nsmall = 1),"%)")
    temdfr <- temdfr[order(temdfr$assay_ND),] # orders the assays based on percentage non detected
    temdfr$sampleNames <- factor(temdfr$sampleNames, levels = temdfr$sampleNames)
    
    mtemdfr <- melt(temdfr, measure.vars = c("highIQR","lowIQR"))
    if(!is.null(splitCol)) mtemdfr <- cbind(mtemdfr, groupNr)
    return(mtemdfr)
    
  }
  
  if(is.null(groupCol)) stop("you need to provide column with group information")
  
  splitfr <- split(x = normFluidCt, f = normFluidCt[ ,groupCol])
  templist <- lapply(X = splitfr,FUN =  function(x) create_melted_rangeDfr(x, groupCol))
  
  finalmdfr <- do.call("rbind",templist)
  finalmdfr <- finalmdfr[order(finalmdfr$rawSampleNames, finalmdfr$groupNr), ]
  
  uniqGenes <- length(unique(finalmdfr$rawSampleNames))
  uniqGroups <- length(unique(finalmdfr$groupNr))
  
  # create y-values that are close to eachother for each assay
  numbering <- NULL
  for(number in 1:uniqGenes) {
    
    x <- rep(seq(from = number,to = (number + (0.1 * (uniqGroups - 1))), by = 0.1),each = 2)
    numbering <- c(numbering, x)
  }
  finalmdfr$yval <- numbering
  
  # create unique names for all conditions
  finalmdfr$sep <- paste0(finalmdfr$rawSampleNames,"_",finalmdfr$groupNr)
  
  # create graph
  tplot <- ggplot(finalmdfr, aes(x = value, y = yval, group = sep, col = rawSampleNames))
  tplot <- tplot + geom_line() + theme_bw() + labs(list(y = "", x = "log2 Expression"))
  tplot <- tplot + scale_y_continuous(labels = unique(finalmdfr$rawSampleNames)) + scale_color_discrete(guide=F)
  if(print_graphs) print(tplot)
  
  rangeplot <- arrangeGrob(tplot)
  
  return(list(data = finalmdfr, rangeplot = rangeplot ))
  
}

