
#' meanExpression_vs_standardDeviation_SC
#'
#' Function to get and plot mean abundance versus standard deviation of expression in single cells
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "lin2Ex"
#' @param selected_genes subset of genes
#' @param log_base transform values with log
#' @param label_deviation jitter on the y-axis for the labels
#' @param show_labels boolean to show labels, defaults to TRUE
#' @param log2Ex_maximum_COV boolean to draw line with hypothetical maxium Coefficient of Variation (COV)
#' @param log2Ex_COV_abline value to draw a line with constant COV on plot, defaults to 0.25
#' @param print_graphs print graphs as output
#' @return returns a list with plot and corresponding data
#' @export
#' @details NA
#' @examples
#' meanExpression_vs_standardDeviation_SC()
#' 
meanExpression_vs_standardDeviation_SC <- function(fluidSCproc,  based_on_values = "lin2Ex",selected_genes = NULL, log_base = NULL, label_deviation = 0.7, show_labels = T,
                                                   log2Ex_maximum_COV = T, log2Ex_COV_abline = 0.25, print_graphs = T) {
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  library(grid)
  library(gridExtra)
  
  if(based_on_values %in% c("log2Ex","log2ExNorm","log2ExNormMerge") & !is.null(log_base)) warning("Creating log scaled values from already log values!")
  
  normFluidCt <- fluidSCproc$data
  if(!is.null(selected_genes)) normFluidCt <- normFluidCt[normFluidCt$Assays %in% selected_genes, ]
  
  meanGeneExpression <- aggregate(normFluidCt[ ,based_on_values], by = list(normFluidCt$Assays), FUN = function(x) mean(x))
  if(!is.null(log_base)) meanGeneExpression$x <- log(meanGeneExpression$x, base = log_base)
  meanGeneStandardDeviation <- aggregate(normFluidCt[ ,based_on_values], by = list(normFluidCt$Assays), FUN = function(x) sd(x))
  if(!is.null(log_base)) meanGeneStandardDeviation$x <- log(meanGeneStandardDeviation$x, base = log_base)
  
  
  testdfr <- data.frame(cbind(meanGeneExpression, meanGeneStandardDeviation))
  testdfr <- testdfr[order(testdfr$x), ]
  colnames(testdfr)[c(2,4)] <- c("abundance","standardDeviation")
  testdfr$jitterHeigth <- testdfr$standardDeviation + rep(c(label_deviation,-label_deviation), (length(testdfr$standardDeviation)))[1:length(testdfr$standardDeviation)]
  
  
  # calculate maximum theoretic COV see example: "Single-cell transcriptomics reveals bimodality in expression and splicing in immune cells"
  nrCells <- length(fluidSCproc$unique_samples)
  my_intercept <- sqrt(nrCells / (nrCells-1)) * log(2)/2
  my_rico <- sqrt(nrCells / (nrCells - 1))
  
  
  # scatterplot
  pl <- ggplot(testdfr, aes(abundance, standardDeviation, label = Group.1))
  pl <- pl + geom_point(size = 4, alpha = 0.8, col = "grey")
  pl <- pl + theme_bw()
  pl <- pl + xlim(c(0, max(testdfr$abundance)))
  
  if(based_on_values %in% c("log2Ex","log2ExNorm","log2ExNormMerge")) {
    pl <- pl + geom_abline(intercept= 0, slope = log2Ex_COV_abline, colour = "red", linetype = 2)
    if(log2Ex_maximum_COV) pl <- pl + geom_abline(intercept= my_intercept,slope = my_rico, colour = "blue", linetype = 2)
  } else warning("line can only be drawn for log2Ex or log2ExNorm values")
  
  
  if(show_labels) pl <- pl + geom_text(aes(abundance, jitterHeigth))
  if(!is.null(log_base)) {
    pl <- pl + labs(list(x = sprintf("Average expression in single cells (log %d scaled)", log_base), y = sprintf("Single cell standard deviation (log %d scaled)", log_base)))
  } else pl <- pl + labs(list(x = "Average expression in single cells", y ="Single cell standard deviation "))
  if(print_graphs) print(pl)
  
  plot <- arrangeGrob(pl)
  
  return(list(plot = plot, data = testdfr))
  
}
