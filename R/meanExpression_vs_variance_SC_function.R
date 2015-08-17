
#' meanExpression_vs_variance_SC
#'
#' Function to get and plot mean abundance versus variance of expression in single cells
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "lin2Ex"
#' @param log_base transform values with log
#' @param label_deviation jitter on the y-axis for the labels
#' @param show_labels boolean to show labels, defaults to TRUE
#' @param print_graphs print graphs as output
#' @return returns a list with plot and corresponding data
#' @export
#' @details NA
#' @examples
#' meanExpression_vs_variance_SC()

meanExpression_vs_variance_SC <- function(fluidSCproc,  based_on_values = "lin2Ex", log_base = NULL, label_deviation = 0.7, show_labels = T, print_graphs = T) {
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  if(based_on_values %in% c("log2Ex","log2ExNorm") & !is.null(log_base)) warning("Creating log scaled values from already log values!")
  
  library(grid)
  library(gridExtra)
  
  normFluidCt <- fluidSCproc$data
  
  meanGeneExpression <- aggregate(normFluidCt[ ,based_on_values], by = list(normFluidCt$Assays), FUN = function(x) mean(x))
  if(!is.null(log_base)) meanGeneExpression$x <- log(meanGeneExpression$x, base = log_base)
  meanGeneVariance <- aggregate(normFluidCt[ ,based_on_values], by = list(normFluidCt$Assays), FUN = function(x) var(x))
  if(!is.null(log_base)) meanGeneVariance$x <- log(meanGeneVariance$x, base = log_base)
  
  
  testdfr <- data.frame(cbind(meanGeneExpression, meanGeneVariance))
  testdfr <- testdfr[order(testdfr$x), ]
  colnames(testdfr)[c(2,4)] <- c("abundance","variance")
  testdfr$jitterHeigth <- testdfr$variance + rep(c(label_deviation,-label_deviation), (length(testdfr$variance)))[1:length(testdfr$variance)]
  print(testdfr)
  
  pl <- ggplot(testdfr, aes(abundance, variance, label = Group.1))
  pl <- pl + geom_point(size = 4, alpha = 0.8, col = "grey")
  pl <- pl + theme_bw()
  pl <- pl + xlim(c(0, max(testdfr$abundance)))
  
  if(show_labels) pl <- pl + geom_text(aes(abundance, jitterHeigth))
  if(!is.null(log_base)) {
    pl <- pl + labs(list(x = sprintf("Average expression in single cells (log %d scaled)", log_base), y = sprintf("Single cell variance (log %d scaled)", log_base)))
  } else pl <- pl + labs(list(x = "Average expression in single cells", y ="Single cell variance "))
  
  if(print_graphs) print(pl)
  
  plot <- arrangeGrob(pl)
  
  return(list(plot = plot, data = testdfr))
  
}
