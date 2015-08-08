
#' Function to identify possible outlier samples
#'
#' This function searches for possible outliers based on an Inter Array Correlation (IAC) score and highlight these based on quantile probability score
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use for IAC score, defaults to "log2Ex"
#' @param color_outliers_probability cutoff to color outliers based on quantile probability distribution
#' @return returns a boxplot with sample names and possible outliers colored in red
#' @export
#' @details 
#' 
#' @references Michael C. Oldham et al. "Functional Organization of the Transcriptome in Human Brain" Supplement
#' @examples
#' identifyOutlierSamples()

identifyOutlierSamples <- function(fluidSCproc, based_on_values = "log2Ex", color_outliers_probability = 0.05) {
  
  library(reshape2)
  library(ggplot2)
  
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  # calculate Inter Array (Cell) Correlation
  IAC = cor(t(normMatrix))
  
  ## VISUALIZE OUTLIERS ##
  # histogram
  hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
  
  # dendrogram
  cluster1=hclust(as.dist(1-IAC),method="average")
  plot(cluster1,cex=0.7,labels=dimnames(t(normMatrix))[[2]])
  
  # boxplot of IAC Z-scores
  meanIAC=apply(IAC,2,mean)
  sdCorr=sd(meanIAC)
  numbersd=(meanIAC-mean(meanIAC))/sdCorr
  
  outlierDfr <- as.data.frame(numbersd)
  outlierDfr$Samples <- names(numbersd)
  
  # add jitter to the data
  outlierDfrJit <- outlierDfr
  outlierDfrJit$xj <- jitter(rep(1,nrow(outlierDfrJit)), amount = 0.3)
  # identify possible outliers based on probability
  outlierDfrJit$suspect <- outlierDfrJit$numbersd < quantile(outlierDfrJit$numbersd, probs = color_outliers_probability)
  
  library(ggplot2)
  pl <- ggplot(data = outlierDfrJit, aes(x = xj, y = numbersd, label = Samples))
  pl <- pl + geom_boxplot()
  pl <- pl + geom_point()
  pl <- pl + theme_bw()
  pl <- pl + geom_text(aes(col = suspect))
  pl <- pl + scale_color_manual(values = c("darkgrey","red"))
  pl <- pl + labs(list(x = "", y = "Z-score deviation from average Inter Cell Correlation score"))
  pl
  
}
