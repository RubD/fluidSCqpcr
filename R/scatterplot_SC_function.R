
#' Function to create scatterplot between 2 assays
#'
#' Function creates scatterplot and shows correlation score between 2 assays (genes)
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param gene1 first selected gene, X-axis
#' @param gene2 second selected gene, Y-axis
#' @param my_formula formula to create smoothed line
#' @param se_display display standard error range around smoothed line
#' @param sample_display add sample names on plot
#' @param cor_use data to use to calculate correlation score, defaults to "complete.obs"
#' @param cor_method method for correlation score algorithm, defaults to "pearson"
#' @param smooth_method method to create smoothed line
#' @return returns a scatterplot with smoothed line and (always linear) correlation score
#' @export
#' @details NA
#' @examples
#' scatterplot_SC()

scatterplot_SC <- function(fluidSCproc,  based_on_values = "log2Ex", gene1, gene2, my_formula = "y ~ x", se_display = F, sample_display = F,
                           cor_use = "complete.obs", cor_method = "pearson", smooth_method = "lm") {
  
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  subset_normMatrix <- normMatrix[,colnames(normMatrix) %in% c(gene1, gene2) ]
  subset_normMatrix$labels <- rownames(subset_normMatrix)
  
  
  ## CORRELATION ##
  gene_x <- subset_normMatrix[,gene1]
  gene_y <- subset_normMatrix[,gene2]
  cortxt <- cor(gene_x, gene_y, use = cor_use, method = cor_method)
  
  pointx <- mean(gene_x, na.rm = T)
  pointy <- mean(gene_y, na.rm = T)
  
  
  ## SCATTERPLOT AND LINEAR MODEL ##
  sc_pl <- ggplot(subset_normMatrix, aes_string(x = gene1, y = gene2, label = "labels"))
  sc_pl <- sc_pl + geom_point(shape = 1, size = 4) + geom_smooth(method = smooth_method, formula = my_formula, se = se_display)
  if(sample_display)  sc_pl <- sc_pl + geom_text() 
  sc_pl <- sc_pl + labs(list(x = gene1, y = gene2))
  sc_pl <- sc_pl + theme_bw()
  
  if (sign(cortxt) == 1) 
    sc_pl <- sc_pl + annotate("text", x = pointx, y = pointy, 
                              label = paste("correlation: ", signif(cortxt, 2)), 
                              size = 6, color = "green")
  else sc_pl <- sc_pl + annotate("text", x = pointx, y = pointy, 
                                 label = paste("correlation: ", signif(cortxt, 2)), size = 6, 
                                 color = "red")
  
  print(sc_pl)
  
}
