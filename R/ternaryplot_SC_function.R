

#' ternaryplot_SC function
#'
#' Function creates a ternary or triangle plot to visualize relations between 3 variables
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param sorted_selected_genes subset of genes
#' @param not_detected_value value that corresponds to not detected expression
#' @param remove_not_detected boolean to remove not detected data
#' @return returns a ternary (triangle) plot
#' @export
#' @details NA
#' @examples
#' ternaryplot_SC()

ternaryplot_SC <- function(fluidSCproc,  based_on_values = "log2Ex",sorted_selected_genes = c("Pou5f1","Psma3","Trp53"),
                           not_detected_value = 0, remove_not_detected = F) {
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  library(ggtern)
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  subMatrix <- normMatrix[,sorted_selected_genes]
  
  if(remove_not_detected) {
    
    not_detected_rows <- apply(subMatrix, 1, FUN = function(x) any(x == not_detected_value))
    subMatrix <- subMatrix[!not_detected_rows, ]
    
  }
  
  trplot <- ggtern(data = subMatrix, aes_string(x = sorted_selected_genes[1], y =  sorted_selected_genes[2], z= sorted_selected_genes[3]))
  trplot <- trplot + geom_point() + theme_tern_bw()
  return(trplot)
  
  
}
