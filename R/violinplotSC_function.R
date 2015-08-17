#' violinplotSC()
#'
#' Function to visualize distribution profile of assays in single cells with violin plots
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param assayCol column name with the assays, defaults to "Assays"
#' @param assays vector with assays to show
#' @param groupColumn column name with group information, will be visualized in different facets
#' @return returns a violinplot
#' @export
#' @details NA
#' @examples
#' violinplotSC()


violinplotSC <- function(fluidSCproc,  based_on_values = "log2Ex", assayCol = "Assays", assays = c("Rpl13a", "Pou5f1"), groupColumn = NULL) {
  
  
  library(ggplot2)
  
  
  normFluidCt <- fluidSCproc$data
  normFluidCt <- normFluidCt[normFluidCt[ ,assayCol] %in% assays,]
  
  if(!is.null(groupColumn) && !groupColumn %in% colnames(normFluidCt)) stop("column name for groupColumn does not exist")
  
  if(!is.null(groupColumn)) number_rows <- length(unique(normFluidCt[,groupColumn]))
  
  
  plot <- ggplot(normFluidCt, aes_string(x = "Assays",y = "log2ExNorm"))
  plot <- plot + geom_violin()
  if(!is.null(groupColumn)) plot <- plot + facet_wrap(as.formula(paste("~", groupColumn)), nrow = number_rows)
  plot <- plot + theme_bw()
  plot
  
  
  
}

