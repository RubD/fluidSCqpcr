
#' dotplotSC
#'
#' Function to visualize distribution profile of assays in single cells with a dotplot
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param assayCol column name with the assays, defaults to "Assays"
#' @param assays vector with assays to show
#' @param missingCol optional column name with information for missing values (or other features)
#' @param coord_ratio/binwidth graphical parameters that need to be adjusted for proper visualization (should have the same value ideally)
#' @return returns a dotplot
#' @export
#' @details 
#' @examples
#' dotplotSC()

dotplotSC <- function (fluidSCproc,  based_on_values = "log2Ex", assayCol = "Assays", assays = c("Rpl13a", "Pou5f1"), missingCol = NULL, coord_ratio = 0.6, binwidth = 0.6) 
{
  
  normFluidCt <- fluidSCproc$data
  normFluidCt <- normFluidCt[normFluidCt[ ,assayCol] %in% assays,]
  max_x_limit <- max(normFluidCt[ ,based_on_values]) + 2
  
  require(ggplot2)
  if (!is.null(missingCol)) {
    plot <- ggplot(normFluidCt, aes_string(x = based_on_values, fill = missingCol))
    plot <- plot + geom_dotplot(binwidth = binwidth)
    plot <- plot + ylim(0, 25) + xlim(0, max_x_limit)
    plot <- plot + coord_fixed(ratio = coord_ratio)
    plot <- plot + facet_wrap(as.formula(paste("~", assayCol)))
    plot <- plot + theme_bw()
    plot <- plot + theme(axis.title = element_text(size = 22))
    plot <- plot + labs(x = based_on_values, y = "number of cells")
  }
  else if (is.null(missingCol)) {
    plot <- ggplot(normFluidCt, aes_string(x = based_on_values))
    plot <- plot + geom_dotplot(binwidth = binwidth)
    plot <- plot + scale_fill_identity(guide = guide_legend(title = "Expression"))
    plot <- plot + ylim(0, 25) + xlim(0, max_x_limit)
    plot <- plot + coord_fixed(ratio = coord_ratio)
    plot <- plot + facet_wrap(as.formula(paste("~", assayCol)))
    plot <- plot + theme_bw()
    plot <- plot + theme(axis.title = element_text(size = 22))
    plot <- plot + labs(x = based_on_values, y = "number of cells")
  }
  print(plot)
  
}

