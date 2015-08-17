
#' dotplotSC
#'
#' Function to visualize distribution profile of assays in single cells with a dotplot
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2Ex"
#' @param assayCol column name with the assays, defaults to "Assays"
#' @param assays vector with assays to show
#' @param missingCol optional column name with information for missing values (or other features)
#' @param facetCol column to create facets
#' @param coord_ratio/binwidth graphical parameters that need to be adjusted for proper visualization (should have the same value ideally)
#' @param ylimit maximum height of y-axis
#' @return returns a dotplot
#' @export
#' @details NA
#' @examples
#' dotplotSC()

dotplotSC <- function(fluidSCproc,  based_on_values = "log2Ex", assayCol = "Assays", assays = c("Rpl13a", "Pou5f1"),
                      missingCol = NULL, facetCol = NULL, coord_ratio = 0.6, binwidth = 0.6, ylimit = 25) {
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  normFluidCt <- fluidSCproc$data
  normFluidCt <- normFluidCt[normFluidCt[ ,assayCol] %in% assays,]
  max_x_limit <- max(normFluidCt[ ,based_on_values]) + 2
  
  require(ggplot2)
  
  if(is.null(facetCol)) {
    
    if (!is.null(missingCol)) {
      myplot <- ggplot(normFluidCt, aes_string(x = based_on_values, fill = missingCol))
      myplot <- myplot + geom_dotplot(binwidth = binwidth)
      myplot <- myplot + ylim(0, ylimit) + xlim(0, max_x_limit)
      myplot <- myplot + coord_fixed(ratio = coord_ratio)
      myplot <- myplot + facet_wrap(as.formula(paste("~", assayCol)))
      myplot <- myplot + theme_bw()
      myplot <- myplot + theme(axis.title = element_text(size = 22))
      myplot <- myplot + labs(x = based_on_values, y = "number of cells")
    }
    else if (is.null(missingCol)) {
      myplot <- ggplot(normFluidCt, aes_string(x = based_on_values))
      myplot <- myplot + geom_dotplot(binwidth = binwidth)
      myplot <- myplot + scale_fill_identity(guide = guide_legend(title = "Expression"))
      myplot <- myplot + ylim(0, ylimit) + xlim(0, max_x_limit)
      myplot <- myplot + coord_fixed(ratio = coord_ratio)
      myplot <- myplot + facet_wrap(as.formula(paste("~", assayCol)))
      myplot <- myplot + theme_bw()
      myplot <- myplot + theme(axis.title = element_text(size = 22))
      myplot <- myplot + labs(x = based_on_values, y = "number of cells")
    }
    
  }
  
  else if (!is.null(facetCol)) {
    
    if (!is.null(missingCol)) {
      myplot <- ggplot(normFluidCt, aes_string(x = based_on_values, fill = missingCol))
      myplot <- myplot + geom_dotplot(binwidth = binwidth)
      myplot <- myplot + ylim(0, ylimit) + xlim(0, max_x_limit)
      myplot <- myplot + coord_fixed(ratio = coord_ratio)
      myplot <- myplot + facet_grid(as.formula(paste(facetCol,"~", assayCol)))
      myplot <- myplot + theme_bw()
      myplot <- myplot + theme(axis.title = element_text(size = 22))
      myplot <- myplot + labs(x = based_on_values, y = "number of cells")
    }
    else if (is.null(missingCol)) {
      myplot <- ggplot(normFluidCt, aes_string(x = based_on_values))
      myplot <- myplot + geom_dotplot(binwidth = binwidth)
      myplot <- myplot + scale_fill_identity(guide = guide_legend(title = "Expression"))
      myplot <- myplot + ylim(0, ylimit) + xlim(0, max_x_limit)
      myplot <- myplot + coord_fixed(ratio = coord_ratio)
      myplot <- myplot + facet_grid(as.formula(paste(facetCol,"~", assayCol)))
      myplot <- myplot + theme_bw()
      myplot <- myplot + theme(axis.title = element_text(size = 22))
      myplot <- myplot + labs(x = based_on_values, y = "number of cells")
    }
    
  }
  
  print(myplot)
  
}



