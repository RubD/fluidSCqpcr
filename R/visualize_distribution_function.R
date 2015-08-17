
#' Function to visualize distribution of raw or normalized data
#' 
#' This function visualized distribution of raw or normalized data with the use of a boxplot, histogram or density plot
#' @param fluidSCraw fluidSCraw S3 object
#' @param based_on_values values to use for distribution visualization
#' @param mode visualization method, defaults to boxplot
#' @param groupColumn create facets based on column with grouping information
#' @return creates a visualization of the distribution
#' @export
#' @examples
#' visualize_distribution()

visualize_distribution <- function(fluidSCraw, based_on_values = "Value", mode = c("boxplot","histogram","density"), groupColumn = NULL) {
  
  library(ggplot2)
  
  rawFluidCt <- fluidSCraw$data
  
  if(!is.null(groupColumn)) { # different groups
    
    if(mode == "boxplot") {
      pl <- ggplot(rawFluidCt[rawFluidCt[ ,based_on_values] < 999, ], aes_string(x = groupColumn,y = based_on_values))
      pl <- pl + geom_boxplot() + theme_bw()
      return(pl)
    }
    
    
    if(mode == "histogram") {
      pl <- ggplot(rawFluidCt[rawFluidCt[ ,based_on_values] < 999, ], aes_string(x = based_on_values, fill = groupColumn))
      pl <- pl + geom_histogram(alpha = 0.4, position = "identity") + theme_bw() 
      return(pl)
    }
    
    
    if(mode == "density") {
      pl <- ggplot(rawFluidCt[rawFluidCt[ ,based_on_values] < 999, ], aes_string(x = based_on_values, fill = groupColumn))
      pl <- pl + geom_density(alpha = 0.4, position = "identity") + theme_bw() 
      return(pl)
    }
    
    
  } else { # no different groups / conditions / ...
    
    if(mode == "boxplot") {
      pl <- ggplot(rawFluidCt[rawFluidCt[ ,based_on_values] < 999, ], aes_string(x = 1, y = based_on_values))
      pl <- pl + geom_boxplot() + theme_bw()
      return(pl)
    }
    
    
    if(mode == "histogram") {
      pl <- ggplot(rawFluidCt[rawFluidCt[ ,based_on_values] < 999, ], aes(x = based_on_values))
      pl <- pl + geom_histogram(alpha = 0.4, position = "identity") + theme_bw() 
      return(pl)
    }
    
    
    if(mode == "density") {
      pl <- ggplot(rawFluidCt[rawFluidCt[ ,based_on_values] < 999, ], aes(x = based_on_values))
      pl <- pl + geom_density(alpha = 0.4, position = "identity") + theme_bw() 
      return(pl)
    }
    
    
  }
  
  
  
  
}
