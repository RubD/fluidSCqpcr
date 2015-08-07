#' Generic function to remove data 
#' 
#' This function removes Samples or Assays from fluidSCraw or fluidSCproc object
#' @param x S3 Object fluidSCraw or fluidSCproc
#' @param samples samples to remove
#' @param assays to remove
#' @keywords remove
#' @return Object of type fluidSCraw or fluidSCproc
#' @export
#' @examples
#' remove_data(x, samples, assays)
remove_data <- function(x, samples, assays) UseMethod("remove_data")


#' Method for new generic function remove_data() 
#' 
#' This function removes Samples or Assays from fluidSCraw object
#' @param x S3 Object fluidSCraw
#' @param samples samples to remove
#' @param assays to remove
#' @keywords remove fluidSCraw
#' @return Object of type fluidSCraw
#' @export
#' @examples
#' remove_data(x, samples, assays)
remove_data.fluidSCraw <- function(x, samples = NULL, assays = NULL) {
  
  dfr <- x$data
  dfr <- dfr[!dfr$Assays %in% assays, ]
  dfr <- dfr[!dfr$Samples %in% samples, ]
  
  fluidSCraw(dfr)
  
}


#' Method for new generic function remove_data() 
#' 
#' This function removes Samples or Assays from fluidSCproc object
#' @param x S3 Object fluidSCproc
#' @param samples samples to remove
#' @param assays to remove
#' @keywords remove fluidSCproc
#' @return Object of type fluidSCproc
#' @export
#' @examples
#' remove_data(x, samples, assays)
remove_data.fluidSCproc <- function(x, samples = NULL, assays = NULL) {
  
  usedLoD <- x$proc_info$LoD
  dfr <- x$data
  dfr <- dfr[!dfr$Assays %in% assays, ]
  dfr <- dfr[!dfr$Samples %in% samples, ]
  
  fluidSCproc(dfr, usedLoD)
  
}
