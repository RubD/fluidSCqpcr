
#' Function to merge fluidSCraw objects
#'
#' This function merges multiple fluidSCraw S3 objects 
#' @param fluidSCraw_list list of fluidSCraw S3 objects, need to have the same column names
#' @return returns a fluidSCraw S3 object
#' @export
#' @examples
#' merge_fluidSCraw()

merge_fluidSCraw <- function(fluidSCraw_list) {
  
  # merge different fluidSCproc S3 objects with rbind
  datalist <- lapply(fluidSCraw_list, FUN = function(x) x$data)
  mergedFluidSCraw <- do.call("rbind", datalist)
  
  return(fluidSCraw(mergedFluidSCraw))
}