#' Function to create a SOM model
#'
#' This function creates a SOM (Self-Organizing Map) model. Can be plotted in various ways
#' @param fluidSCproc fluidSCproc S3 object
#' @param based_on_values values to use, defaults to "log2ExNorm"
#' @param xdim,ydim dimensions of SOM map
#' @param topo topology of SOM map
#' @param rlen the number of times the complete data set will be presented to the network
#' @param alpha learning vector, see kohonen package
#' @param toroidal join edges of map, see kohonen package
#' @param n.hood shape of neighbourhood
#' @param keep.data save data in object
#' @param NAvalues what to do with NA values
#' @param not_detected_value value to replace NA values with
#' @return returns a SOM object
#' @export
#' @details NA
#' @references http://stackoverflow.com/questions/19858729/r-package-kohonen-how-to-plot-hexagons-instead-of-circles-as-in-matlab-som-too
#' #https://github.com/geoss/som_visualization_r
#' http://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/
#' @examples
#' create_SOM_model()

create_SOM_model <- function(fluidSCproc, based_on_values = 'log2ExNorm', scale_data = T,
                             xdim = 8, ydim=6, topo=c("hexagonal","rectangular"),
                             rlen = 100, alpha = c(0.05, 0.01), toroidal = F, n.hood = c("square","circular"), keep.data = TRUE,
                             NAvalues = c("remove_assays","replace_with_not_detected"), not_detected_value = 0) {
  
  # reference: http://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/
  
  # load libraries & unload monocle and ggplot2, because otherwise error
  if("monocle" %in% (.packages()))  detach("package:monocle", unload = T) 
  if("ggplot2" %in% (.packages()))  detach("package:ggplot2", unload = T)
  
  library(kohonen)
  
  if(nargs() == 0) stop(paste0("you need to provide parameters, for more info see ?",sys.call()))
  
  topo = match.arg(topo); print(topo)
  n.hood = match.arg(n.hood); print(n.hood)
  
  
  
  # create dataframe
  normFluidCt <- fluidSCproc$data
  normMatrix <- dcast(normFluidCt, formula = Samples ~ Assays, value.var = based_on_values); rownames(normMatrix) <- normMatrix$Samples; normMatrix <- normMatrix[, -1]
  
  ## for merged fluidSCproc objects, what to do with NA values?
  
  NAvalues <- match.arg(NAvalues)
  
  # remove genes with NA values
  if(NAvalues == "remove_assays") {
    NAgenes <- colnames(normMatrix)[apply(normMatrix, 2, FUN = function(x) any(is.na(x)))]
    if(identical(NAgenes, character(0))) {
      normMatrix <- normMatrix
    }
    else normMatrix <- normMatrix[,-(which(colnames(normMatrix) %in% NAgenes))]
  } 
  # replace NA values with 0
  else if(NAvalues == "replace_with_not_detected") {
    normMatrix[is.na(normMatrix)] <- not_detected_value
  } else stop("Choose one of the 2 options to take care of NA values")
  
  
  # remove genes with zero variance!
  keepgenes <- !apply(normMatrix, MARGIN = 2, function(x) var(x,na.rm = T)) == 0
  normMatrix <- normMatrix[, keepgenes]
  
  
  # Create a training data set (rows are samples, columns are variables)
  expr_matrix <- t(normMatrix)
  
  # Change the data frame with training data to a matrix
  # Also center and scale all variables to give them equal importance during
  # the SOM training process.
  if(scale_data) expr_matrix <- scale(expr_matrix)
  data_train_matrix <- as.matrix(normMatrix)
  
  # Create the SOM Grid - you generally have to specify the size of the 
  # training grid prior to training the SOM. Hexagonal and Circular 
  # topologies are possible
  som_grid <- somgrid(xdim = xdim, ydim = ydim, topo = topo)
  
  # Finally, train the SOM, options for the number of iterations,
  # the learning rates, and the neighbourhood are available
  som_model <- som(data = data_train_matrix, 
                   grid = som_grid, 
                   rlen = rlen, 
                   alpha = alpha, 
                   keep.data = keep.data,
                   n.hood = n.hood) 
  
  return(som_model)
}
