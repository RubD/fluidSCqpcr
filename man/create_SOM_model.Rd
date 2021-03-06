% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/create_SOM_model_function.R
\name{create_SOM_model}
\alias{create_SOM_model}
\title{Function to create a SOM model}
\usage{
create_SOM_model(fluidSCproc, based_on_values = "log2ExNorm",
  scale_data = T, xdim = 8, ydim = 6, topo = c("hexagonal",
  "rectangular"), rlen = 100, alpha = c(0.05, 0.01), toroidal = F,
  n.hood = c("square", "circular"), keep.data = TRUE,
  NAvalues = c("remove_assays", "replace_with_not_detected"),
  not_detected_value = 0)
}
\arguments{
\item{fluidSCproc}{fluidSCproc S3 object}

\item{based_on_values}{values to use, defaults to "log2ExNorm"}

\item{xdim,ydim}{dimensions of SOM map}

\item{topo}{topology of SOM map}

\item{rlen}{the number of times the complete data set will be presented to the network}

\item{alpha}{learning vector, see kohonen package}

\item{toroidal}{join edges of map, see kohonen package}

\item{n.hood}{shape of neighbourhood}

\item{keep.data}{save data in object}

\item{NAvalues}{what to do with NA values}

\item{not_detected_value}{value to replace NA values with}
}
\value{
returns a SOM object
}
\description{
This function creates a SOM (Self-Organizing Map) model. Can be plotted in various ways
}
\details{
NA
}
\examples{
create_SOM_model()
}
\references{
http://stackoverflow.com/questions/19858729/r-package-kohonen-how-to-plot-hexagons-instead-of-circles-as-in-matlab-som-too
#https://github.com/geoss/som_visualization_r
http://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/
}

