% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/tSNE_SC_function.R
\name{tSNE_SC}
\alias{tSNE_SC}
\title{Function to create a t-SNE map from single cell data}
\usage{
tSNE_SC(fluidSCproc, based_on_values = "log2Ex", scaleData = F,
  perplexity = 5, clusterColumn = NULL, NAvalues = c("remove_assays",
  "replace_with_not_detected"), not_detected_value = 0)
}
\arguments{
\item{fluidSCproc}{fluidSCproc S3 object~}

\item{based_on_values}{values to use, defaults to "log2Ex"}

\item{scaleData}{boolean to scale data before tSNE clustering, default = F}

\item{perplexity}{measure for information about the Shannon entropy, see reference}

\item{clusterColumn}{name of column with cluster information}

\item{NAvalues}{how to process NA values, remove or replace with not_detected_value (default = 0)}

\item{not_detected_value}{value to replace NA values with}
}
\value{
returns a t-SNE plot
}
\description{
This function creates a heatmap?
}
\details{
NA
}
\examples{
tSNE_SC()
}
\references{
http://lvdmaaten.github.io/tsne/
}

