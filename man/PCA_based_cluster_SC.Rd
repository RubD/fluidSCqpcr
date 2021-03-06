% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/PCA_based_cluster_SC_function.R
\name{PCA_based_cluster_SC}
\alias{PCA_based_cluster_SC}
\title{Function to cluster single cell data}
\usage{
PCA_based_cluster_SC(fluidSCproc, based_on_values = "log2ExNorm",
  scaleData = T, PCAscaled = F, clustMethod = c("Kmeans", "Hierarchical",
  "Correlation"), firstPC = "PC1", secondPC = "PC2", nrClust = 2,
  hclustMethod = c("ward.D", "ward.D2", "single", "complete", "average",
  "centroid"), cluster_column_name = "PCAclust", return_sign_loadMatrix = F,
  NAvalues = c("remove_assays", "replace_with_not_detected"),
  not_detected_value = 0, ...)
}
\arguments{
\item{fluidSCproc}{fluidSCproc S3 object~}

\item{based_on_values}{values to use, defaults to "log2Ex"}

\item{scaleData}{boolean to scale data before clustering, default = F}

\item{clustMethod}{choose method for clustering samples}

\item{nrClust}{number of clusters you want to detect with Kmeans}

\item{hclustMethod}{choose method to cluster distance matrix}

\item{cluster_column_name}{name of column that will be created with the clustering results}

\item{return_sign_loadMatrix}{boolean to return loadings multiplied by their sign for all PCs, default to FALSE}

\item{NAvalues}{how to handle NA values, remove or replace with not_detected_value}

\item{not_detected_value}{value to replace NA values}

\item{...}{additional parameters for distance methods}

\item{firstPC/secondPC}{principal components to use}
}
\value{
returns a fluidSCproc S3 object appended with the cluster results
}
\description{
This function merges several clustering methods and appends the result to the fluidSCproc object
}
\details{
NA
}
\examples{
PCA_based_cluster_SC()
}

