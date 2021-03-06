% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Kmeans_cluster_estimation_functions.R
\name{sparseKmeansCBI}
\alias{sparseKmeansCBI}
\title{Helper function for Kmeans_cluster_statistics()}
\usage{
sparseKmeansCBI(inputMatrix, sparseK, wbound = 2)
}
\arguments{
\item{sparseK}{number of clusters for sparse Kmeans clustering}

\item{wbound}{weight for kmeans clustering bound to features (genes)}

\item{inputMatriw}{matrix with scores}
}
\value{
returns list that can be used with clusterboot function
}
\description{
Helper function for Kmeans_cluster_statistics()
}

