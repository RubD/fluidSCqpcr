% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/decide_CallMeltTemp_function.R
\name{decide_Call_meltTemp}
\alias{decide_Call_meltTemp}
\title{Function to call quality of values based on list with melting temparatures for all assays}
\usage{
decide_Call_meltTemp(fluidSCraw, meltTempRange = c(2, 2),
  named_listMeltTemp = NULL)
}
\arguments{
\item{fluidSCraw}{fluidSCraw S3 object}

\item{meltTempRange}{range around melt temperature peak}

\item{named_listMeltTemp}{list with melting temperatures and names that correspond to the assays}
}
\value{
fluidSCraw S3 object were new Call column has been added or replaced
}
\description{
This function takes as input a list with melt temperatures for all assays and calculates whether data points are within range
}
\examples{
decide_Call_meltTemp()
}
\keyword{fluidSCraw}
\keyword{melt}
\keyword{temperature}

