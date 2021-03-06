% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/scatterplot_SC_function.R
\name{scatterplot_SC}
\alias{scatterplot_SC}
\title{Function to create scatterplot between 2 assays}
\usage{
scatterplot_SC(fluidSCproc, based_on_values = "log2Ex", gene1, gene2,
  my_formula = "y ~ x", se_display = F, sample_display = F,
  cor_use = "complete.obs", cor_method = "pearson", smooth_method = "lm",
  remove_not_detected = F, not_detected_value = 0)
}
\arguments{
\item{fluidSCproc}{fluidSCproc S3 object}

\item{based_on_values}{values to use, defaults to "log2Ex"}

\item{gene1}{first selected gene, X-axis}

\item{gene2}{second selected gene, Y-axis}

\item{my_formula}{formula to create smoothed line}

\item{se_display}{display standard error range around smoothed line}

\item{sample_display}{add sample names on plot}

\item{cor_use}{data to use to calculate correlation score, defaults to "complete.obs"}

\item{cor_method}{method for correlation score algorithm, defaults to "pearson"}

\item{smooth_method}{method to create smoothed line}

\item{remove_not_detected}{boolean to remove or keeop not detected values}

\item{not_detected_value}{define what values are considered not detected}
}
\value{
returns a scatterplot with smoothed line and (always linear) correlation score
}
\description{
Function creates scatterplot and shows correlation score between 2 assays (genes)
}
\details{
NA
}
\examples{
scatterplot_SC()
}

