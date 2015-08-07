
#' Function to create fluidSCraw object from fluidigm .csv output
#' 
#' This function takes as input the output table format from Fluidigm Data analysis and returns an S3 object 'fluidSCraw'
#' @param path_to_csv_table_results path to the stored fluidigm results table in .csv format (e.g. "/path/to/data.csv")
#' @param skiplines number of lines to skip, defaults to 11
#' @keywords read fluidSCraw
#' @return Object of type fluidSCraw
#' @export
#' @examples
#' read_Fluidigm()

read_Fluidigm <- function(path_to_csv_table_results, skiplines = 11) {
  
  table_results <- read.csv(file = path_to_csv_table_results, header = T, skip = skiplines)
  colnames(table_results)[[2]] <- "Samples"
  colnames(table_results)[[3]] <- "Sample_type"
  colnames(table_results)[[5]] <- "Assays"
  colnames(table_results)[[6]] <- "Assay_type"
  colnames(table_results)[[11]] <- "Melt_temp"
  
  return(fluidSCraw(table_results))
  
  
}
