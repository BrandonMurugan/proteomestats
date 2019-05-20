#' Valid value filter for label-free proteomic LC-MS experiments
#'
#' This function takes a file with protein names 1st column) and quantitation columns, and applies a user defined valid value filter with per-group constraints.
#' The filter will assess each group for the minimum number of valid values, and return the rows where this filter passes in \bold{all} groups
#' The output file is the expression matrix with rows that do \bold{not} pass the filter removed
#'
#'
#' @param input.df Input dataframe. Must have one column with protein/gene identifier and rest of columns with expression values for all samples
#' @param experimental.groups List of experimental groups with replicate suffix excluded (e.g.: c("Control", "Treatment"), \bold{not} c("Control_1", "Control_2", etc). This list must match at least part of the intensity column names in the expression matrix)
#' @param valid.values Number of valid values to filter for per-group
#' @param zero.as.na Returns dataframe with zeros converted to NAs (Default = T)
#'
#'
#' @author Brandon D. Murugan

# input.df <- test
# valid.values <- 2
# experimental.groups <- experiments2$experiments

filter_valid_values <- function(input.df, experimental.groups, valid.values = 2, zero.as.na=T){

  require(dplyr)
  require(splitstackshape)
  require(here)
  require("miscTools")
  require(preprocessCore)
  require(XML)
  require(dplyr)
  require(reshape2)
  require(tidyr)

  input.df[input.df==0] <- NA
  temporary <- data.frame(NA)
  colnames(temporary) <- "NA"
  test2 <- cbind(NA,input.df[1,])
  test2[1,] <- NA

  for (i in rownames(input.df)){
    temporary <- data.frame(NA)
    colnames(temporary) <- "NA"
    for (j in experimental.groups){
      # if ((sum(input.df[i,grep(pattern = j, x = colnames(input.df))] != 0) > (valid.values-1))) {
      if (sum(is.na(test[i,grep(pattern = j, x = colnames(test))])) <= (length(test[i,grep(pattern = j, x = colnames(test))]) - valid.values)) {
        temporary <- cbind(temporary,input.df[i,grep(pattern = j, x = colnames(input.df))])}
    }
    if (dim(temporary)[2] == dim(input.df)[2]){
      temporary <- cbind(IDname=input.df[i,1], temporary)
      colnames(temporary)[1] <- colnames(input.df)[1]
      temporary <- temporary[names(test2)]
      test2 <- rbind(test2,temporary)
    }
  }

  test2[,"NA"] <- NULL
  test2 <- test2[-1,]
  test3 <- test2
# removes rows where there are ALL NAs
  for (i in rownames(test2)){
    currentrow <- as.numeric(i)
    if (sum(is.na(test2[currentrow,-1])) == (dim(input.df)[2]-1)){
      test2 <- test2[-currentrow,]
    }
  }

  if (zero.as.na == F){
    test2[test2==NA] <- 0
  }
return(test2)
}
