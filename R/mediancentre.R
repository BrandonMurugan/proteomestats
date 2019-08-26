#' Mediancentre function
#'
#' This function reads in a dataframe with protein ID column as the 1st column
#' 
#' median centre normalisation:
#' data distributions are adjusted such that each individual sample's medians are aligned to the average of medians of all the samples 
#'
#'   
#'
#' 


mediancentre <- function(x){
  library(dplyr)
  library(splitstackshape)
  library("miscTools")
  library(preprocessCore)
  X[x == 0] <- NA
  ProteinQuant.norm <- x[,-1]
  ColumnMedians <- colMedians(ProteinQuant.norm, na.rm = T)
  ColumnMedians <- as.data.frame(ColumnMedians)
  meanMeans <- colMeans(ColumnMedians)
  normFactor <- meanMeans/ColumnMedians
  normFactor <- t(normFactor)
  ProteinQuant.norm <- data.frame(mapply("*", ProteinQuant.norm, normFactor))
  ProteinQuant.norm <- cbind(Protein.ID = x[,1], ProteinQuant.norm)
  return(ProteinQuant.norm)
}
