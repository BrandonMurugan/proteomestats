#' Hello, world!
#'
#' This is an example function named 'hello'
#' which prints 'Hello, world!'.
#'
#' You can learn more about package authoring with RStudio at:
#'
#'   http://r-pkgs.had.co.nz/
#'
#' Some useful keyboard shortcuts for package authoring:
#'
#'   Build and Reload Package:  'Ctrl + Shift + B'
#'   Check Package:             'Ctrl + Shift + E'
#'   Test Package:              'Ctrl + Shift + T'

mediancentre <- function(x = ProteinQuant){
  library(dplyr)
  library(splitstackshape)
  library(here)
  library("miscTools")
  library(preprocessCore)
  ProteinQuant[ProteinQuant == 0] <- NA
  ProteinQuant.norm <- ProteinQuant[,-1]
  ColumnMedians <- colMedians(ProteinQuant.norm, na.rm = T)
  ColumnMedians <- as.data.frame(ColumnMedians)
  meanMeans <- colMeans(ColumnMedians)
  normFactor <- meanMeans/ColumnMedians
  normFactor <- t(normFactor)
  ProteinQuant.norm <- data.frame(mapply("*", ProteinQuant.norm, normFactor))
  ProteinQuant.norm <- cbind(Protein.ID = ProteinQuant[,1], ProteinQuant.norm)
  return(ProteinQuant.norm)
}
