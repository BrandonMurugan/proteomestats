#' Protein ID Filter
#'
#' This function takes a file with protein IDs and returns the 1st \bold{sp} (SwissProt) identifier. The function will return the 1st \bold{tr} (trembl) entry if no sp entries are found.
#'
#' The output file is the input file with an additional \bold{user defined} ProteinID column (Default = \bold{ProteinID_main})
#'
#'
#' @param input.df Input dataframe (e.g.: MaxQuant proteinGroups.txt file) (environmental variable)}
#' @param ProteinID.column Protein ID column to be filtered through to find 1st instance of an \bold{sp} identifier (string)
#' @param delimiter Delimiter within \bold{ProteinID.column} (string)
#' @param new.ID.column Name of new ID column. Default is \bold{ProteinID_main} (string)
#'
#'
#' @author Brandon D. Murugan

# ProteinID.column <- "Protein.IDs"
# input.df <- proteinGroups
# delimiter <- ";"
# new.ID.column <- "ProteinID_main"

# Extract string following [sp|]. Will take whole string if [sp|] not found. keep 1st protein
sp.select <- function(input.df, ProteinID.column, new.ID.column = "ProteinID_main", delimiter){
  library(splitstackshape)
  library(dplyr)
  for (i in rownames(input.df)){
    temp <- grep("sp\\|", unlist(stringr::str_split(input.df[i,grep(paste0("^",ProteinID.column,"$"), colnames(input.df))], delimiter)))
    if (length(temp) == 0){
      temp <- 1
      if (grepl("Uncharacterized", unlist(stringr::str_split(input.df[i,grep(paste0("^",ProteinID.column,"$"), colnames(input.df))], delimiter))[temp]))
        temp <- temp+1
    }
    if (length(temp) != 0)
      temp <- temp[1]
    input.df[i,new.ID.column] <- unlist(stringr::str_split(input.df[i,grep(paste0("^Protein.IDs$"), colnames(input.df))], delimiter))[temp]
  }
  return(input.df)
}
