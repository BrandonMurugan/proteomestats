#' Protein ID Filter
#'
#' This function takes a file with protein IDs and returns the 1st \bold{sp} (SwissProt) identifier. The function will return the 1st \bold{tr} (trembl) entry if no sp entries are found.
#'
#' The output file is the input file with an additional \italic{user defined} ProteinID column (Default = \bold{ProteinID_main})
#'
#'
#' @param input.df Input dataframe (e.g.: MaxQuant proteinGroups.txt file) \italic{(environmental variable)}
#' @param ProteinID.column Protein ID column to be filtered through to find 1st instance of an \bold{sp} identifier \italic{(string)}
#' @param delimiter Delimiter within \bold{ProteinID.column} \italic{(string)}
#' @param new.ID.column Name of new ID column. Default is \bold{ProteinID_main} \italic{(string)}
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
  sp_proteins <- data.frame(sp_proteins=substring(input.df[,grep(paste0("^",ProteinID.column,"$"), colnames(input.df))],regexpr("sp\\|", input.df[,grep(paste0("^",ProteinID.column,"$"), colnames(input.df))])))

  sp_proteins <- cSplit(indt = sp_proteins, splitCols = "sp_proteins", sep = delimiter)
  sp_proteins <- sp_proteins[,1:grep("sp_proteins_01", colnames(sp_proteins))]
  colnames(sp_proteins) <- sub("_01","",colnames(sp_proteins))
  sp_proteins$sp_proteins <- as.character(sp_proteins$sp_proteins)
  colnames(sp_proteins) <- new.ID.column
  input.df <- cbind(sp_proteins, input.df, deparse.level = 1)
  input.df <- data.frame(input.df)
  return(input.df)
}
