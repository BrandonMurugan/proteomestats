#' Imports and preprocesses data files for label-free proteomic LC-MS experiments
#'
#' This function imports MaxQuant proteingroups.txt, peptides.txt, Phospho (STY)Sites.txt (if it exists), and mqpar.xml files.
#' Proteingroups.txt file is then filtered and normalised according to chosen quantitation and normlisation types
#'
#' @param txtpath MaxQuant txt folder
#' @param mqparPath MaxQuant mqpar.xml file path
#' @param quant_type Quanitiation type to select. Possible choices are "Intensity" (default), or "LFQ".
#' @param norm_type Type of normalisation to perform. Possible choices are "MedianCentre" (default), "Quantile", or "None".
#' @param valid_values Numeric. Removes entire protein if number of valid (non-zero) values are less than stated amount for \bold{any} group
#' @param imputation Logical. If \code{TRUE}, will impute with lowest value of all intensities found in dataset. Default is FALSE
#'
#'
#' @author Brandon D. Murugan

data_preprocess <- function(txtpath, mqparfile, quant_type = "Intensity", norm_type = "MedianCentre", valid_values = 3, imputation = F){
  library(grid)
  library(XML)
  library(xml2)
  library(dplyr)
  library(reshape2)
  library(tidyr)
  library(splitstackshape)
  library(here)
  library("miscTools")
  library(preprocessCore)

  direct <- txtpath

  # direct <- "C:/Users/Brandon/Documents/Data Analysis/Kim_reprocess/H_timecourse"
  # load(paste0(direct,"/images/imported_Data.RData"))

  # Rip info from mqpar file (also moved into text folder; less work on user's part) ####
          # mqpar <- xmlTreeParse(file = paste0(direct,"//mqpar.xml"))
  mqpar <- xmlTreeParse(file = mqparfile)
  # mqpar <- read_xml(paste0(direct,"\\mqpar.xml"))
  topMqpar <- xmlRoot(mqpar)
  # topMqpar <- xml_root(mqpar)
  topMqpar <- xmlApply(topMqpar, function(x) xmlSApply(x, xmlValue))
  # topMqpar <- xml2::(topMqpar)
  # mqpar_df <- data.frame(t(topMqpar), row.names = NULL)
  mqpar_df <- data.frame(t(topMqpar), row.names = NULL)

  expDesign <- data.frame(mqpar_df$filePaths)
  expDesign$experiments <- data.frame(mqpar_df$experiments)
  experiments <<- expDesign$experiments
  Intensity.experiments <- paste("Intensity",experiments$experiments,sep=".")
  expType <- data.frame(mqpar_df$restrictMods)

  # set experiment title based on raw file folder name (for report filename)
  parentFolder <- as.character(expDesign$filePaths[1])
  parentFolder <- strsplit(parentFolder, split = "\\\\")
  parentFolder <- unlist(parentFolder)
  parentFolder <- parentFolder[length(parentFolder)-1]
  versionNumber <- data.frame(mqpar_df$maxQuantVersion)
  # dir.create("/temp", showWarnings = T)   #made in python
  save(versionNumber, file = paste0(direct,"/temp/versionNumber"), compress = T)

  # Import txt files ####
  peptides <- read.delim(file = paste0(direct,"//peptides.txt"))
  proteinGroups <- read.delim(file = paste0(direct,"//proteinGroups.txt"))
  summary <- read.delim(file = paste0(direct,"//summary.txt"))
  if (grepl("Phos", expType[1,]) == TRUE){
    PhosphoSTYsites <- read.delim(file = paste0(direct,"//Phospho (STY)Sites.txt"))
  }
  # dir.create("/images", showWarnings = T)   #made in python
  # save.image(file = paste0(direct,"/images/imported_Data"), compress = T)
  # save(list = ls(all.names = TRUE), file = paste0(direct,"/images/imported_Data.RData"))

  # txtpath <- "C:/Users/Brandon/Documents/Data Analysis/Kim_reprocess/H_timecourse"

  # load(paste0(direct,"/images/imported_Data.RData"))

  # Remove Contaminants, Reverse hits ####
  # proteinGroups <- read.delim("~/example_txt/proteinGroups.txt", sep = "\t")
  PGs <- proteinGroups
  PGs <- PGs[PGs$Potential.contaminant != "+" & PGs$Reverse != "+",]

  # Make Quantitation file (choose LFQ/Intensities from python GUI**) ####
  uniprot.IDs <- data.frame(Protein.ID = PGs$Majority.protein.IDs, stringsAsFactors = F)
  uniprot.IDs$Protein.ID <- as.character(uniprot.IDs$Protein.ID)
  uniprot.IDs <- cSplit(indt = uniprot.IDs, splitCols = "Protein.ID", sep = ";")
  uniprot.IDs <- data.frame(uniprot.IDs[,1])
  colnames(uniprot.IDs) <- "Protein.ID"

  if (quant_type == "Intensity"){
    intensity.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "Intensity.")])
    intensity.PG$Protein.ID <- as.character(intensity.PG$Protein.ID)
    # Make quant file with common name despite LFQ/Intensity (makes downstream tasks easier)
    ProteinQuant <- intensity.PG
  }

  if (quant_type == "LFQ"){
    LFQ.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "LFQ.")])
    LFQ.PG$Protein.ID <- as.character(LFQ.PG$Protein.ID)
    # Make quant file with common name despite LFQ/Intensity (makes downstream tasks easier)
    ProteinQuant <- LFQ.PG
  }
  # Normalise data (usually for Intensity data, Choose in python GUI**) ####
  if (norm_type == "MedianCentre"){
    ProteinQuant[ProteinQuant == 0] <- NA
    ProteinQuant.norm <- ProteinQuant[,-1]
    ColumnMedians <- colMedians(ProteinQuant.norm, na.rm = T)
    ColumnMedians <- as.data.frame(ColumnMedians)
    meanMeans <- colMeans(ColumnMedians)
    normFactor <- meanMeans/ColumnMedians
    normFactor <- t(normFactor)
    ProteinQuant.norm <- data.frame(mapply("*", ProteinQuant.norm, normFactor))
    ProteinQuant.norm <- cbind(Protein.ID = ProteinQuant[,1], ProteinQuant.norm)
  }

  if (norm_type == "None"){
    ProteinQuant.norm <- ProteinQuant
  }

  if (norm_type == "Quantile"){
    ProteinQuant[ProteinQuant == 0] <- NA
    ProteinQuant.norm <- normalize.quantiles(as.matrix(ProteinQuant[,-1]))
    ProteinQuant.norm <- as.data.frame(ProteinQuant.norm)
    ProteinQuant.norm <- cbind(ProteinQuant[,1], ProteinQuant.norm)
    # colnames(ProteinQuant.norm) <- colnames(ProteinQuant)
  }

  load("ProteinQuant.norm.rda")

  # # VALID VALUE HANDLING
  #   ProteinQuant.norm <- ProteinQuant.norm[complete.cases(ProteinQuant.norm) ,]

  # IMPUTATION IN PROGRESS
  if (imputation == T){
    for (i in colnames(ProteinQuant.norm[,-1])){
      ProteinQuant.norm[,i] <- tidyr::replace_na(ProteinQuant.norm[,i], min(ProteinQuant.norm[,i], na.rm = T))
    }
    }



}
