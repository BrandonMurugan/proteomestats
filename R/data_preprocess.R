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

data_preprocess <- function(txtfolder = txtpath, mqparPath = mqpath, Quantitation = QuantType){
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

  # Rip info from mqpar file (also moved into text folder; less work on user's part) ####
          # mqpar <- xmlTreeParse(file = paste0(direct,"//mqpar.xml"))
  mqpar <- xmlTreeParse(file = mqpath)
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

  load(paste0(direct,"/images/imported_Data.RData"))

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

  if (QuantType == "Intensity"){
    intensity.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "Intensity.")])
    intensity.PG$Protein.ID <- as.character(intensity.PG$Protein.ID)
    # Make quant file with common name despite LFQ/Intensity (makes downstream tasks easier)
    ProteinQuant <- intensity.PG
  }

  if (QuantType == "LFQ"){
    LFQ.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "LFQ.")])
    LFQ.PG$Protein.ID <- as.character(LFQ.PG$Protein.ID)
    # Make quant file with common name despite LFQ/Intensity (makes downstream tasks easier)
    ProteinQuant <- LFQ.PG
  }
}
