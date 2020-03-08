#' Gene Name Extractor
#'
#' This function Extracts the Gene Name and Uniprot Identifiers from the Fasta.header column of the ProteinGroups file.
#' Will also extract the appropriate Uniprot ID.
#' Can select first SwissProt entry in protein group if required
#'
#'
#' @param inputDF Input MQ proteomic dataframe with Fasta.headers column
#' @param SwissProt.only Will select first SwissProt entry in ProteinGroups. (Default = TRUE)
#'
#'
#' @author Brandon D. Murugan

gene.name.extractor <- function(inputDF, SwissProt.only = T){
  library(EnsDb.Hsapiens.v86)
  # library(SetRank)
  # source('~/Data Analysis/proteomestats/R/createIDconverter2.R')
  inputDF$Fasta.headers <- sub(";.*$","",inputDF$Fasta.headers)
  rownames(inputDF) <- inputDF$seqid <- seq(1,length(inputDF$Fasta.headers))
  unkGene.names <- inputDF
  unkGene.names <- data.frame(GeneNames = unkGene.names$Fasta.headers, seqid=unkGene.names$seqid, stringsAsFactors = F)
  unkGene.names$GeneNames <- as.character(unkGene.names$GeneNames)

  GN <- unkGene.names[grep("GN=[[:alnum:]]+",unkGene.names$GeneNames),]
  rest <- setdiff(unkGene.names$seqid,GN$seqid)
  nonGN <- unkGene.names[unkGene.names$seqid %in% rest,]

  GN$GeneNames <- sub("(^.*.GN=)","",GN$GeneNames)
  GN$GeneNames <- sub("( PE.*.$)","",GN$GeneNames)
  nonGN$GeneNames <- sub("^sp\\|.*\\|","",nonGN$GeneNames)
  nonGN$GeneNames <- sub("^tr\\|.*\\|","",nonGN$GeneNames)
  nonGN$GeneNames <- sub("_HU.*$","",nonGN$GeneNames)
  uniprot2symbol <- proteomestats::createIDconverter2(annotationPackageName = "EnsDb.Hsapiens.v86", from = "UNIPROTID", to = "SYMBOL")
  nonGN$GeneNames <- uniprot2symbol(nonGN$GeneNames, na.rm=FALSE, drop.ambiguous=TRUE)
  GeneNames_2 <- rbind(GN,nonGN)

  temp <- inputDF[inputDF$seqid %in% GeneNames_2$seqid,]
  inputDF_2 <- inputDF[!inputDF$seqid %in% GeneNames_2$seqid,]
  temp$Gene.names <- NULL

  temp <- dplyr::inner_join(GeneNames_2,temp, by = "seqid")
  colnames(temp)[1] <- "Gene.names"
  temp <- temp[,colnames(inputDF_2)]
  inputDF_2 <- rbind(inputDF_2,temp)
  inputDF_2 <- inputDF_2[!is.na(inputDF_2$Gene.names),]
  inputDF_2$Gene.names <- sub(";.*$","",inputDF_2$Gene.names)
  inputDF_2$seqid <- NULL
  inputDF_2$uniprot <- sub("^(tr|sp)\\|","",inputDF_2$Fasta.headers)
  inputDF_2$uniprot <- sub("\\|.*$","",inputDF_2$uniprot)
  if (!is.null(inputDF_2$id)) inputDF_2$id <- as.numeric(inputDF_2$id)
  return(inputDF_2)
}
