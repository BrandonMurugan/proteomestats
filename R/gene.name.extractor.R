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
  # library(EnsDb.Hsapiens.v86)
  data(Uniprot2Genename_swissprot_20200224)
  data(NCBI_200210)
  # library(SetRank)
  # source('~/Data Analysis/proteomestats/R/createIDconverter2.R')
  inputDF <- sp.select(input.df = inputDF, ProteinID.column = "Fasta.headers", new.ID.column = "Protein.ID_spFilter", delimiter = ";")
  rownames(inputDF) <- inputDF$seqid <- seq(1,length(inputDF$Fasta.headers))
  # inputDF$Fasta.headers <- sub(";.*$","",inputDF$Fasta.headers)
  unkGene.names <- inputDF
  unkGene.names <- data.frame(GeneNames = unkGene.names$Protein.ID_spFilter, GeneName_MQ = unkGene.names$Gene.names, seqid=unkGene.names$seqid, stringsAsFactors = F)
  unkGene.names$GeneNames <- as.character(unkGene.names$GeneNames)
  unkGene.names$GeneName_MQ <- as.character(unkGene.names$GeneName_MQ)

  # GN from Fasta Header
  GN <- unkGene.names[grep("PE=[[:alnum:]]+",unkGene.names$GeneNames),]  #Fasta header until PE at least
  GN <- GN[grep("GN=[[:alnum:]]+",GN$GeneNames),]
  GN$GeneNames <- sub("(^.*.GN=)","",GN$GeneNames)
  GN$GeneNames <- sub("( PE.*.$)","",GN$GeneNames)

  # from MQ genename column
  rest <- setdiff(unkGene.names$seqid,GN$seqid)
  nonGN <- unkGene.names[unkGene.names$seqid %in% rest,]
  nonGN$GeneName_MQ <- sub(";.*$","",nonGN$GeneName_MQ)
  nonGN$GeneNames <- nonGN$GeneName_MQ

  GeneNames_2 <- rbind(GN, nonGN)
  # temp <- inputDF[inputDF$seqid %in% GeneNames_2$seqid,]
  # inputDF_2 <- inputDF[!inputDF$seqid %in% GeneNames_2$seqid,]
  temp <- dplyr::inner_join(GeneNames_2,temp, by = "seqid")
  temp$Gene.names <- NULL
  colnames(temp)[1] <- "Gene.names"
  inputDF_2 <- temp
  inputDF_2 <- inputDF_2[!is.na(inputDF_2$Gene.names),]
  # inputDF_2$Gene.names <- sub(";.*$","",inputDF_2$Gene.names)
  inputDF_2$seqid <- NULL
  inputDF_2$uniprot <- sub("^(tr|sp)\\|","",inputDF_2$Protein.ID_spFilter)
  inputDF_2$uniprot <- sub("\\|.*$","",inputDF_2$uniprot)
  if (!is.null(inputDF_2$id)) inputDF_2$id <- as.numeric(inputDF_2$id)

  # Sort out blanks
  inputDF_2[grep("^$",inputDF_2$Gene.names),"Gene.names"] <- inputDF_2[grep("^$",inputDF_2$Gene.names),"uniprot"]
  inputDF_2[grep("^$",inputDF_2$Gene.names),"GeneName_MQ"] <- inputDF_2[grep("^$",inputDF_2$Gene.names),"uniprot"]

  # Adding reviewed uniprot IDs from genenames
  inputDF_2 <- dplyr::left_join(inputDF_2, Uniprot2Genename_swissprot_20200224[,c("Gene.names","Entry")], "Gene.names")
  inputDF_2$Entry <- as.character(inputDF_2$Entry)
  inputDF_2$Entry[is.na(inputDF_2$Entry)] <- inputDF_2$uniprot[is.na(inputDF_2$Entry)]


  # Add ENTREZID from NCBI200210
  colnames(NCBI_200210)[colnames(NCBI_200210) == "Symbol"] <- "Gene.names"
  inputDF_2 <- dplyr::left_join(inputDF_2, NCBI_200210[,c("Gene.names", "GeneID")], "Gene.names")

  # uniprot2entrez <- proteomestats::createIDconverter2(annotationPackageName = "EnsDb.Hsapiens.v86", from = "UNIPROTID", to = "ENTREZID")
  # nonGN$GeneNames <- uniprot2symbol(nonGN$GeneNames, na.rm=FALSE, drop.ambiguous=TRUE)

  colnames(inputDF_2)[colnames(inputDF_2)=="Entry"] <- "UniProt_sp"
  return(inputDF_2)
}
