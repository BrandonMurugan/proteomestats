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
  library(proteomestats)
  data("Uniprot2Genename_swissprot_20200313")
  # data(NCBI_assembly_200313)
  data("NCBI_genes_200210_unique")
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
  temp <- dplyr::inner_join(GeneNames_2, inputDF, by = "seqid")
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

  inputDF_2 <- dplyr::left_join(inputDF_2, Uniprot2Genename_swissprot_20200313[,c("Gene.names","Entry")], "Gene.names")

  # Adding reviewed uniprot IDs from genenames
  inputDF_2 <- dplyr::left_join(inputDF_2, Uniprot2Genename_swissprot_20200313[,c("Gene.names","Entry")], "Gene.names")
  inputDF_2$Entry <- as.character(inputDF_2$Entry)
  inputDF_2$Entry[is.na(inputDF_2$Entry)] <- inputDF_2$uniprot[is.na(inputDF_2$Entry)]

  # # Add ENTREZID from NCBI_200210
  # colnames(NCBI_genes_200210_unique)[colnames(NCBI_genes_200210_unique) == "Symbol"] <- "Gene.names"
  # inputDF_3 <- dplyr::left_join(inputDF_2, NCBI_genes_200210_unique[,c("Gene.names", "GeneID")], "Gene.names")
  # noEntrez <- inputDF_3[is.na(inputDF_3$GeneID),]
  # inputDF_3 <- inputDF_3[!is.na(inputDF_3$GeneID),]
  #
  # symbol2entrez <- proteomestats::createIDconverter2(annotationPackageName = "EnsDb.Hsapiens.v86", from = "SYMBOL", to = "ENTREZID")
  # noEntrez$GeneID <- symbol2entrez(noEntrez$Gene.names, na.rm=FALSE, drop.ambiguous=TRUE)
  # noEntrez_2 <- noEntrez[is.na(noEntrez$GeneID),]
  # noEntrez <- noEntrez[!is.na(noEntrez$GeneID),]
  #
  #
  # noEntrez_2$GeneID <- NULL
  # # Add ENTREZID from uniprot_mapping
  # noEntrez_2 <- dplyr::left_join(noEntrez_2, Uniprot2Genename_swissprot_20200313[,c("Entry", "Cross.reference..GeneID.")], "Entry")
  # colnames(noEntrez_2)[colnames(noEntrez_2) == "Cross.reference..GeneID."] <- "GeneID"
  # noEntrez_2$GeneID <- as.character(noEntrez_2$GeneID)
  # noEntrez_2$GeneID <- sub(";","",noEntrez_2$GeneID)
  # noEntrez_2$GeneID <- sub("^$",NA,noEntrez_2$GeneID)
  # noEntrez_3 <- noEntrez_2[is.na(noEntrez_2$GeneID),]
  # noEntrez_2 <- noEntrez_2[!is.na(noEntrez_2$GeneID),]
  #


  # colnames(noEntrez)[colnames(noEntrez)=="Aliases"] <- "Gene.names"
  # colnames(noEntrez_2)[colnames(noEntrez_2)=="Aliases"] <- "Gene.names"
  #
  # inputDF_4 <- rbind(inputDF_3, noEntrez, noEntrez_2)



  colnames(inputDF_2)[colnames(inputDF_2)=="Entry"] <- "UniProt_sp"
  return(inputDF_2)
}
