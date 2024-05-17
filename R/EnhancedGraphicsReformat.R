#' Enhanced Graphics to add modifications to WikiPathways in Cy3 ####
#'
#' This function takes a file with ___ columns, and reformats it for import into Cytoscape for use with the \bold{EnhancedGraphic} package.
#'
#'
#'
#' @param input.df Input dataframe. Must have columns as follows:
#' \bold{ID:} JUN-S73
#' \bold{Symbols:} JUN
#' \bold{Sites:} S73
#' \bold{Effect:} NA
#' \bold{log2FC:} 0.5555
#'
#'
#' @author Brandon D. Murugan

EnhancedGraphicsReformat <- function(input.df){
# require(dplyr)
# require(miscTools)
# require(reshape2)
# require(tidyr)
temp <- input.df
temp[temp==""] <- NA
cardinalpoints <-   c('"northwest"', '"north"', '"northeast"', '"southwest"', '"south"',  '"southeast"', '"west"', '"east"')
final_df <- data.frame()

for (prot in unique(na.omit(temp$Symbols))){
  temp2 <- temp[grep(paste0("^",prot,"$"), temp$Symbols),]
  if (any(duplicated(temp2$ID))) next
  enhancedProt <- temp2[1,-1]
  enhancedProt$Symbols <- prot
  if (any(temp2$ID==prot)){enhancedProt$FC <- temp2[temp2$ID==prot,"FC"]}
  if (all(temp2$ID!=prot)){enhancedProt$FC <- NA}
  enhancedProt$Sites  <- NA
  temp2$info <- paste0((temp2$ID), " ", (temp2$FC))
  enhancedProt$info <- paste(unlist(temp2$info), collapse = "|")
  enhancedProt <- cbind(enhancedProt, mod1 = 'label: padding="0.05" labelsize="8" outline="false" bgColor=varA outlineWidth="3.0" background="true" outlineColor=varB label=varD position=varC')
  cnames <- colnames(enhancedProt)[1:length(enhancedProt)-1]
  enhancedProt <- cbind(enhancedProt, rep(enhancedProt["mod1"],7))
  colnames(enhancedProt) <- c(cnames, paste0("mod",seq(1, 8)))
  if (any(is.na(temp2[, "Sites"])) & nrow(temp2) > 1){temp2 <- temp2[-c(temp2$ID==prot),]}
  nsite <- 0
  for (site in na.omit(temp2$Sites)){
    nsite <- nsite + 1
    if (is.null(site)) next
    nSites <- length(na.omit(temp2$Sites))
    if (temp2[grep(site, temp2$Sites),"FC"] < 0){varA <- '"#4292C6"'}
    if (temp2[grep(site, temp2$Sites),"FC"] > 0){varA <- '"#EF3B2C"'}
    varB <- "#00FFFFFF"
    varC <- cardinalpoints[nsite]
    varD <- paste0("'",site,"'")
    modstring <- 'label: padding="0.05" labelsize="8" outline="false" bgColor=varA outlineWidth="3.0" background="true" outlineColor=varB label=varD position=varC'
      modstring <- sub("varA",varA,modstring)
      modstring <- sub("varB",varB,modstring)
      modstring <- sub("varC",varC,modstring)
      modstring <- sub("varD",varD,modstring)
      modstring <- gsub("\'","\"",modstring)
      enhancedProt[,grep(paste0("mod",nsite), colnames(enhancedProt))] <- modstring
  }
  final_df <- rbind(final_df, enhancedProt)
}
final_df[apply(final_df, 1:2, function(i) grepl('\\varC', i))] <- ""

return(final_df)
}
