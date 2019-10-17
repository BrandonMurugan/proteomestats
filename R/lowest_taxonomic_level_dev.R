#' Lowest Taxonomic Level for each OTU
#'
#' This function takes in a phyloseq object and adds a column in the tax_table with the lowest taxonomic level for each OTU/ASV.
#'
#'
#' @param phyloseq Input phyloseq object
# @param combine.genus.species Combine Genus and species (Default = T)
# @param disambiguation Keep only the first entry for multiple taxa possibilites (Default = F)
#'
#'
#' @author Brandon D. Murugan


lowest.tax.level <- function(phyloseq.file)
  # , combine.genus.species = T, disambiguation = F){
  # load("../Sim_16S/ps_uchoose.rdata")

  library(phyloseq)
  library("stringr")

  # replace any "__" following taxa type in taxa entries (e.g.: k__Bacteria,	p__Firmicutes)
  phyloseq.file@tax_table@.Data <- apply(phyloseq.file@tax_table@.Data, MARGIN = 2, FUN = function(x){gsub(".*__","",x)})
  # phyloseq.file@tax_table@.Data <- apply(phyloseq.file@tax_table@.Data[],2, as.character)
  phyloseq.file@tax_table@.Data[,][phyloseq.file@tax_table@.Data[,] == "NA"] <- ""
  phyloseq.file@tax_table@.Data <- DescTools::Rev(phyloseq.file@tax_table@.Data, margin = 2)
  phyloseq.file@tax_table@.Data[is.na(phyloseq.file@tax_table@.Data)] <- ""

  tax_fin <- NULL
  flag <-  TRUE
  for (i in 1:dim(phyloseq.file@tax_table@.Data)[1]){
    for (j in 1:dim(phyloseq.file@tax_table@.Data)[2]){
      if((phyloseq.file@tax_table@.Data[i,j] != "") & (flag == TRUE)){
        if(colnames(phyloseq.file@tax_table@.Data)[j] == "Species"){
          temp <- paste0(strtrim(phyloseq.file@tax_table@.Data[i,j+1],1), ".", phyloseq.file@tax_table@.Data[i,j])
        }
        else{temp <- phyloseq.file@tax_table@.Data[i,j]}
        tax_fin <- c(tax_fin,temp)
        flag <- FALSE
      }
    }
    flag <-  TRUE
  }

  phyloseq.file@tax_table@.Data <- DescTools::Rev(phyloseq.file@tax_table@.Data, margin = 2)
  phyloseq.file@tax_table@.Data <- cbind(phyloseq.file@tax_table@.Data, lowest_tax_level=tax_fin)
}
