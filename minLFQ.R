library(dplyr)
library(splitstackshape)
library(here)
library("miscTools")
library(preprocessCore)
library(XML)
library(dplyr)
library(reshape2)
library(tidyr)
library(here)
library(xtable)


txtpath <- "C:/Users/Brandon/Documents/Data Analysis/Kim_reprocess/H_timecourse/"

direct <- txtpath

setwd(direct)

peptides <- read.delim(file = paste0(direct,"/txt/peptides.txt"))
modSpecPeptides <- read.delim(file = paste0(direct,"/txt/modificationSpecificPeptides.txt"))
proteinGroups <- read.delim(file = paste0(direct,"/txt/proteinGroups.txt"))
summary <- read.delim(file = paste0(direct,"/txt/summary.txt"))

PGs <- proteinGroups

uniprot.IDs <- data.frame(Protein.ID = PGs$Majority.protein.IDs, stringsAsFactors = F)
uniprot.IDs$Protein.ID <- as.character(uniprot.IDs$Protein.ID)
uniprot.IDs <- cSplit(indt = uniprot.IDs, splitCols = "Protein.ID", sep = ";")
uniprot.IDs <- data.frame(uniprot.IDs[,1])
colnames(uniprot.IDs) <- "Protein.ID"

intensity.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "Intensity.")])
intensity.PG$Protein.ID <- as.character(intensity.PG$Protein.ID)

peptides_mod <- cbind(peptides[grep(colnames(peptides), pattern = "^Sequence$")],
                      peptides[grep(colnames(peptides), pattern = "^Proteins$")],
                      peptides[grep(colnames(peptides), pattern = "Intensity.")])

modSpecPeptides_mod <- cbind(modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Sequence$")],
                             modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Modifications$")],
                             modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Proteins$")],
                             modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Charges$")],
                             modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "Intensity.")])

# Filter by protein
#
proteins_1st <- data.frame(First.protein=modSpecPeptides_mod$Proteins)
proteins_1st$First.protein <- (substr(proteins_1st$First.protein,1,9))
modSpecPeptides_mod <- cbind(proteins_1st,modSpecPeptides_mod)
modSpecPeptides_mod_2 <- modSpecPeptides_mod


start_time <- Sys.time()###
shared_all <- list()
test_df <- data.frame(matrix(ncol = length(colnames(peptide_membership[,-1])), nrow = length(colnames(peptide_membership[,-1]))))
colnames(test_df) <- sub("Intensity.","",colnames(peptide_membership[,-1]))
rownames(test_df) <- sub("Intensity.","",colnames(peptide_membership[,-1]))

# for (prot in unique(modSpecPeptides_mod_2$First.protein)){
for (prot in c("sp|O60341", "sp|Q86U42", "sp|Q96JP5")){
  # prot = "sp|O60341"
  modSpecPeptides_mod <- modSpecPeptides_mod_2[modSpecPeptides_mod_2$First.protein == prot,]

  # Peptide membership matrix (FIG. 2.C)
  # Each peptide different (modded peptide seen as separate peptide species)
  peptide_membership <- cbind(modSpecPeptides_mod[grep(colnames(modSpecPeptides_mod), pattern = "^Sequence$")],
                              modSpecPeptides_mod[grep(colnames(modSpecPeptides_mod), pattern = "Intensity.")])
  peptide_membership$Sequence <- as.character(peptide_membership$Sequence)

  peptide_membership_int <- peptide_membership[,-1]

  peptide_membership_int[peptide_membership_int > 0] <- 1
  peptide_membership_int[peptide_membership_int == 0] <- 0

  peptide_membership <- cbind(Sequence=peptide_membership$Sequence, peptide_membership_int)

  # Ratio count (validity cut-off = 2 peptides shared; (FIG. 2.D))

  # i = "Intensity.HT1BME1"
  # j = "Intensity.HT1BME2"

for (i in colnames(peptide_membership[,-1])){
  for (j in colnames(peptide_membership[,-1])){
    if (i != j){
      temp <- cbind(sequence=peptide_membership[,1],peptide_membership[,c(i,j)])
      temp2 <- rowSums(temp[,-1])
      temp <- cbind(temp, temp2)
      temp$temp2[temp$temp2 < 2] <- 0
      temp$temp2[temp$temp2 == 2] <- 1
      shared_peps <- sum(temp$temp2)
      i_2 <- sub("Intensity.","",i)
      j_2 <- sub("Intensity.","",j)
      # assign(paste0(i,"_",j),temp)
      # shared <- cbind(prot,shared_peps, exp=paste0(i_2,"_",j_2))
      test_df[i_2,j_2] <- shared_peps
      # assign(prot,test_df)
      test_df_half <- test_df
        # This section removes the top half of the matrix. Saves memory but runs slower (-50%)
          # test_df_half[upper.tri(test_df_half)]<-NA
          # test_df_half <- as.data.frame(test_df_half)
      shared_all[[prot]] <- test_df_half
      }
    }
  }
}
end_time <- Sys.time()###
print(end_time - start_time)



test <- shared_all


test_df_half[upper.tri(test_df_half)]<-NA
test_df_half <- as.data.frame(test_df_half)






#
#
#
# Seq <- data.frame(Seq = peptides$Sequence, stringsAsFactors = F)
# Seq$Seq <- as.character(Seq$Seq)
# # uniprot.IDs <- cSplit(indt = uniprot.IDs, splitCols = "Protein.ID", sep = ";")
# uniprot.IDs <- data.frame(uniprot.IDs[,1])
# colnames(uniprot.IDs) <- "Protein.ID"
#
# intensity.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "Intensity.")])
# intensity.PG$Protein.ID <- as.character(intensity.PG$Protein.ID)
#
# # Make quant file with common name despite LFQ/Intensity (makes downstream tasks easier)
# ProteinQuant.norm <- intensity.PG
#
#
#
#
#
#
#
#
#
#
#
#

