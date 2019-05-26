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
## NOTE: This algorithm does NOT use FastLFQ (clarify)
## Optional Settings for future dev:
# require MS/MS for LFQ comparisons: For each pairwise peptide intensity comparison, at least one of the two peptides ahve been identified by MS/MS (set to ON by default here)



txtpath <- "C:/Users/Brandon/Documents/Data Analysis/Kim_reprocess/H_timecourse/"

direct <- txtpath

setwd(direct)

peptides <- read.delim(file = paste0(direct,"/txt/peptides.txt"))
modSpecPeptides <- read.delim(file = paste0(direct,"/txt/modificationSpecificPeptides.txt"))
proteinGroups <- read.delim(file = paste0(direct,"/txt/proteinGroups.txt"))
summary <- read.delim(file = paste0(direct,"/txt/summary.txt"))

peptides<- peptides[peptides$Potential.contaminant != "+" & peptides$Reverse != "+",]
modSpecPeptides <- modSpecPeptides[modSpecPeptides$Potential.contaminant != "+" & modSpecPeptides$Reverse != "+",]

PGs <- proteinGroups

uniprot.IDs <- data.frame(Protein.ID = PGs$Majority.protein.IDs, stringsAsFactors = F)
uniprot.IDs$Protein.ID <- as.character(uniprot.IDs$Protein.ID)
uniprot.IDs <- cSplit(indt = uniprot.IDs, splitCols = "Protein.ID", sep = ";")
uniprot.IDs <- data.frame(uniprot.IDs[,1])
colnames(uniprot.IDs) <- "Protein.ID"

intensity.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "Intensity.")])
intensity.PG$Protein.ID <- as.character(intensity.PG$Protein.ID)

peptides <- cbind(peptides[grep(colnames(peptides), pattern = "^Sequence$")],
                      peptides[grep(colnames(peptides), pattern = "^Protein.Groups$")],
                      peptides[grep(colnames(peptides), pattern = "^Proteins$")],
                      peptides[grep(colnames(peptides), pattern = "Intensity.")])

modSpecPeptides_mod <- cbind(modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Sequence$")],
                             modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Protein.Groups$")],
                             # modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Modifications$")],
                             modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Proteins$")],
                             # modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "^Charges$")],
                             modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "Identification.type")],
                             modSpecPeptides[grep(colnames(modSpecPeptides), pattern = "Intensity.")])


# Filter by protein
#
proteins_1st <- data.frame(First.protein=modSpecPeptides_mod$Proteins)
proteins_1st$First.protein <- (substr(proteins_1st$First.protein,1,9))
modSpecPeptides_mod <- cbind(proteins_1st,modSpecPeptides_mod)
modSpecPeptides_mod_2 <- modSpecPeptides_mod
modSpecPeptides_mod_2$Protein.Groups <- as.character(modSpecPeptides_mod_2$Protein.Groups)
proteinGroups$id <- as.character(proteinGroups$id)


start_time <- Sys.time()###
shared_all <- list()
# test_df <- data.frame(matrix(ncol = length(colnames(peptide_membership[,-1])), nrow = length(colnames(peptide_membership[,-1]))))
# colnames(test_df) <- sub("Intensity.","",colnames(peptide_membership[,-1]))
# rownames(test_df) <- sub("Intensity.","",colnames(peptide_membership[,-1]))

# for (prot in unique(modSpecPeptides_mod_2$First.protein)){
# for (prot in c("sp|O60341", "sp|Q86U42", "sp|Q96JP5")){
# for (prot in proteinGroups$id){
  for (prot in c("997", "2983")){
  # prot = "sp|O60341"
  # modSpecPeptides_mod <- modSpecPeptides_mod_2[modSpecPeptides_mod_2$Protein.Groups == prot,]
    modSpecPeptides_mod <- modSpecPeptides_mod_2[grep(paste0("\\b",prot,"\\b"),modSpecPeptides_mod_2$Protein.Groups, ),]

  # Peptide membership matrix (FIG. 2.C)
  # Each peptide different (modded peptide seen as separate peptide species)
  peptide_membership <- cbind(modSpecPeptides_mod[grep(colnames(modSpecPeptides_mod), pattern = "^Sequence$")],
                              modSpecPeptides_mod[grep(colnames(modSpecPeptides_mod), pattern = "Intensity.")])
  peptide_membership$Sequence <- as.character(peptide_membership$Sequence)

  match_type <- cbind(modSpecPeptides_mod[grep(colnames(modSpecPeptides_mod), pattern = "^Sequence$")],
                              modSpecPeptides_mod[grep(colnames(modSpecPeptides_mod), pattern = "Identification.type.")])
  match_type[,] <- data.frame(apply(match_type[,], 2, function(x) as.character(x)), stringsAsFactors = F)


  test_df <- data.frame(matrix(ncol = length(colnames(peptide_membership[,-1])), nrow = length(colnames(peptide_membership[,-1]))))
  colnames(test_df) <- sub("Intensity.","",colnames(peptide_membership[,-1]))
  rownames(test_df) <- sub("Intensity.","",colnames(peptide_membership[,-1]))

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
      i_2 <- sub("Intensity.","",i)
      j_2 <- sub("Intensity.","",j)
      temp <- cbind(sequence=peptide_membership[,1],peptide_membership[,c(i,j)])
      colnames(temp) <- c("sequence",i_2,j_2)
      temp_matchType <- cbind(match_type[,1],
                              match_type[,grep(i_2,colnames(match_type))],
                              match_type[,grep(j_2, colnames(match_type))])
      colnames(temp_matchType) <- c("sequence",i_2,j_2)

      for (t in temp_matchType){


      }





      apply(temp_matchType, 1, function(x) if (any(temp_matchType[1,] %in% "By matching")){temp_matchType[x,length(colnames(temp_matchType))+1] <- "yes"})


      temp2 <- rowSums(temp[,-1])
      temp <- cbind(temp, temp2)
      temp$temp2[temp$temp2 < 2] <- 0
      temp$temp2[temp$temp2 == 2] <- 1
      shared_peps <- sum(temp$temp2)

      test_df[i_2,j_2] <- shared_peps
      # assign(prot,test_df)
       shared_all[[prot]] <- test_df
      }
    }
  }
}
end_time <- Sys.time()###
print(end_time - start_time)

# removed unassigned peptides
sum(names(shared_all) == "")
shared_all_2 <- shared_all[names(shared_all) != ""]

# removed rev peptides
sum(grepl(x = names(shared_all_2), "REV"))
shared_all_2 <- shared_all_2[!grepl(x = names(shared_all_2), "REV")]

# removed con peptides
sum(grepl(x = names(shared_all_2), "CON"))
shared_all_2 <- shared_all_2[!grepl(x = names(shared_all_2), "CON")]


start_time <- Sys.time()###
shared_all_half <- list()
# for (i in 1:length(shared_all)){
for (i in names(shared_all_2)){
  tempTest <- shared_all_2[[i]]
  tempTest[upper.tri(tempTest)]<- NA
  # tempTest[is.na(tempTest)] <- ""
  tempTest <- as.data.frame(tempTest)
  shared_all_half[[i]] <- tempTest
}
end_time <- Sys.time()###
print(end_time - start_time)


start_time <- Sys.time()###
ratios <- list()
# for (i in 1:length(shared_all)){
for (i in names(shared_all_half)){
  tempTest <- shared_all_half[[i]]
  tempTest[tempTest < 2] <- 0
  tempTest[tempTest >= 2] <- 1
  ratios[[i]] <- tempTest
}
end_time <- Sys.time()###
print(end_time - start_time)


medianRatios <- list()
# for (k in names(ratios)){
for (k in c("997", "2983")){
  modSpecPeptides_mod <- modSpecPeptides_mod_2[grep(paste0("\\b",k,"\\b"),modSpecPeptides_mod_2$Protein.Groups),]
  colnames(modSpecPeptides_mod) <- sub("Intensity.","",colnames(modSpecPeptides_mod))
  ratioDF <- ratios[[k]]
  for (x in rownames(ratioDF)){
    for (y in colnames(ratioDF)){
        if (!is.na(ratioDF[x,y]) & ratioDF[x,y] == 1){
          temp <- modSpecPeptides_mod[,c(1,grep(pattern = x, x = colnames(modSpecPeptides_mod)), grep(pattern = y, x = colnames(modSpecPeptides_mod)))]
          temp[temp==0] <- NA
          temp <- temp[complete.cases(temp),]
          temp$ratio <- temp[,x]/temp[,y]
          med <- median(temp$ratio)
          ratioDF[x,y] <- med
        }
    }
  }
  # ratioDF[ratioDF == 0] <- NA
  medianRatios[[k]] <- ratioDF
}

equations <- list()
# for (prot in names(ratios)){
for (prot in c("997", "2983")){    for (k in colnames(medianRatios[[prot]])){
      if (mediansDF[j,k] != 0 & !is.na(mediansDF[j,k])){
        temp <- rbind(temp,cbind(y=mediansDF[j,k], a=j, b=k))
      }
    # log(mediansDF[j,k]) - log(j) + log(k)
    # log(mediansDF[j,k]) = j/k
    }

  temp <- NULL
  mediansDF <- medianRatios[[prot]]
  for (j in rownames(medianRatios[[prot]])){
  }
  temp <- data.frame(temp)
  temp$y <- as.character(temp$y)
  vars <- unique(c(as.character(temp$a), as.character(temp$b)))
  # names(vars) <- paste0("x",seq(1, length(vars), 1))
  # vars <- cbind(vars, paste0("x",seq(1, length(vars), 1)))
  xvars <- paste0("x",seq(1, length(vars), 1))
  names(xvars) <- vars
  temp$a_var <- xvars[match(temp$a, names(xvars))]
  temp$b_var <- xvars[match(temp$b, names(xvars))]
  equations[[prot]] <- temp
}

x_vector <- paste0("x",seq(1, length(vars), 1))

temp23 <- equations[["2983"]]












log(temp$ratio) - log(x) + log(y)
temp$ratio/(x*y)
