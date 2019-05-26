x_matrix <- data.frame(matrix(ncol = length(x_vector), nrow = dim(temp)[1]))
colnames(x_matrix) <- x_vector

for (i in rownames(x_matrix)){
  for (j in colnames(x_matrix)){
    if (temp[i,"a_var"] == j){
      x_matrix[i,j] <- 1
    }
    if (temp[i,"b_var"] == j){
      x_matrix[i,j] <- 1
    }
  }
}
x_matrix[is.na(x_matrix)] <- 0

testQR <- qr.solve(x_matrix,temp$y)

names(testQR) <- names(xvars)

###

PGs <- proteinGroups
PGs <- PGs[PGs$Potential.contaminant != "+" & PGs$Reverse != "+",]


# Make Quantitation file (choose LFQ/Intensities from python GUI**) ####
uniprot.IDs <- data.frame(Protein.ID = PGs$Majority.protein.IDs, stringsAsFactors = F)
uniprot.IDs$Protein.ID <- as.character(uniprot.IDs$Protein.ID)
uniprot.IDs <- cSplit(indt = uniprot.IDs, splitCols = "Protein.ID", sep = ";")
uniprot.IDs <- data.frame(uniprot.IDs[,1])
colnames(uniprot.IDs) <- "Protein.ID"

intensity.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "Intensity.")])
# intensity.PG <- cbind(uniprot.IDs, PGs[grep(colnames(PGs), pattern = "LFQ.")])
intensity.PG$Protein.ID <- as.character(intensity.PG$Protein.ID)
# Make quant file with common name despite LFQ/Intensity (makes downstream tasks easier)
ProteinQuant.norm <- intensity.PG

# intensity column naming discrepancy
colnames(ProteinQuant.norm) <- gsub(pattern = "_", replacement = "", x = colnames(ProteinQuant.norm))
colnames(ProteinQuant.norm) <- sub(pattern = "T", replacement = "HT", x = colnames(ProteinQuant.norm))

# random NA assignments
ProteinQuant.norm <- ProteinQuant.norm[complete.cases(ProteinQuant.norm),]



normed <- ProteinQuant.norm[grep(x = ProteinQuant.norm$Protein.ID, "Q96JP5"),]

rowSums(normed[,-1])*(exp(1)^(0.3430506))






