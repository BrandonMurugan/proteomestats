library(NMF)
# nmf.options(grid.patch=TRUE)#set to avoid blank first pdf page being created
# nmf.options(grid.patch=TRUE)
library(splitstackshape)
library(gplots)
library(tidyr)
require("GMD")
require(grid)
library(colorspace)
library(dtplyr)

library("colorspace")
pal <- choose_palette()


raw <- ProteinQuant.norm_2
annot <- colnames(ProteinQuant.norm[-1])
# annot <- cbind(annot, colnames(ProteinQuant.norm[-1]))
# annot <- t(annot)
annot <- gsub(x = annot, pattern = "Intensity.", replacement = "")
annot <- gsub(x = annot, pattern = "_[0-9]$", replacement = "")
# annot <- factor(annot)
annot <- data.frame(cbind(annot,annot))
annot$annot <- factor(annot)
annot$annot.1 <- NULL

raw_Tat0_Tat24 <- raw[colnames(raw) == grep(x = raw, pattern = "")]

raw <- raw[complete.cases(raw[,-1]),]

# raw$Protein.ID <- NULL

data <- raw

data <- log(data[,-1], 10)
#data[,][data[,] == -Inf] <- NA



data_matrix <- data.matrix(data)
rownames(data_matrix) <- raw[,1]
# data_matrix <- data_matrix[,-1]
#data_matrix <- data_matrix[,-44]
data_matrix_t <- t(data_matrix)

# data_matrix <- data_matrix[-1,]

# sitecol <- c("darksalmon", "firebrick4")
# sitecolList <- setNames(sitecol, c("JHB", "CPT"))

# #Color Brewer Palette selector:
# library("colorspace")
# pal <- choose_palette()
# #show hardcoded palette: pal(7) [if looking at 7 colours - can look at subset of seletced palette]

# agecol <- c("#023FA5", "#5465AB", "#7D87B9", "#A1A6C8", "#BEC1D4", "#D6D7DD", "#E2E2E2")
# cols <- c("#0000CD", "#66CD00", "#8B0A50", "#8B0000", "#FF4500", "#030303", "#FFFF00", "#228B22", "#8B4513", "#27408B")
cols <- c("#023FA5", "#6A76B2", "#A1A6C8", "#CBCDD9","#086D30", "#317B46", "#508C5E", "#76A480", "#000000", "#6C0000")
samplecolList <- setNames(cols, unique(annot$annot))
ann_colors = list(annot = samplecolList)

#par(mfrow=c(1,3))

pdf("C:/Users/Brandon/Documents/Data Analysis/Tat_timecourse/newColReruns(sepReps_txt)/proteome_heat.pdf", width = 40,height = 5,pointsize =2, paper = "special")
par(mar=c(0.5, 0.5, 0.5, 0.5))
#mat <- matrix(c(0,0,1,2,3,4),1)
#layout(mat, widths=c(0.1,0.1,0.1,0.1,0.1,0.18), heights=c(1,1,1,1,1,1))
#aheatmap(Matrix_baseline, Rowv = NA, Colv = NA, legend = T, labRow = NA, cellheight = 12, cellwidth = 0)
aheatmap(data_matrix_t, #fontsize = 10,
         # Rowv = T,
         Colv = T,
         # width = 100,
         # height = 50,
         distfun = "euclidean",
         hclustfun = "complete",
         #breaks = 0,
         annRow = annot,
         border_color = "black",
         #color = "heat",
         legend = T, #info = T,
         #main = "BV Data",
         annColors = ann_colors,
         # labRow = rownames(data_matrix_t),
         cellheight = 2,
         cellwidth = 2)
dev.off()


pdf("C:/Users/Brandon/Documents/Data Analysis/Tat_timecourse/newColReruns(sepReps_txt)/proteome_heat_long.pdf", height = 25, width = 10,pointsize =2, paper = "special")
par(mar=c(0.5, 0.5, 0.5, 0.5))
#mat <- matrix(c(0,0,1,2,3,4),1)
#layout(mat, widths=c(0.1,0.1,0.1,0.1,0.1,0.18), heights=c(1,1,1,1,1,1))
#aheatmap(Matrix_baseline, Rowv = NA, Colv = NA, legend = T, labRow = NA, cellheight = 12, cellwidth = 0)
aheatmap(data_matrix, #fontsize = 10,
         Rowv = T,cexCol = 4,
         Colv = T,
         # width = 100,
         # height = 50,
         distfun = "euclidean",
         hclustfun = "complete",
         #breaks = 0,
         annCol = annot,
         border_color = "black",
         #color = "heat",
         legend = T, #info = T,
         #main = "BV Data",
         annColors = ann_colors,
         # labRow = rownames(data_matrix_t),
         cellheight = 2,
         cellwidth = 10)
dev.off()



library(ggfortify)

autoplot(prcomp(raw[,-1]), data = annot, colour = 'annot',
         loadings = F, loadings.colour = 'blue',scale = T,
         loadings.label = F, loadings.label.size = 3)
