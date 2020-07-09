#################################################################################################################################################################
########################################################
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figures_Revision/Figure3")
########################################################
library(sem)
library(pheatmap)
library(dendextend)
library(ape)
library(tidyverse)


makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
########################################################
x <- scan("distance_matrix.cran.txt")
dims <- floor(sqrt(length(x) * 2))
m <- matrix(NA, dims, dims)
m[upper.tri(m, diag = TRUE)] <- x
m <- t(m)
m

full_data <- makeSymm(m)
head(full_data)
dim(full_data)
names <- read.delim(file = "distance_matrix.cran.rows.txt", header = FALSE)
names
names <- c(
  "chr10A centromere",
  "chr10A non-centromere",
  "chr10B centromere",
  "chr10B  non-centromere",
  "chr11A centromere",
  "chr11A  non-centromere",
  "chr11B centromere",
  "chr11B  non-centromere",
  "chr12A centromere",
  "chr12A  non-centromere",
  "chr12B centromere",
  "chr12B non-centromere",
  "chr13A centromere",
  "chr13A non-centromere",
  "chr13B centromere",
  "chr13B non-centromere",
  "chr14A centromere",
  "chr14A non-centromere",
  "chr14B centromere",
  "chr14B non-centromere",
  "chr15A centromere",
  "chr15A non-centromere",
  "chr15B centromere",
  "chr15B non-centromere",
  "chr16A centromere",
  "chr16A non-centromere",
  "chr16B centromere",
  "chr16B non-centromere",
  "chr17A centromere",
  "chr17A non-centromere",
  "chr17B centromere",
  "chr17B non-centromere",
  "chr18A centromere",
  "chr18A non-centromere",
  "chr18B centromere",
  "chr18B non-centromere",
  "chr1A centromere",
  "chr1A non-centromere",
  "chr1B centromere",
  "chr1B non-centromere",
  "chr2A centromere",
  "chr2A non-centromere",
  "chr2B centromere",
  "chr2B non-centromere",
  "chr3A centromere",
  "chr3A non-centromere",
  "chr3B centromere",
  "chr3B non-centromere",
  "chr4A centromere",
  "chr4A non-centromere",
  "chr4B centromere",
  "chr4B non-centromere",
  "chr5A centromere",
  "chr5A non-centromere",
  "chr5B centromere",
  "chr5B non-centromere",
  "chr6A centromere",
  "chr6A non-centromere",
  "chr6B centromere",
  "chr6B non-centromere",
  "chr7A centromere",
  "chr7A non-centromere",
  "chr7B centromere",
  "chr7B non-centromere",
  "chr8A centromere",
  "chr8A non-centromere",
  "chr8B centromere",
  "chr8B non-centromere",
  "chr9A centromere",
  "chr9A non-centromere",
  "chr9B centromere",
  "chr9B non-centromere"
)
########################################################
rownames(full_data) <- names
colnames(full_data) <- names
head(full_data)

# Compute distances and hierarchical clustering
hc <- hclust(dd, method = "ward.D2")
dend <- as.dendrogram(hc)
########################################################
png("Dendogram_Flipped.png", height = 14, width = 10, units = 'in', res = 300)
par(mar = c(4,1,1,15))
dend %>%
  set("labels_col", value = c("#253494", "#b30000"), k=2) %>%
  set("branches_k_color", value = c("#253494", "#b30000"), k = 2) %>%
  set("labels_cex", 1) %>%
  plot(horiz=TRUE, xlim=c(1.0,0.0), xlab = "Mash distance")
dev.off()
########################################################
png("Dendogram.png", height = 8, width = 15, units = 'in', res = 300)
par(mar = c(15,4,1,1))
dend %>%
  set("labels_col", value = c("#253494", "#b30000"), k=2) %>%
  set("branches_k_color", value = c("#253494", "#b30000"), k = 2) %>%
  set("labels_cex", 1) %>%
  plot(horiz=FALSE, ylab = "Mash distance")
dev.off()

