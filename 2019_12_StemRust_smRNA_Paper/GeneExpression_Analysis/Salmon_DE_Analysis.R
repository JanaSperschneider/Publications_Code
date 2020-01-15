#--------------------------------------------
#--------------------------------------------
library(DESeq2)
library(ggplot2)
library(tximportData)
library(readr)
library(tximport)
#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
samples <- c('HSRNA1_quant', 'HSRNA2_quant', 'HSRNA3_quant',
             'PGTGS1_quant', 'PGTGS2_quant', 'PGTGS3_quant',
             'RW2I1.fastq.gz_quant', 'RW2I2.fastq.gz_quant', 'RW2I3.fastq.gz_quant',
             'RW3I1.fastq.gz_quant', 'RW3I2.fastq.gz_quant', 'RW3I3.fastq.gz_quant',
             'RW4I1.fastq.gz_quant', 'RW4I2.fastq.gz_quant', 'RW4I3.fastq.gz_quant',
             'RW5I1.fastq.gz_quant', 'RW5I2.fastq.gz_quant', 'RW5I3.fastq.gz_quant',
             'RW6I1.fastq.gz_quant', 'RW6I2.fastq.gz_quant', 'RW6I3.fastq.gz_quant',
             'RW7I1.fastq.gz_quant', 'RW7I2.fastq.gz_quant', 'RW7I3.fastq.gz_quant'
             )
samples

files <- file.path('.', "quants", samples, "quant.sf")
files
all(file.exists(files))

tx2gene <- read_delim("tx2gene.txt", delim='\t', col_names=FALSE)
head(tx2gene)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)

tail(txi$counts)
########################################################################################
write.table(txi$counts, file="Tx2_ReadCounts.txt", sep="\t")
write.table(txi$abundance, file="Tx2_Abundance.txt", sep="\t")
########################################################################################
sampleTable <- data.frame(condition = factor(rep(c("HS", "GS", "dpi2", "dpi3", 
                                                   "dpi4", "dpi5", "dpi6", "dpi7"), each = 3)))

rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
# dds is now ready for DESeq() see DESeq2 vignette
########################################################################################
dds <- dds[ rowSums(counts(dds)) > 1, ]
########################################################################################
dds <- DESeq(dds)
########################################################################################
res <- results(dds, contrast=c("condition","HS","GS"))
res
summary(res)
resSig <- subset(res, padj < 0.1)
write.csv(as.data.frame(resSig), 
          file="DEseq2_HaustoriavsGerminatedSpores.csv")
########################################################################################
res <- results(dds, contrast=c("condition","dpi2","GS"))
res
summary(res)
resSig <- subset(res, padj < 0.1)
write.csv(as.data.frame(resSig), 
          file="DEseq2_2dpivsGerminatedSpores.csv")
########################################################################################
res <- results(dds, contrast=c("condition","dpi3","GS"))
res
summary(res)
resSig <- subset(res, padj < 0.1)
write.csv(as.data.frame(resSig), 
          file="DEseq2_3dpivsGerminatedSpores.csv")
########################################################################################
res <- results(dds, contrast=c("condition","dpi4","GS"))
res
summary(res)
resSig <- subset(res, padj < 0.1)
write.csv(as.data.frame(resSig), 
          file="DEseq2_4dpivsGerminatedSpores.csv")
########################################################################################
res <- results(dds, contrast=c("condition","dpi5","GS"))
res
summary(res)
resSig <- subset(res, padj < 0.1)
write.csv(as.data.frame(resSig), 
          file="DEseq2_5dpivsGerminatedSpores.csv")
########################################################################################
res <- results(dds, contrast=c("condition","dpi6","GS"))
res
summary(res)
resSig <- subset(res, padj < 0.1)
write.csv(as.data.frame(resSig), 
          file="DEseq2_6dpivsGerminatedSpores.csv")
########################################################################################
res <- results(dds, contrast=c("condition","dpi7","GS"))
res
summary(res)
resSig <- subset(res, padj < 0.1)
write.csv(as.data.frame(resSig), 
          file="DEseq2_7dpivsGerminatedSpores.csv")
########################################################################################
########################################################################################
########################################################################################
rlogs <- rlog(dds, blind=FALSE)
write.table(assay(rlogs), file="DESeq2_rlogsFull.txt", sep="\t")

pcaData <- DESeq2::plotPCA(rlogs, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

head(pcaData)

p1 <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  scale_shape_manual(values=1:nlevels(pcaData$group)) +
  geom_point(aes(fill=pcaData$group),colour="black",pch=21, size=5) +
  scale_fill_brewer("Sample", palette="Spectral") 

p1

png("DeSeq2_PCAplot.png", height = 6, width = 8, units = 'in', res = 300)
p1
dev.off()
########################################################################################
