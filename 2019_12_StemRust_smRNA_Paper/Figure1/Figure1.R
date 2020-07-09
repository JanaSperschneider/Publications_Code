#################################################################################################################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("karyoploteR")
#BiocManager::install("GenomicFeatures")

library(karyoploteR)
library(GenomicFeatures)
set.seed(42)
########################################################
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figures_Revision/Figure1")
########################################################
########################################################
# Read in the data
########################################################
########################################################
centromeres <- toGRanges(
  data.frame(
    chromosome = c("chr1A", "chr1B", 
                   "chr2A", "chr2B",
                   "chr3A", "chr3B",
                   "chr4A", "chr4B",
                   "chr5A", "chr5B",
                   "chr6A", "chr6B", 
                   "chr7A", "chr7B",
                   "chr8A", "chr8B",
                   "chr9A", "chr9B",
                   "chr10A", "chr10B", 
                   "chr11A", "chr11B", 
                   "chr12A", "chr12B",
                   "chr13A", "chr13B",
                   "chr14A", "chr14B",
                   "chr15A", "chr15B",
                   "chr16A", "chr16B",
                   "chr17A", "chr17B",
                   "chr18A", "chr18B"),
    start = c(2360000, 2620000,
              1750000, 1650000,
              1950000, 2160000,
              3250000, 3450000,
              4550000, 6150000,
              1500000, 1300000,
              1900000, 2050000,
              1450000, 1450000,
              2050000, 2250000,
              1550000, 2150000,
              3600000, 3850000,
              3000000, 2900000,
              1950000, 2150000,
              2250000, 2350000,
              2700000, 2700000,
              2400000, 2350000,
              2450000, 2300000,
              1900000, 2100000), 
    end = c(2360000, 2620000,
            1750000, 1650000,
            1950000, 2160000,
            3250000, 3450000,
            4550000, 6150000,
            1500000, 1300000,
            1900000, 2050000,
            1450000, 1450000,
            2050000, 2250000,
            1550000, 2150000,
            3600000, 3850000,
            3000000, 2900000,
            1950000, 2150000,
            2250000, 2350000,
            2700000, 2700000,
            2400000, 2350000,
            2450000, 2300000,
            1900000, 2100000),    
    labels = c("Centromere", "Centromere", 
               "Centromere", "Centromere",
               "Centromere", "Centromere",
               "Centromere", "Centromere",
               "Centromere", "Centromere",
               "Centromere", "Centromere",
               "Centromere", "Centromere", 
               "Centromere", "Centromere",
               "Centromere", "Centromere",
               "Centromere", "Centromere", 
               "Centromere", "Centromere", 
               "Centromere", "Centromere",
               "Centromere", "Centromere",
               "Centromere", "Centromere", 
               "Centromere", "Centromere",
               "Centromere", "Centromere",
               "Centromere", "Centromere",
               "Centromere", "Centromere")))
########################################################
RNAseq_7dpi_coverage <- read.delim("../WI7_RNAseq_salignments.readcounts.bed", header=FALSE)
WI7_RPMs <- RNAseq_7dpi_coverage
WI7_RPMs[4] <- apply(WI7_RPMs[4],2,function(x){1000000*x/sum(x)})
head(WI7_RPMs)
WI7_RPMs <- WI7_RPMs[WI7_RPMs$V4 < 100,]
WI7_RPMs <- toGRanges(data.frame(chr=WI7_RPMs$V1, start=WI7_RPMs$V2, end=WI7_RPMs$V3, y=WI7_RPMs$V4))
WI7_RPMs

RNAseq_GS_coverage <- read.delim("../GS_RNAseq_alignments.readcounts.bed", header=FALSE)
GS_RPMs <- RNAseq_GS_coverage
GS_RPMs[4] <- apply(GS_RPMs[4],2,function(x){1000000*x/sum(x)})
head(GS_RPMs)
GS_RPMs <- GS_RPMs[GS_RPMs$V4 < 100,]
GS_RPMs <- toGRanges(data.frame(chr=GS_RPMs$V1, start=GS_RPMs$V2, end=GS_RPMs$V3, y=GS_RPMs$V4))
GS_RPMs
########################################################
########################################################
# Chromosome plots
########################################################
########################################################
custom.genome_haplotypeA <- toGRanges(data.frame(chr=c("chr1A", "chr2A", "chr3A", "chr4A", "chr5A", "chr6A", "chr7A", "chr8A", "chr9A",
                                                       "chr10A", "chr11A", "chr12A", "chr13A", "chr14A", "chr15A", "chr16A", "chr17A", "chr18A"), 
                                                 start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
                                                 end=c(6156719, 6062785, 6034716, 5967513, 5558200, 5554679, 5184328, 5113503, 4787822, 
                                                       4648557, 4640345, 3976801, 3570070, 3568213, 3494038, 3430719, 3318234, 2873395)))
custom.genome_haplotypeB <- toGRanges(data.frame(chr=c("chr1B", "chr2B", "chr3B", "chr4B", "chr5B", "chr6B", "chr7B", "chr8B", "chr9B",
                                                       "chr10B", "chr11B", "chr12B", "chr13B", "chr14B", "chr15B", "chr16B", "chr17B", "chr18B"), 
                                                 start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
                                                 end=c(6527790, 6111596, 4933499, 6361177, 7278493, 5249173, 5504590, 2822573, 5144719,
                                                       4975925, 4948285, 3939593, 3305433, 3562274, 3444477, 5892386, 2935969, 3060192
                                                 )))
################################
png("Chromosome_TranscriptionLevels_HaplotypeA.png", height = 20, width = 20, units = 'in', res = 300)

kp <- plotKaryotype(genome = custom.genome_haplotypeA, plot.type=2, cex=2.0)#, zoom=zoom.region)

kpPlotRibbon(kp, data=WI7_RPMs, ymin=0, ymax=100, y1=WI7_RPMs$V4, data.panel=1,
             border="darkred", col="darkred")

kpPlotRibbon(kp, data=GS_RPMs, ymin=0, ymax=100, y1=GS_RPMs$V4, data.panel=2,
             border="darkgreen", col="darkgreen")

kpPlotMarkers(kp, data=centromeres, labels=centromeres$labels, label.color = "darkblue", r1=1.5, text.orientation="horizontal", line.with=3)

dev.off()
################################
png("Chromosome_TranscriptionLevels_HaplotypeB.png", height = 20, width = 20, units = 'in', res = 300)

kp <- plotKaryotype(genome = custom.genome_haplotypeB, plot.type=2, cex=2.0)#, zoom=zoom.region)

kpPlotRibbon(kp, data=WI7_RPMs, ymin=0, ymax=100, y1=WI7_RPMs$V4, data.panel=1,
             border="darkred", col="darkred")

kpPlotRibbon(kp, data=GS_RPMs, ymin=0, ymax=100, y1=GS_RPMs$V4, data.panel=2,
             border="darkgreen", col="darkgreen")

kpPlotMarkers(kp, data=centromeres, labels=centromeres$labels, label.color = "darkblue", r1=1.5, text.orientation="horizontal", line.with=3)

dev.off()