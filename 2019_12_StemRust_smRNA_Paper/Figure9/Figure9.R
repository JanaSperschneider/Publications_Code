#################################################################################################################################################################
library(karyoploteR)
library(genomation)
library(GenomicFeatures)
library(Biostrings)
library(seqinr)
library(SeqinR)
library(biomartr)
library(BSgenome)
########################################################
########################################################
########################################################
# Read in the data
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

########################################################
centromeres <- readBed("centromeres.bed", track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)

WI7coverage <- read.delim("WI7_alignments.readcounts.bed", header=FALSE)
WI7_RPMs <- WI7coverage
WI7_RPMs_scaled <- WI7_RPMs[4]*0.171953
WI7_RPMs[4] <- WI7_RPMs_scaled
head(WI7_RPMs)
WI7_RPMs <- WI7_RPMs[WI7_RPMs$V4 < 1000,]
WI7_RPMs <- toGRanges(data.frame(chr=WI7_RPMs$V1, start=WI7_RPMs$V2, end=WI7_RPMs$V3, y=WI7_RPMs$V4))
WI7_RPMs

GScoverage <- read.delim("GS_alignments.readcounts.bed", header=FALSE)
GS_RPMs <- GScoverage
GS_RPMs_scaled <- GS_RPMs[4]*0.036315
GS_RPMs[4] <- GS_RPMs_scaled
head(GS_RPMs)
GS_RPMs <- GS_RPMs[GS_RPMs$V4 < 1000,]
GS_RPMs <- toGRanges(data.frame(chr=GS_RPMs$V1, start=GS_RPMs$V2, end=GS_RPMs$V3, y=GS_RPMs$V4))
GS_RPMs
########################################################
########################################################
# Chromosome plots
########################################################
########################################################
png("Chromosome_smRNAs_HaplotypeA.png", height = 30, width = 20, units = 'in', res = 300)

kp <- plotKaryotype(genome = custom.genome_haplotypeA, plot.type=2, cex=2.0)#, zoom=zoom.region)
kpPlotRegions(kp, data=centromeres, col="#FFEECC", layer.margin = 0.05, border="#FFCCAA", r0=0, r1=1.0)

kpPlotRibbon(kp, data=WI7_RPMs, ymin=0, ymax=1000, y1=WI7_RPMs$y, data.panel=1, border="darkred", col="darkred", r0=0, r1=1)
kpPlotRibbon(kp, data=GS_RPMs, ymin=0, ymax=1000, y1=GS_RPMs$y, data.panel=2, border="darkgreen", col="darkgreen", r0=0, r1=1)

kpAxis(kp, ymin = 0, ymax=1000, data.panel=1)
kpAxis(kp, ymin = 0, ymax=1000, data.panel=2)

dev.off()
########################################################
png("Chromosome_smRNAs_HaplotypeB.png", height = 30, width = 20, units = 'in', res = 300)

kp <- plotKaryotype(genome = custom.genome_haplotypeB, plot.type=2, cex=1.9)#, zoom=zoom.region)
kpPlotRegions(kp, data=centromeres, col="#FFEECC", layer.margin = 0.05, border="#FFCCAA", r0=0, r1=1.0)

kpPlotRibbon(kp, data=WI7_RPMs, ymin=0, ymax=1000, y1=WI7_RPMs$y, data.panel=1, border="darkred", col="darkred", r0=0, r1=1)
kpPlotRibbon(kp, data=GS_RPMs, ymin=0, ymax=1000, y1=GS_RPMs$y, data.panel=2, border="darkgreen", col="darkgreen", r0=0, r1=1)

kpAxis(kp, ymin = 0, ymax=1000, data.panel=1)
kpAxis(kp, ymin = 0, ymax=1000, data.panel=2)

dev.off()
########################################################
########################################################
