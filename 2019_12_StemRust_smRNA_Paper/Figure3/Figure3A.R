#################################################################################################################################################################
library(karyoploteR)
library(genomation)
library(GenomicFeatures)
library(Biostrings)
library(seqinr)
library(SeqinR)
library(biomartr)
library(BSgenome)
set.seed(42)
########################################################
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figures_Revision/Figure3/")
########################################################
########################################################
# Read in the data
########################################################
########################################################
custom.genome <- toGRanges(data.frame(chr=c("chr1A", "chr2A", "chr3A", "chr4A", "chr5A", "chr6A", "chr7A", "chr8A", "chr9A",
                                            "chr10A", "chr11A", "chr12A", "chr13A", "chr14A", "chr15A", "chr16A", "chr17A", "chr18A",
                                            "chr1B", "chr2B", "chr3B", "chr4B", "chr5B", "chr6B", "chr7B", "chr8B", "chr9B",
                                            "chr10B", "chr11B", "chr12B", "chr13B", "chr14B", "chr15B", "chr16B", "chr17B", "chr18B"), 
                                      start=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
                                      end=c(6156719, 6062785, 6034716, 5967513, 5558200, 5554679, 5184328, 5113503, 4787822, 
                                            4648557, 4640345, 3976801, 3570070, 3568213, 3494038, 3430719, 3318234, 2873395,
                                            6527790, 6111596, 4933499, 6361177, 7278493, 5249173, 5504590, 2822573, 5144719,
                                            4975925, 4948285, 3939593, 3305433, 3562274, 3444477, 5892386, 2935969, 3060192
                                      )))
################################
all.genes.expressed <- toGRanges("../ExpressedGenes_KaryoplotR.txt")
head(all.genes.expressed)

repeats <- toGRanges("../chrs_repet_no_SSR.bed")
head(repeats)
################################
centromeres <- readBed("../centromeres.bed", track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)
################################
synteny <- read.delim("../map_chr_B_to_chr_A.primary_alignments.bed", header=FALSE)
head(synteny)

synteny.filtered <- subset(synteny, synteny$V3-synteny$V2 > 20000 & synteny$V6-synteny$V5 > 20000)

head(synteny.filtered)

start.regs <- toGRanges(data.frame(synteny.filtered$V1, synteny.filtered$V2, synteny.filtered$V3))
end.regs <- toGRanges(data.frame(synteny.filtered$V4, synteny.filtered$V5, synteny.filtered$V6))
########################################################
########################################################
# Chromosome plots
########################################################
########################################################
plot.params <- getDefaultPlotParams(plot.type=4)
plot.params$ideogramheight <- 40
plot.params$bottommargin <- 50
plot.params$ideogramlateralmargin <- 0.05
########################################################

for (x in c("1","2","3","4","5","6","7","8","9",
            "10","11","12","13", "14","15","16","17","18"))
{
  filename = paste("Chromosomes_Synteny_Chr",x,".png",sep="")
  
  png(filename, height = 10, width = 20, units = 'in', res = 300)
  
  chromosomeA = paste("chr", x, "A",sep="")
  chromosomeB = paste("chr", x, "B",sep="")
  
  kp <- plotKaryotype(genome = custom.genome, plot.type=4, chromosomes = c(chromosomeA, chromosomeB),
                      plot.params=plot.params, cex=1.8)#, zoom=zoom.region)
  kpAddBaseNumbers(kp, tick.dist=1000000, tick.len=5, add.units=TRUE, digits=2, minor.ticks=TRUE, 
                   minor.tick.dist=100000, minor.tick.len=2,  cex=1, tick.col=NULL, minor.tick.col=NULL, clipping=TRUE)
  kpPlotRegions(kp, data=centromeres, data.panel="ideogram", col="#FFEECC", layer.margin = 0.05, border="#FFCCAA", r0=0, r1=1.0)
  
  kpPlotLinks(kp, data=start.regs, data2=end.regs, col="#a1d99b")

  
dev.off()
}
########################################################
