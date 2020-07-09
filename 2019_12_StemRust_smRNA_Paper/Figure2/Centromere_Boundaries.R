#################################################################################################################################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("karyoploteR")
#BiocManager::install("GenomicFeatures")

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
########################################################
########################################################
# Chromosome plots
########################################################
plot.params <- getDefaultPlotParams(plot.type=2)
plot.params$ideogramheight <- 40
plot.params$data1height <- 60
plot.params$data1inmargin <- 1
plot.params$data2height <- 60
plot.params$data2inmargin <- 30

bases <- c("C", "G", "A", "T")
########################################################
########################################################
########################################################
library(BSgenome)
########################################################
# Zoom into the centromeric regions
centromeres <- readBed("centromeres.bed", track.line = FALSE, remove.unusual = FALSE, zero.based = TRUE)

all.genes.expressed <- toGRanges("ExpressedGenes_KaryoplotR.txt")
head(all.genes.expressed)

repeat_coverage <- read.delim("chr_A_B_unassigned_w1000_repeatcoverage.bed", header=FALSE)
repeat_coverage <- toGRanges(data.frame(chr=repeat_coverage$V1, start=repeat_coverage$V2, end=repeat_coverage$V3, y=repeat_coverage$V4))
head(repeat_coverage)

RNAseq_7dpi_coverage <- read.delim("WI7_RNAseq_alignments.readcounts.bed", header=FALSE)
WI7_RPMs <- RNAseq_7dpi_coverage
WI7_RPMs[4] <- apply(WI7_RPMs[4],2,function(x){1000000*x/sum(x)})
head(WI7_RPMs)
WI7_RPMs <- WI7_RPMs[WI7_RPMs$V4 < 500,]
WI7_RPMs <- toGRanges(data.frame(chr=WI7_RPMs$V1, start=WI7_RPMs$V2, end=WI7_RPMs$V3, y=WI7_RPMs$V4))
WI7_RPMs

RNAseq_GS_coverage <- read.delim("GS_RNAseq_alignments.readcounts.bed", header=FALSE)
GS_RPMs <- RNAseq_GS_coverage
GS_RPMs[4] <- apply(GS_RPMs[4],2,function(x){1000000*x/sum(x)})
head(GS_RPMs)
GS_RPMs <- GS_RPMs[GS_RPMs$V4 < 500,]
GS_RPMs <- toGRanges(data.frame(chr=GS_RPMs$V1, start=GS_RPMs$V2, end=GS_RPMs$V3, y=GS_RPMs$V4))
GS_RPMs
########################################################
png("chr2A_centromere.png", height = 10, width = 20, units = 'in', res = 300)

chromosome <- "chr2A"
x<- "2A"
# Zoom into region of interest
zoom.region <- toGRanges(data.frame(chromosome, 1e6, 2.5e6))

# A/T content graph
scaffold <- paste("Chromosome_FASTAs/chr",tolower(x),".fasta",sep="")
scaffold_id <- chromosome
x <- readDNAStringSet(scaffold)
dna_sequence <- x[[1]]
seqlength = c(id=length(dna_sequence))
window_size <- 1000 # configure the window size for GC content
custom.genome <- GRanges(seqnames=c(scaffold_id), ranges=IRanges(start=1, end=seqlength), strand="+")
tiles <- tile(x=custom.genome, width = window_size)
tiles <- unlist(tiles)
seqs <- getSeq(x, tiles)
counts <- alphabetFrequency(seqs, baseOnly=TRUE)
freqs <- counts/rowSums(counts)
mcols(tiles) <- DataFrame(freqs[,bases])
content <- lapply(seq_len(length(tiles)), 
                  function(i) {
                    return(as.numeric(data.frame(mcols(tiles))[i,c(1:4)]))
                  })
content <- do.call(rbind, content)
content <- data.frame(content)
names(content) <- bases
head(content)

kp <- plotKaryotype(genome = custom.genome, plot.type=2, chromosomes = c(chromosome),
                    plot.params=plot.params, cex=1.2, zoom=zoom.region, labels.plotter = NULL)
kpAddBaseNumbers(kp, tick.dist=100000, tick.len=5, add.units=TRUE, digits=2, minor.ticks=TRUE, 
                 minor.tick.dist=10000, minor.tick.len=2,  cex=1, tick.col=NULL, minor.tick.col=NULL, clipping=TRUE)
kpAddMainTitle(kp, main=chromosome, cex=2)
kpDataBackground(kp, data.panel = 1)

kpAxis(kp, ymin=0, ymax=500, data.panel = 1, r0=0.5, r1=1.0)
kpAddLabels(kp, labels="RPM             ", data.panel = 1, r0=0.5, r1=1.0)
kpLines(kp, data=GS_RPMs, ymin=0, ymax=500, y=GS_RPMs$V4, data.panel=1, col="darkgreen", r0=0.5, r1=1.0, lwd=2)
kpLines(kp, data=WI7_RPMs, ymin=0, ymax=500, y=WI7_RPMs$V4, data.panel=1, col="darkred", r0=0.5, r1=1.0, lwd=2)

kpPlotRegions(kp, data=centromeres, data.panel="ideogram", col="#FFEECC", layer.margin = 0.05, border="#FFCCAA", r0=0, r1=1.0)
kpPlotDensity(kp, all.genes.expressed, window.size = 1000, data.panel="ideogram", col="#67a9cf", border="#67a9cf", r0=0.5, r1=1)
kpAddLabels(kp, labels="Gene density", r0=0.5, r1=1, data.panel = "ideogram")

kpPlotRibbon(kp, data=repeat_coverage, data.panel=1, border="#ef8a62", col="#ef8a62", y0=0, y1=repeat_coverage$y, r0=0.0, r1=0.4)
kpAxis(kp, ymin=0.0, ymax=1.0, data.panel = 1, r0=0.0, r1=0.4)
kpAddLabels(kp, labels="Repeat        \ncoverage        ", data.panel = 1, r0=0, r1=0.4)

# Plot A/T content
kpLines(kp, data=tiles, y=content$C+content$G, data.panel = "ideogram" , lwd=2, r0=0.5, r1=0, col="black", ymin=0.2,ymax=0.8)
kpAbline(kp, h=0.435, col="red", r0=0.5, r1=0, ymin=0.2, ymax=0.8, data.panel = "ideogram")
kpAddLabels(kp, labels="% GC              ", r0=0.5, r1=0, data.panel = "ideogram")
kpAxis(kp, ymin=0.2,ymax=0.8, r0=0.5, r1=0, data.panel = "ideogram")

dev.off()

########################################################
png("chr1A_centromere.png", height = 10, width = 20, units = 'in', res = 300)

chromosome <- "chr1A"
x<- "1A"
# Zoom into region of interest
zoom.region <- toGRanges(data.frame(chromosome, 1.5e6, 3e6))

# A/T content graph
scaffold <- paste("Chromosome_FASTAs/chr",tolower(x),".fasta",sep="")
scaffold_id <- chromosome
x <- readDNAStringSet(scaffold)
dna_sequence <- x[[1]]
seqlength = c(id=length(dna_sequence))
window_size <- 1000 # configure the window size for GC content
custom.genome <- GRanges(seqnames=c(scaffold_id), ranges=IRanges(start=1, end=seqlength), strand="+")
tiles <- tile(x=custom.genome, width = window_size)
tiles <- unlist(tiles)
seqs <- getSeq(x, tiles)
counts <- alphabetFrequency(seqs, baseOnly=TRUE)
freqs <- counts/rowSums(counts)
mcols(tiles) <- DataFrame(freqs[,bases])
content <- lapply(seq_len(length(tiles)), 
                  function(i) {
                    return(as.numeric(data.frame(mcols(tiles))[i,c(1:4)]))
                  })
content <- do.call(rbind, content)
content <- data.frame(content)
names(content) <- bases
head(content)

kp <- plotKaryotype(genome = custom.genome, plot.type=2, chromosomes = c(chromosome),
                    plot.params=plot.params, cex=1.2, zoom=zoom.region, labels.plotter = NULL)
kpAddBaseNumbers(kp, tick.dist=100000, tick.len=5, add.units=TRUE, digits=2, minor.ticks=TRUE, 
                 minor.tick.dist=10000, minor.tick.len=2,  cex=1, tick.col=NULL, minor.tick.col=NULL, clipping=TRUE)
kpAddMainTitle(kp, main=chromosome, cex=2)
kpDataBackground(kp, data.panel = 1)

kpAxis(kp, ymin=0, ymax=500, data.panel = 1, r0=0.5, r1=1.0)
kpAddLabels(kp, labels="RPM             ", data.panel = 1, r0=0.5, r1=1.0)
kpLines(kp, data=GS_RPMs, ymin=0, ymax=500, y=GS_RPMs$V4, data.panel=1, col="darkgreen", r0=0.5, r1=1.0, lwd=2)
kpLines(kp, data=WI7_RPMs, ymin=0, ymax=500, y=WI7_RPMs$V4, data.panel=1, col="darkred", r0=0.5, r1=1.0, lwd=2)

kpPlotRegions(kp, data=centromeres, data.panel="ideogram", col="#FFEECC", layer.margin = 0.05, border="#FFCCAA", r0=0, r1=1.0)
kpPlotDensity(kp, all.genes.expressed, window.size = 1000, data.panel="ideogram", col="#67a9cf", border="#67a9cf", r0=0.5, r1=1)
kpAddLabels(kp, labels="Gene density", r0=0.5, r1=1, data.panel = "ideogram")

kpPlotRibbon(kp, data=repeat_coverage, data.panel=1, border="#ef8a62", col="#ef8a62", y0=0, y1=repeat_coverage$y, r0=0.0, r1=0.4)
kpAxis(kp, ymin=0.0, ymax=1.0, data.panel = 1, r0=0.0, r1=0.4)
kpAddLabels(kp, labels="Repeat        \ncoverage        ", data.panel = 1, r0=0, r1=0.4)

# Plot A/T content
kpLines(kp, data=tiles, y=content$C+content$G, data.panel = "ideogram" , lwd=2, r0=0.5, r1=0, col="black", ymin=0.2,ymax=0.8)
kpAbline(kp, h=0.435, col="red", r0=0.5, r1=0, ymin=0.2, ymax=0.8, data.panel = "ideogram")
kpAddLabels(kp, labels="% GC              ", r0=0.5, r1=0, data.panel = "ideogram")
kpAxis(kp, ymin=0.2,ymax=0.8, r0=0.5, r1=0, data.panel = "ideogram")

dev.off()
########################################################
