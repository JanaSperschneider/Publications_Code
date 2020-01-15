#################################################################################################################################################################
#################################################################################################################################################################
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.2)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
#################################################################################################################################################################
library(karyoploteR)
library(GenomicFeatures)
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

all.genes.expressed <- toGRanges("ExpressedGenes_KaryoplotR.txt")
head(all.genes.expressed)

repeats <- toGRanges("chrs_repeatmasker_no_simple_repeats.txt")
head(repeats)

gypsys <- toGRanges("chrs_repeatmasker_gypsys.txt")
head(gypsys)

DNAs <- toGRanges("chrs_repeatmasker_DNAtransposons.txt")
head(DNAs)

smRNAs_late = read.delim(file = "smRNAs_up_late_infection_karyplotR.txt", sep="\t", header=FALSE)
smRNAs_late <- toGRanges(smRNAs_late)

smRNAs_early = read.delim(file = "smRNAs_up_early_infection_karyplotR.txt", sep="\t", header=FALSE)
smRNAs_early <- toGRanges(smRNAs_early)

smRNAs_spores = read.delim(file = "smRNAs_upSpores_karyplotR.txt", sep="\t", header=FALSE)
smRNAs_spores <- toGRanges(smRNAs_spores)

smRNAs_noDE = read.delim(file = "smRNAs_noDE_karyplotR.txt", sep="\t", header=FALSE)
smRNAs_noDE <- toGRanges(smRNAs_noDE)

genes_of_interest <- toGRanges(
  data.frame(
    chromosome = c("chr14B", "chr14A", 
                   "chr9A", "chr9B", 
                   "chr4A", "chr4A", 
                   "chr4B", "chr4B"),
    start = c(2553060, 2504240, 
              1423057, 1502987, 
              4199566, 4202372, 
              4409661, 4413251),
    end = c(2553458, 2507288, 
            1425317, 1505069, 
            4202185, 4204669, 
            4413155, 4415546),
    gene = c("AvrSr50", "AvrSr35", 
             "STE3.3", "STE3.2", 
             "bW3", "bE3", 
             "bW1", "bE1")))
head(genes_of_interest)

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
########################################################
# Chromosome plots
########################################################
bam1 <- "merged_alignments.bam"
bamGS <- "GS_alignments.bam"
bamWI7 <- "WI7_alignments.bam"

rnaseqWi7 <- "RW7I.bam"
rnaseqGS <- "PGTGS.bam"
rnaseqHS <- "HSRNA.bam"
rnaseqWi3 <- "RW3I.bam"
########################################################
plot.params <- getDefaultPlotParams(plot.type=2)
plot.params$ideogramheight <- 40
plot.params$data1height <- 40
plot.params$data1inmargin <- 1
plot.params$data2height <- 40
plot.params$data2inmargin <- 1
########################################################
for (x in c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B","8A","8B","9A","9B",
            "10A","10B","11A","11B","12A","12B","13A","13B","14A","14B","15A","15B","16A","16B","17A","17B","18A","18B"))
{
  filename = paste("Chromosome",x,".png",sep="")
  
  png(filename, height = 10, width = 20, units = 'in', res = 300)
  
  chromosome = paste("chr",x,sep="")
  chromosome
  
  kp <- plotKaryotype(genome = custom.genome, plot.type=2, chromosomes = c(chromosome),
                      plot.params=plot.params, cex=1.2, labels.plotter = NULL)#, zoom=zoom.region)
  
  kpAddBaseNumbers(kp, tick.dist=500000, tick.len=5, add.units=TRUE, digits=2, minor.ticks=TRUE, 
                   minor.tick.dist=100000, minor.tick.len=2,  cex=1, tick.col=NULL, minor.tick.col=NULL, clipping=TRUE)
  
  kpDataBackground(kp, data.panel = 1)
  kpDataBackground(kp, data.panel = 2)
  
  kp <- kpPlotBAMDensity(kp, data=rnaseqWi7, window.size=20000, normalize=TRUE, data.panel=2, r0=0.45, r1=0.0, col="darkblue")
  kpText(kp, chr=chromosome, x=0, y=0.225, data.panel = 2, col="darkblue", label="GS RNAseq\nread density", cex=1.2, pos=2)
  
  kp <- kpPlotBAMDensity(kp, data=rnaseqGS, window.size=20000, normalize=TRUE, data.panel=2, r0=0.95, r1=0.5, col="darkgreen")
  kpText(kp, chr=chromosome, x=0, y=0.725, data.panel = 2, col="darkgreen", label="7dpi RNAseq\nread density", cex=1.2, pos=2)
  
  kp <- kpPlotBAMDensity(kp, data=bamGS, window.size=20000, normalize=TRUE, data.panel=1, r0=0, r1=0.45, col="darkgreen")
  kpText(kp, chr=chromosome, x=0, y=0.225, data.panel = 1, col="darkgreen", label="GS sRNA\nread density", cex=1.2, pos=2)
  
  kp <- kpPlotBAMDensity(kp, data=bamWI7, window.size=20000, normalize=TRUE, data.panel=1, r0=0.5, r1=0.95, col="darkred")
  kpText(kp, chr=chromosome, x=0, y=0.725, data.panel = 1, col="darkred", label="7dpi sRNA\nread density", cex=1.2, pos=2)
  
  kpText(kp, chr=seqlevels(kp$genome), y=0.0, x=0, col="black", label="Genes", data.panel="ideogram", r0=0.75, cex=1.2, pos=2)
  kpPlotDensity(kp, all.genes.expressed, window.size = 20000, data.panel="ideogram", col="#67a9cf", border="#67a9cf", r0=0.5, r1=1)
  kpText(kp, chr=seqlevels(kp$genome), y=0.0, x=0, col="black", label="Repeats", data.panel="ideogram", r0=0.25, cex=1.2, pos=2)
  kpPlotDensity(kp, data = repeats, window.size = 20000, data.panel="ideogram", col="#ef8a62", border="#ef8a62", r0=0.5, r1=0)
  
  kpPlotMarkers(kp, data=genes_of_interest, labels=genes_of_interest$gene, label.color = "#333333", r1=1.5)
  
  kpPlotMarkers(kp, data=centromeres, labels=centromeres$labels, label.color = "#333333", r1=1.5)
  
  dev.off()
}
########################################################
