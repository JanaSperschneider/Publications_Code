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
            #legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.title = element_text(face="italic"),
            #plot.margin=unit(c(10,5,5,5),"mm"),
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

########################################################
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figure6/")
########################################################
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(GenomicFeatures)
library(karyoploteR)

seqlength = c(chr1A=6156719,chr2A=6062785,chr3A=6034716,chr4A=5967513,chr5A=5558200,chr6A=5554679,
              chr7A=5184328,chr8A=5113503, chr9A=4787822,chr10A=4648557,chr11A=4640345,chr12A=3976801,
              chr13A=3570070,chr14A=3568213,chr15A=3494038,chr16A=3430719,chr17A=3318234,chr18A=2873395)
windows <- tileGenome(seqlength, tilewidth = 200000, cut.last.tile.in.chrom = TRUE)

all.genes.expressed <- toGRanges("ExpressedGenes_KaryoplotR.txt")
head(all.genes.expressed)

repeats <- toGRanges("chrs_repeatmasker_no_simple_repeats.txt")
head(repeats)

smRNAs_late = read.delim(file = "smRNAs_up_late_infection_karyplotR.txt", sep="\t", header=FALSE)
smRNAs_late <- toGRanges(smRNAs_late)

smRNAs_spores = read.delim(file = "smRNAs_upSpores_karyplotR.txt", sep="\t", header=FALSE)
smRNAs_spores <- toGRanges(smRNAs_spores)

dens_genes <- countOverlaps(windows, all.genes.expressed)
dens_repeats <- countOverlaps(windows, repeats)
dens_smRNA_GS <- countOverlaps(windows, smRNAs_spores)
dens_smRNA_G7dpi <- countOverlaps(windows, smRNAs_late)

df <- data.frame(seqnames=seqnames(windows),
                 starts=start(windows)-1,
                 ends=end(windows),
                 dens_smRNA_G7dpi)
head(df)
options(scipen=10)
write.table(df, file = "dens_smRNA_G7dpi.txt", quote=FALSE, row.names = FALSE)

df_scatter <- data.frame(dens_genes, dens_repeats, 
                         dens_smRNA_GS, dens_smRNA_G7dpi)
smRNA_late_points <- subset(df_scatter, dens_smRNA_G7dpi > 5)
smRNA_GS <- subset(df_scatter, dens_smRNA_GS > 5)

df_scatter <- data.frame(dens_genes, dens_repeats)
head(df_scatter)

g1 <- ggplot(df_scatter, aes(dens_genes, dens_repeats)) + 
  stat_density2d(aes(fill=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="#91bfdb",high="#fc8d59") +
  geom_point(smRNA_late_points, mapping= aes(x=dens_genes,y=dens_repeats,size=dens_smRNA_G7dpi)) +
  xlab("Gene density (# of genes per 200 kb)") + 
  ylab("Repeat density (# of repeats per 200 kb)") +
  ggtitle("sRNAs up-regulated during late infection")  + 
  guides(size=FALSE)
g1

g1 <- g1 + theme_Publication() + theme(legend.position='bottom') + ylim(50, 270) + xlim(0, 60)
g1 <- g1 + labs(fill = "Density")
g1 <- g1 + theme(
  legend.key.size = unit(1.5, "cm"),
  legend.key.height = unit(0.5,"cm"))
g1

png("smRNAs_late_infection.png", height = 6, width = 10, units = 'in', res = 300)
g1
dev.off()

g2 <- ggplot(df_scatter, aes(dens_genes, dens_repeats)) + 
  stat_density2d(aes(fill=..level..),geom='polygon',colour='black') + 
  scale_fill_continuous(low="#91bfdb",high="#fc8d59") +
  geom_point(smRNA_GS, mapping= aes(x=dens_genes,y=dens_repeats, size=dens_smRNA_GS)) +
  xlab("Gene density (# of genes per 200 kb)") + 
  ylab("Repeat density (# of repeats per 200 kb)") +
  ggtitle("sRNAs up-regulated in germinated spores") + 
  guides(size=FALSE)

g2 <- g2 + theme_Publication() + theme(legend.position='bottom')+ ylim(50, 270) + xlim(0, 60)
g2 <- g2 + labs(fill = "Density")
g2 <- g2 + theme(
  legend.key.size = unit(1.5, "cm"),
  legend.key.height = unit(0.5,"cm"))
g2

png("smRNAs_germinated_spores.png", height = 6, width = 10, units = 'in', res = 300)
g2
dev.off()