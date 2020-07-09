#################################################################################################################################################################
#################################################################################################################################################################
theme_Publication <- function(base_size=12, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.5), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(face = "bold",size = rel(1)), 
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
  discrete_scale("fill","Publication",manual_pal(values = c("#dfc27d","#018571","#dfc27d","#018571",
                                                            "#dfc27d","#018571","#dfc27d","#018571",
                                                            "#dfc27d","#018571","#dfc27d","#018571",
                                                            "#dfc27d","#018571","#dfc27d","#018571",
                                                            "#dfc27d","#018571","#dfc27d","#018571",
                                                            "#dfc27d","#018571","#dfc27d","#018571",
                                                            "#dfc27d","#018571","#dfc27d","#018571")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#dfc27d","#018571","#dfc27d","#018571",
                                                              "#dfc27d","#018571","#dfc27d","#018571",
                                                              "#dfc27d","#018571","#dfc27d","#018571",
                                                              "#dfc27d","#018571","#dfc27d","#018571",
                                                              "#dfc27d","#018571","#dfc27d","#018571",
                                                              "#dfc27d","#018571","#dfc27d","#018571",
                                                              "#dfc27d","#018571","#dfc27d","#018571")), ...)
  
}

#################################################################################################################################################################
library("RColorBrewer")
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library("ggsignif")
library("reshape2")
library(textshape)
library(pheatmap)
########################################################
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figures_Revision/Figure4")
########################################################
set.seed(42)
########################################################
########################################################
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)
########################################################
df <- read.delim("Centromeric_RepeatCoverage.txt", header = TRUE, row.names=1)
head(df)

df_filtered <- subset(df, rowSums(df) > 100.0)
df_transposed <- as.data.frame(t(df_filtered))

breaksList = seq(0, 60, by = 1)

png("Centromeric_RepeatCoverage.png", height = 15, width = 10, units = 'in', res = 300)
pheatmap(df_transposed, cluster_cols=FALSE,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         main = "Centromeric region covered by repeat superfamily",
         fontsize_col = 16, fontsize_row = 16,
         display_numbers = TRUE, fontsize_number=16, cellheight=16) # Sets the breaks of the color scale as in breaksList
dev.off()
################################################
df <- read.delim("NonCentromeric_RepeatCoverage.txt", header = TRUE, row.names=1)
head(df)

df_filtered <- subset(df, rowSums(df) > 100.0)
df_transposed <- as.data.frame(t(df_filtered))

breaksList = seq(0, 60, by = 1)

png("NonCentromeric_RepeatCoverage.png", height = 15, width = 10, units = 'in', res = 300)
pheatmap(df_transposed, cluster_cols=FALSE,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         main = "Non-centromeric region covered by repeat superfamily",
         fontsize_col = 16, fontsize_row = 16,
         display_numbers = TRUE, fontsize_number=16, cellheight=16) # Sets the breaks of the color scale as in breaksList
dev.off()
################################################
################################################
df <- read.delim("centromeres_repeat_coverage.txt", header = FALSE, row.names=1)
head(df)

chromosomes = rownames(df)
centromere_coverage = df$V7

df <- read.delim("not_centromeres_repeat_coverage.txt", header = FALSE)
head(df)

x <- "chr1A"
chr1A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr1B"
chr1B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr2A"
chr2A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr2B"
chr2B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr3A"
chr3A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr3B"
chr3B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr4A"
chr4A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr4B"
chr4B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr5A"
chr5A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr5B"
chr5B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr6A"
chr6A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr6B"
chr6B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr7A"
chr7A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr7B"
chr7B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr8A"
chr8A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr8B"
chr8B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr9A"
chr9A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr9B"
chr9B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr10A"
chr10A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr10B"
chr10B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr11A"
chr11A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr11B"
chr11B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr12A"
chr12A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr12B"
chr12B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr13A"
chr13A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr13B"
chr13B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr14A"
chr14A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr14B"
chr14B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr15A"
chr15A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr15B"
chr15B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr16A"
chr16A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr16B"
chr16B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr17A"
chr17A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr17B"
chr17B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr18A"
chr18A_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)
x <- "chr18B"
chr18B_coverage = (subset(df, df$V1 == x)[1,]$V5 + subset(df, df$V1 == x)[2,]$V5)/(subset(df, df$V1 == x)[1,]$V6 + subset(df, df$V1 == x)[2,]$V6)

not_centromere_coverage = c(chr1A_coverage, chr1B_coverage,
                            chr2A_coverage, chr2B_coverage,
                            chr3A_coverage, chr3B_coverage,
                            chr4A_coverage, chr4B_coverage,
                            chr5A_coverage, chr5B_coverage,                            
                            chr6A_coverage, chr6B_coverage,
                            chr7A_coverage, chr7B_coverage,
                            chr8A_coverage, chr8B_coverage,
                            chr9A_coverage, chr9B_coverage,
                            chr10A_coverage, chr10B_coverage,
                            chr11A_coverage, chr11B_coverage,
                            chr12A_coverage, chr12B_coverage,
                            chr13A_coverage, chr13B_coverage,
                            chr14A_coverage, chr14B_coverage,
                            chr15A_coverage, chr15B_coverage,
                            chr16A_coverage, chr16B_coverage,
                            chr17A_coverage, chr17B_coverage,
                            chr18A_coverage, chr18B_coverage
                                                              )

data <- data.frame(chromosomes, centromere_coverage, not_centromere_coverage)
head(data)

mdf <- melt(data)
head(mdf)

g <- ggplot(mdf, aes(x = variable, y = value, fill=variable)) +
  geom_boxplot() + geom_jitter(width = 0.2) +
  xlab("") +   
  ylab("Fraction of bases covered\nby repetitive elements") +   
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(legend.title=element_blank(),legend.position="none") +
  #coord_flip()+
  geom_signif(comparisons = list(c("centromere_coverage","not_centromere_coverage")),
              map_signif_level=TRUE, test = "t.test")
g

g <- g + scale_fill_Publication() + scale_colour_Publication() + scale_x_discrete(breaks=c("centromere_coverage","not_centromere_coverage"),
                                                                                  labels=c("Centromeric\nregions", "Non-centromeric\nregions"))
g

png("Repeat_Coverage.png", height = 6, width = 10, units = 'in', res = 300)
g
dev.off()
