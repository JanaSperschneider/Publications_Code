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
library(pheatmap)
library(RColorBrewer)
library("ggsignif")
library("reshape2")
########################################################
########################################################
# Now look at methylated repeats and nearby genes
########################################################
DFplus <- read.delim("germinated_spores.CG.methylatedsites_repeats_genes_methylated_TPMs.txt", header = FALSE)
head(DFplus)
# The GS column 
logs_plus <- as.numeric(DFplus$V3)
distances <- as.numeric(DFplus$V10)
logs_plus <- log1p(logs_plus)
DFplus <- data.frame(logs_plus,distances)
head(DFplus)

DFminus <- read.delim("germinated_spores.CG.methylatedsites_repeats_genes_notmethylated_TPMs.txt", header = FALSE)
head(DFminus)
# The GS column 
logs_minus <- as.numeric(DFminus$V3)
distances <- as.numeric(DFminus$V10)
logs_minus <- log1p(logs_minus)
# The distance column
DFminus <- data.frame(logs_minus,distances)
head(DFminus)
################################################
################################################
plot.data <- rbind(
  data.frame(group = "Methylated TE overlapping gene", value = DFplus[DFplus$distances==0,]$logs),
  data.frame(group = "Non-methylated TE overlapping gene", value = DFminus[DFminus$distances==0,]$logs),
  
  data.frame(group = "Methylated TE near gene (0-500 bp)", value = DFplus[DFplus$distances < 500 & DFplus$distances > 0,]$logs),
  data.frame(group = "Non-methylated TE near gene (0-500 bp)", value = DFminus[DFminus$distances < 500 & DFminus$distances > 0,]$logs),
  
  data.frame(group = "Methylated TE near gene (500-1000 bp)", value = DFplus[DFplus$distances < 1000 & DFplus$distances > 500,]$logs),
  data.frame(group = "Non-methylated TE near gene (500-1000 bp)", value = DFminus[DFminus$distances < 1000 & DFminus$distances > 500,]$logs)
)

xlabs <- paste(levels(plot.data$group),"\n(N=",table(plot.data$group),")",sep="")

g1 <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() + 
  coord_flip()+
  geom_signif(comparisons = list(c("Methylated TE overlapping gene","Non-methylated TE overlapping gene")),
              map_signif_level=TRUE, test = "t.test")+    
  geom_signif(comparisons = list(c("Methylated TE near gene (0-500 bp)","Non-methylated TE near gene (0-500 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  geom_signif(comparisons = list(c("Methylated TE near gene (500-1000 bp)","Non-methylated TE near gene (500-1000 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  xlab("")+
  ylab("Gene expression levels in germinated spores") +  
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") + 
  scale_x_discrete(labels=xlabs) +
  ggtitle("Germinated spores methylation") + 
  ylim(0,13)

g1
################################################
DFplus <- read.delim("infected_leaves.CG.methylatedsites_repeats_genes_methylated_TPMs.txt", header = FALSE)
head(DFplus)
logs_plus <- as.numeric(DFplus$V9)
distances <- as.numeric(DFplus$V10)
logs_plus <- log1p(logs_plus)
DFplus <- data.frame(logs_plus,distances)
head(DFplus)

DFminus <- read.delim("infected_leaves.CG.methylatedsites_repeats_genes_notmethylated_TPMs.txt", header = FALSE)
head(DFminus)
logs_minus <- as.numeric(DFminus$V9)
distances <- as.numeric(DFminus$V10)
logs_minus <- log1p(logs_minus)
# The distance column
DFminus <- data.frame(logs_minus,distances)
head(DFminus)
################################################
################################################
plot.data <- rbind(
  data.frame(group = "Methylated TE overlapping gene", value = DFplus[DFplus$distances==0,]$logs),
  data.frame(group = "Non-methylated TE overlapping gene", value = DFminus[DFminus$distances==0,]$logs),
  
  data.frame(group = "Methylated TE near gene (0-500 bp)", value = DFplus[DFplus$distances < 500 & DFplus$distances > 0,]$logs),
  data.frame(group = "Non-methylated TE near gene (0-500 bp)", value = DFminus[DFminus$distances < 500 & DFminus$distances > 0,]$logs),
  
  data.frame(group = "Methylated TE near gene (500-1000 bp)", value = DFplus[DFplus$distances < 1000 & DFplus$distances > 500,]$logs),
  data.frame(group = "Non-methylated TE near gene (500-1000 bp)", value = DFminus[DFminus$distances < 1000 & DFminus$distances > 500,]$logs)
)

xlabs <- paste(levels(plot.data$group),"\n(N=",table(plot.data$group),")",sep="")

g2 <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() + 
  coord_flip()+
  geom_signif(comparisons = list(c("Methylated TE overlapping gene","Non-methylated TE overlapping gene")),
              map_signif_level=TRUE, test = "t.test")+    
  geom_signif(comparisons = list(c("Methylated TE near gene (0-500 bp)","Non-methylated TE near gene (0-500 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  geom_signif(comparisons = list(c("Methylated TE near gene (500-1000 bp)","Non-methylated TE near gene (500-1000 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  xlab("")+
  ylab("Gene expression levels in late infection") +  
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") + 
  scale_x_discrete(labels=xlabs) +
  ggtitle("Late infection methylation") + 
  ylim(0,13)

g2
################################################
require(gridExtra)
png("Figure11A.png", height = 6, width = 18, units = 'in', res = 300)
grid.arrange(g2, g1, ncol=2)
dev.off()
########################################################
########################################################
DFplus <- read.delim("infected_leaves.CG_repeats_sRNAplus_TPMs.txt", header = FALSE)
head(DFplus)
# The 7dpi column 
logs <- as.numeric(DFplus$V9)
logs <- log1p(logs)
DFplus <- data.frame(logs)
head(DFplus)

DFminus <- read.delim("infected_leaves.CG_repeats_sRNAminus_TPMs.txt", header = FALSE)
head(DFminus)
# The 7dpi column 
logs <- as.numeric(DFminus$V9)
logs <- log1p(logs)
DFminus <- data.frame(logs)
head(DFminus)
################################################
plot.data <- rbind(
                   data.frame(group = "Methylated TE+sRNA overlapping gene", value = DFplus$logs),
                   data.frame(group = "Methylated TE-sRNA overlapping gene", value = DFminus$logs))

xlabs <- paste(levels(plot.data$group),"\n(N=",table(plot.data$group),")",sep="")

g7 <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() + 
  coord_flip()+
  geom_signif(comparisons = list(c("Methylated TE+sRNA overlapping gene","Methylated TE-sRNA overlapping gene")),
              map_signif_level=TRUE, test = "t.test")+    
  xlab("")+
  ylab("Gene expression levels at 7 dpi") +  
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") + 
  scale_x_discrete(labels=xlabs) +
  ggtitle("Late infection methylation") + 
  ylim(0,13)

g7
################################################
DFplus <- read.delim("infected_leaves.CG_YoungRepeats_sRNAplus_TPMs.txt", header = FALSE)
head(DFplus)
# The 7dpi column 
logs <- as.numeric(DFplus$V9)
logs <- log1p(logs)
DFplus <- data.frame(logs)
head(DFplus)

DFminus <- read.delim("infected_leaves.CG_YoungRepeats_sRNAminus_TPMs.txt", header = FALSE)
head(DFminus)
# The 7dpi column 
logs <- as.numeric(DFminus$V9)
logs <- log1p(logs)
DFminus <- data.frame(logs)
head(DFminus)
################################################
plot.data <- rbind(
  data.frame(group = "Methylated young TE+sRNA overlapping gene", value = DFplus$logs),
  data.frame(group = "Methylated young TE-sRNA overlapping gene", value = DFminus$logs))

xlabs <- paste(levels(plot.data$group),"\n(N=",table(plot.data$group),")",sep="")

g7_young <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() + 
  coord_flip()+
  geom_signif(comparisons = list(c("Methylated young TE+sRNA overlapping gene","Methylated young TE-sRNA overlapping gene")),
              map_signif_level=TRUE, test = "t.test")+    
  xlab("")+
  ylab("Gene expression levels at 7 dpi") +  
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") + 
  scale_x_discrete(labels=xlabs) +
  ggtitle("Late infection methylation")

g7_young
################################################
################################################
DFplus <- read.delim("infected_leaves.CG_OldRepeats_sRNAplus_TPMs.txt", header = FALSE)
head(DFplus)
# The 7dpi column 
logs <- as.numeric(DFplus$V9)
logs <- log1p(logs)
DFplus <- data.frame(logs)
head(DFplus)

DFminus <- read.delim("infected_leaves.CG_OldRepeats_sRNAminus_TPMs.txt", header = FALSE)
head(DFminus)
# The 7dpi column 
logs <- as.numeric(DFminus$V9)
logs <- log1p(logs)
DFminus <- data.frame(logs)
head(DFminus)
################################################
plot.data <- rbind(
  data.frame(group = "Methylated old TE+sRNA overlapping gene", value = DFplus$logs),
  data.frame(group = "Methylated old TE-sRNA overlapping gene", value = DFminus$logs))

xlabs <- paste(levels(plot.data$group),"\n(N=",table(plot.data$group),")",sep="")

g7_old <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() + 
  coord_flip()+
  geom_signif(comparisons = list(c("Methylated old TE+sRNA overlapping gene","Methylated old TE-sRNA overlapping gene")),
              map_signif_level=TRUE, test = "t.test")+    
  xlab("")+
  ylab("Gene expression levels at 7 dpi") +  
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") + 
  scale_x_discrete(labels=xlabs) +
  ggtitle("Late infection methylation") 
g7_old
########################################################
png("Figure11B.png", height = 6, width = 18, units = 'in', res = 300)
grid.arrange(g7_young, g7_old, ncol=2)
dev.off()
