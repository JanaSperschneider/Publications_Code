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
set.seed(42)
########################################################
########################################################
########################################################
DFplus <- read.delim("smRNAs_up_late_infection_genes_TE_associated_plus_TPMs.txt", header = FALSE)
head(DFplus)
# The 7dpi column 
logs <- as.numeric(DFplus$V9)
logs <- log1p(logs)
# The distance column
distances <- DFplus$V10
DFplus <- data.frame(logs,distances)
head(DFplus)

DFminus <- read.delim("smRNAs_up_late_infection_genes_TE_associated_minus_TPMs.txt", header = FALSE)
head(DFminus)
# The 7dpi column 
logs <- as.numeric(DFminus$V9)
logs <- log1p(logs)
# The distance column
distances <- DFminus$V10
DFminus <- data.frame(logs,distances)
head(DFminus)
################################################
plot.data <- rbind(
                   data.frame(group = "sRNA+TE overlapping gene", value = DFplus[DFplus$distances==0,]$logs),
                   data.frame(group = "sRNA-TE overlapping gene", value = DFminus[DFminus$distances==0,]$logs),
                   
                   data.frame(group = "sRNA+TE near gene (0-250 bp)", value = DFplus[DFplus$distances < 250 & DFplus$distances >0,]$logs),
                   data.frame(group = "sRNA-TE near gene (0-250 bp)", value = DFminus[DFminus$distances < 250 & DFminus$distances >0,]$logs),
                   
                   data.frame(group = "sRNA+TE near gene (250-500 bp)", value = DFplus[DFplus$distances < 500 & DFplus$distances >= 250,]$logs),
                   data.frame(group = "sRNA-TE near gene (250-500 bp)", value = DFminus[DFminus$distances < 500 & DFminus$distances >= 250,]$logs),
                   
                   data.frame(group = "sRNA+TE near gene (500-1000 bp)", value = DFplus[DFplus$distances < 1000 & DFplus$distances >= 500,]$logs),
                   data.frame(group = "sRNA-TE near gene (500-1000 bp)", value = DFminus[DFminus$distances < 1000 & DFminus$distances >= 500,]$logs))

xlabs <- paste(levels(plot.data$group),"\n(N=",table(plot.data$group),")",sep="")

g7 <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() + 
  coord_flip()+
  geom_signif(comparisons = list(c("sRNA+TE overlapping gene","sRNA-TE overlapping gene")),
              map_signif_level=TRUE, test = "t.test")+    
  geom_signif(comparisons = list(c("sRNA+TE near gene (0-250 bp)","sRNA-TE near gene (0-250 bp)")),
              map_signif_level=TRUE, test = "t.test")+  
  geom_signif(comparisons = list(c("sRNA+TE near gene (250-500 bp)","sRNA-TE near gene (250-500 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  geom_signif(comparisons = list(c("sRNA+TE near gene (500-1000 bp)","sRNA-TE near gene (500-1000 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  xlab("")+
  ylab("Gene expression levels at 7 dpi") +  
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") + 
  scale_x_discrete(labels=xlabs)
g7

png("Results_7dpi_sRNAs.png", height = 6, width = 10, units = 'in', res = 300)
g7
dev.off()
########################################################
DFplus <- read.delim("smRNAs_upSpores_genes_TE_associated_plus_TPMs.txt", header = FALSE)
head(DFplus)
# The spores column 
logs <- as.numeric(DFplus$V3)
logs <- log1p(logs)
# The distance column
distances <- DFplus$V10
DFplus <- data.frame(logs,distances)
head(DFplus)

DFminus <- read.delim("smRNAs_upSpores_genes_TE_associated_minus_TPMs.txt", header = FALSE)
head(DFminus)
# The spores column 
logs <- as.numeric(DFminus$V3)
logs <- log1p(logs)
# The distance column
distances <- DFminus$V10
DFminus <- data.frame(logs,distances)
head(DFminus)

plot.data <- rbind(
  data.frame(group = "sRNA+TE overlapping gene", value = DFplus[DFplus$distances==0,]$logs),
  data.frame(group = "sRNA-TE overlapping gene", value = DFminus[DFminus$distances==0,]$logs),

  data.frame(group = "sRNA+TE near gene (0-250 bp)", value = DFplus[DFplus$distances < 250 & DFplus$distances >0,]$logs),
  data.frame(group = "sRNA-TE near gene (0-250 bp)", value = DFminus[DFminus$distances < 250 & DFminus$distances >0,]$logs),
  
  data.frame(group = "sRNA+TE near gene (250-500 bp)", value = DFplus[DFplus$distances < 500 & DFplus$distances >= 250,]$logs),
  data.frame(group = "sRNA-TE near gene (250-500 bp)", value = DFminus[DFminus$distances < 500 & DFminus$distances >= 250,]$logs),
  
  data.frame(group = "sRNA+TE near gene (500-1000 bp)", value = DFplus[DFplus$distances < 1000 & DFplus$distances >= 500,]$logs),
  data.frame(group = "sRNA-TE near gene (500-1000 bp)", value = DFminus[DFminus$distances < 1000 & DFminus$distances >= 500,]$logs))

xlabs <- paste(levels(plot.data$group),"\n(N=",table(plot.data$group),")",sep="")

gS <- ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() + 
  coord_flip()+
  geom_signif(comparisons = list(c("sRNA+TE overlapping gene","sRNA-TE overlapping gene")),
              map_signif_level=TRUE, test = "t.test")+    
  geom_signif(comparisons = list(c("sRNA+TE near gene (0-250 bp)","sRNA-TE near gene (0-250 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  geom_signif(comparisons = list(c("sRNA+TE near gene (250-500 bp)","sRNA-TE near gene (250-500 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  geom_signif(comparisons = list(c("sRNA+TE near gene (500-1000 bp)","sRNA-TE near gene (500-1000 bp)")),
              map_signif_level=TRUE, test = "t.test")+
  xlab("")+
  ylab("Gene expression levels in germinated spores") +  
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") + 
  scale_x_discrete(labels=xlabs)
gS

png("Results_GS_sRNAs.png", height = 6, width = 10, units = 'in', res = 300)
gS
dev.off()
########################################################
library(gridExtra)

png("TE_sRNAs.png", height = 6, width = 18, units = 'in', res = 300)
grid.arrange(g7,gS,ncol=2)
dev.off()

