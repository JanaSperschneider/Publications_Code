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
library(ggplot2)
library(reshape2)
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figures_Revision/Figure10")
#################################################################################################################################################################
# Record the % of methylated TEs that are sRNA_plus
df <- data.frame(methylation <- c("Late infection", "Germinated spores"),
                 sRNA_up_late <- c(33.22, 30.52),
                 sRNA_up_late_22 <- c(21.96, 19.82),
                 sRNA_up_late_22_A <- c(16.15, 14.21),
                 sRNA_noDE <- c(8.30, 7.58),
                 sRNAs_22 <- c(22.77,	20.56),
                 sRNAs_22A <- c(16.65, 14.65)
  )
names(df) <- c("Methylation", 
               "sRNAs up-regulated in late infection", 
               "sRNAs up-regulated in late infection (22 nt)", 
               "sRNAs up-regulated in late infection (22 nt, 5' A)", 
               "sRNAs with no differential expression", 
               "22 nt sRNAs", 
               "22 nt sRNAs with 5' A")
               
head(df)
mdf <- melt(df, id="Methylation")
head(mdf)

g1 <- ggplot(mdf, aes(x=variable, y=value, fill=Methylation)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  ggtitle("TEs+sRNA") +
  ylab("% of TEs+sRNA that are methylated") + xlab("") +
  coord_flip() + 
  scale_fill_Publication() + scale_colour_Publication() +
  theme_bw(base_size = 24, base_family = "Helvetica") +
  scale_y_continuous(limits = c(0, 55)) +
  theme(legend.title=element_blank(),legend.position="bottom") 
g1
#################################################################################################################################################################
df <- data.frame(methylation <- c("Late infection", "Germinated spores"),
                 sRNA_up_late <- c(40.58,		37.62),
                 sRNA_up_late_22 <- c(27.65, 25.45),
                 sRNA_up_late_22_A <- c(21.27, 19.30),                 
                 sRNA_noDE <- c(10.45,		9.77),
                 sRNAs_22 <- c(28.69,		26.46),
                 sRNAs_22A <- c(21.88,		19.86)
)

names(df) <- c("Methylation", 
               "sRNAs up-regulated in late infection", 
               "sRNAs up-regulated in late infection (22 nt)", 
               "sRNAs up-regulated in late infection (22 nt, 5' A)", 
               "sRNAs with no differential expression", 
               "22 nt sRNAs", 
               "22 nt sRNAs with 5' A")

head(df)
mdf <- melt(df, id="Methylation")
head(mdf)

g2 <- ggplot(mdf, aes(x=variable, y=value, fill=Methylation)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  ggtitle("Young TEs+sRNA") +
  ylab("% of TEs+sRNA that are methylated") + xlab("") +
  coord_flip() + 
  scale_fill_Publication() + scale_colour_Publication() +
  theme_bw(base_size = 24, base_family = "Helvetica") +
  scale_y_continuous(limits = c(0, 55)) +
  theme(legend.title=element_blank(),legend.position="bottom") 
g2
#################################################################################################################################################################
df <- data.frame(methylation <- c("Late infection", "Germinated spores"),
                 sRNA_up_late <- c(51.34, 49.10),
                 sRNA_up_late_22 <- c(34.22, 33.26),
                 sRNA_up_late_22_A <- c(27.85, 27.52),
                 sRNA_noDE <- c(10.57,		10.09),
                 sRNAs_22 <- c(34.73,		33.46),
                 sRNAs_22A <- c(28.18,		27.72)
)

names(df) <- c("Methylation", 
               "sRNAs up-regulated in late infection", 
               "sRNAs up-regulated in late infection (22 nt)", 
               "sRNAs up-regulated in late infection (22 nt, 5' A)", 
               "sRNAs with no differential expression", 
               "22 nt sRNAs", 
               "22 nt sRNAs with 5' A")

head(df)
mdf <- melt(df, id="Methylation")
head(mdf)

g3 <- ggplot(mdf, aes(x=variable, y=value, fill=Methylation)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  ggtitle("Young TEs+sRNA inside centromeres") +
  ylab("% of TEs+sRNA that are methylated") + xlab("") +
  coord_flip() + 
  scale_fill_Publication() + scale_colour_Publication() +
  theme_bw(base_size = 24, base_family = "Helvetica") +
  scale_y_continuous(limits = c(0, 55)) +
  theme(legend.title=element_blank(),legend.position="bottom") 
g3
#################################################################################################################################################################
df <- data.frame(methylation <- c("Late infection", "Germinated spores"),
                 sRNA_up_late <- c(39.71,		36.66),
                 sRNA_up_late_22 <- c(27.05,		24.63),
                 sRNA_up_late_22_A <- c(20.74,		18.61),
                 sRNA_noDE <- c(10.33,		9.54),
                 sRNAs_22 <- c(28.14, 25.71),
                 sRNAs_22A <- c(21.38, 19.20)
)

names(df) <- c("Methylation", 
               "sRNAs up-regulated in late infection", 
               "sRNAs up-regulated in late infection (22 nt)", 
               "sRNAs up-regulated in late infection (22 nt, 5' A)", 
               "sRNAs with no differential expression", 
               "22 nt sRNAs", 
               "22 nt sRNAs with 5' A")

head(df)
mdf <- melt(df, id="Methylation")
head(mdf)

g4 <- ggplot(mdf, aes(x=variable, y=value, fill=Methylation)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  ggtitle("Young TEs+sRNA outside centromeres") +
  ylab("% of TEs+sRNA that are methylated") + xlab("") +
  coord_flip() + 
  scale_fill_Publication() + scale_colour_Publication() +
  theme_bw(base_size = 24, base_family = "Helvetica") +
  scale_y_continuous(limits = c(0, 55)) +
  theme(legend.title=element_blank(),legend.position="bottom") 
g4
#################################################################################################################################################################
require(gridExtra)
png("Figure10.png", height = 20, width = 30, units = 'in', res = 300)
grid.arrange(g1, g2, g3, g4, ncol=2)
dev.off()
