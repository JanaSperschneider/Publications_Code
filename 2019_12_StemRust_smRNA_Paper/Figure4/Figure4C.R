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
            legend.title = element_blank(),
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
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figures_Revision/Figure4")
#################################################################################################################################################################

df <- data.frame(methylation <- c("Germinated\nspores", "Late\ninfection"),
  young_TEs_centromeres <- c(21.24, 25.07),
  young_TEs_not_in_centromeres <- c(23.64, 23.62),
  old_TEs_centromeres <- c(10.19, 12.97),
  old_TEs_not_in_centromeres <- c(5.38, 4.95) 
)
names(df) <- c("Methylation", "Young TEs in centromeres", "Young TEs not in centromeres",
               "Old TEs in centromeres", "Old TEs not in centromeres")
head(df)
mdf <- melt(df, id="Methylation")
head(mdf)

g1 <- ggplot(mdf, aes(x=variable, y=value, fill=Methylation)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  ylab("% of TEs that are methylated") + xlab("") +
  coord_flip() + 
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() 
g1

png("TE_Age_Methylation.png", height = 8, width = 8, units = 'in', res = 300)
g1
dev.off()

#################################################################################################################################################################
#---------------
bins <- c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95)
age <- c("0", "0", "0", "0", "0", "0", "0", "0", "0", "1")
#---------------
# Centromeres
values <- c(0.0, 0.07, 2.88, 8.69, 5.88, 0.97, 5.74, 21.33, 26.25, 28.19)

df <- data.frame(bins, values, age)
g1 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +
  geom_bar(stat="identity") + ggtitle("Repeats in centromeres") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) + 
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") 
g1
#---------------
# Non-centromeres
values <- c(0.0, 0.05, 4.15, 11.85, 7.79, 1.68, 4.83, 28.57, 22.55, 18.52)

df <- data.frame(bins, values, age)
g2 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + ggtitle("Repeats outside centromeres") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) + 
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") 
g2
#---------------
require(gridExtra)
png("TE_Age.png", height = 6, width = 10, units = 'in', res = 300)
grid.arrange(g1, g2, ncol=2)
dev.off()
#---------------
values <- c(0.0, 0.01, 1.13, 2.55, 1.37, 0.3, 5.05, 17.64, 19.6, 52.36)
df <- data.frame(bins, values, age)
g3 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + 
  ggtitle("Repeats outside centromeres\nMethylation in late infection") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") + scale_y_continuous(limits = c(0, 55))
g3
#---------------
values <- c(0.0, 0.02, 1.23, 2.72, 1.42, 0.27, 4.83, 18.45, 20.86, 50.21)
df <- data.frame(bins, values, age)
g4 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + 
  ggtitle("Repeats outside centromeres\nMethylation in germinated spores") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none")  + scale_y_continuous(limits = c(0, 55))
g4
#---------------
values <- c(0.0, 0.08, 1.9, 4.71, 3.39, 0.08, 6.04, 17.87, 22.83, 43.09)
df <- data.frame(bins, values, age)
g5 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + 
  ggtitle("Repeats in centromeres\nMethylation in late infection") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none")  + scale_y_continuous(limits = c(0, 55))
g5
#---------------
values <- c(0.0, 0.0, 1.32, 3.96, 2.54, 0.2, 6.19, 16.04, 24.26, 45.48)
df <- data.frame(bins, values, age)
g6 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + 
  ggtitle("Repeats in centromeres\nMethylation in germinated spores") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none")  + scale_y_continuous(limits = c(0, 55))
g6
#---------------
require(gridExtra)
png("TE_Age_Methylation.png", height = 10, width = 10, units = 'in', res = 300)
grid.arrange(g3, g4, g5, g6, ncol=2)
dev.off()
#---------------
# smRNAs and repeats
#---------------
values <- c(0.0, 0.0, 1.97, 5.35, 2.54, 0.28, 8.73, 14.65, 25.07, 41.41)
df <- data.frame(bins, values, age)
g7 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + 
  ggtitle("Methylated repeats (GS) in centromeres that overlap\nwith smRNAs up-regulated during late infection") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") 
g7

values <- c(0.0, 0.0, 3.07, 5.04, 2.63, 0.22, 8.33, 14.69, 23.25, 42.76)
df <- data.frame(bins, values, age)
g8 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + 
  ggtitle("Methylated repeats (7dpi) in centromeres that overlap\nwith smRNAs up-regulated during late infection") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") 
g8

values <- c(0.0, 0.0, 2.4, 4.34, 2.54, 0.53, 5.68, 12.49, 12.09, 59.92)
df <- data.frame(bins, values, age)
g9 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + 
  ggtitle("Methylated repeats (GS) outside centromeres that overlap\nwith smRNAs up-regulated during late infection") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") 
g9

values <- c(0.0, 0.0, 1.81, 4.34, 1.93, 0.54, 5.24, 12.29, 11.93, 61.93)
df <- data.frame(bins, values, age)
g10 <- ggplot(data=df, aes(x = bins, y = values, fill=age)) +  geom_bar(stat="identity") + 
  ggtitle("Methylated repeats (7dpi) outside centromeres that overlap\nwith smRNAs up-regulated during late infection") +
  xlab("TE family % identity") + ylab("% of repeats") +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  theme_Publication() + scale_fill_Publication() + scale_colour_Publication() + 
  theme(legend.position="none") 
g10
#---------------
require(gridExtra)
png("TE_Age_smRNAs.png", height = 10, width = 15, units = 'in', res = 300)
grid.arrange(g7, g8, g9, g10, ncol=2)
dev.off()
