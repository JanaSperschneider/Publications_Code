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
#################################################################################################################################################################
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
