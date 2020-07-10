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
