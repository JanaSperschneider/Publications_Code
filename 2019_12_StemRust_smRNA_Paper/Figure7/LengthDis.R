theme_Publication <- function(base_size=0.64, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.5), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1.5)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(face = "bold",size = rel(1.2)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            #panel.grid.major = element_line(colour="#f0f0f0"),
            #panel.grid.minor = element_blank(),
            #legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#8c510a","#80cdc1","#dfc27d","#01665e")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#8c510a","#80cdc1","#dfc27d","#01665e")), ...)
  
}
#--------------------------------------------------------------
library(ggplot2)
library(scales)
library(reshape2)
library(lattice)
library(ggrepel)

samples = c(18,19,20,21,22,23,24,25,26,27,28)
samples <- as.character(samples)

PGT_most_abundant_in_cluster = c(1.1 , 2.2 , 26.5 , 34.7 , 31.8 , 3.6 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0)
Wheat_most_abundant_in_cluster = c(0.0 , 0.0 , 9.0 , 49.9 , 14.8 , 3.9 , 22.4 , 0.0 , 0.0 , 0.0 , 0.0)

PGT_reads_in_clusters = c(3.3 , 8.8 , 23.7 , 26.2 , 25.4 , 9.0 , 1.7 , 0.7 , 0.5 , 0.4 , 0.4)
Wheat_reads_in_clusters = c(8.5 , 14.8 , 13.1 , 29.4 , 11.4 , 6.6 , 11.4 , 1.7 , 0.7 , 1.5 , 0.7)

df <- data.frame(PGT_most_abundant_in_cluster, PGT_reads_in_clusters, Wheat_most_abundant_in_cluster, Wheat_reads_in_clusters,
                 samples=c(samples))

head(df)
#---------------------------------
g1 <- ggplot(df, aes(x=samples, group=1)) + 
  geom_line(aes(y = PGT_most_abundant_in_cluster, colour = "Most abundant rust sRNA for each locus", size=0.4)) + 
  geom_point(aes(y = PGT_most_abundant_in_cluster, colour = "Most abundant rust sRNA for each locus", size=0.8)) + 
  geom_line(aes(y = PGT_reads_in_clusters, colour = "Read composition of rust sRNA loci", size=0.4)) + 
  geom_point(aes(y = PGT_reads_in_clusters, colour = "Read composition of rust sRNA loci", size=0.8)) + 
  scale_y_continuous(name="% of sequences", labels = comma) + 
  xlab("Sequence length") +   
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.major.y = element_blank()) +  
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill = guide_legend(override.aes = list(size = NULL))) +
  labs(title = "") 

g1 <- g1 + scale_fill_Publication() + scale_colour_Publication()
g1 <- g1 + guides(size = FALSE)
g1

g2 <- ggplot(df, aes(x=samples, group=1)) + 
  geom_line(aes(y = Wheat_most_abundant_in_cluster, colour = "Most abundant wheat sRNA for each locus", size=0.4)) +
  geom_point(aes(y = Wheat_most_abundant_in_cluster, colour = "Most abundant wheat sRNA for each locus", size=0.8)) +
  geom_line(aes(y = Wheat_reads_in_clusters, colour = "Read composition of wheat sRNA loci", size=0.4)) + 
  geom_point(aes(y = Wheat_reads_in_clusters, colour = "Read composition of wheat sRNA loci", size=0.8)) + 
  scale_y_continuous(name="% of sequences", labels = comma) + 
  xlab("Sequence length") +   
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.major.y = element_blank()) +  
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill = guide_legend(override.aes = list(size = NULL))) +
  labs(title = "") 

g2 <- g2 + scale_fill_Publication() + scale_colour_Publication()
g2 <- g2 + guides(size = FALSE)
g2
#---------------------------------
#---------------------------------
#---------------------------------
require(gridExtra)
#---------------------------------
png("ReadDistribution.png", height = 6, width = 14, units = 'in', res = 300)
grid.arrange(g1, g2, ncol=2)
dev.off()
#---------------------------------

