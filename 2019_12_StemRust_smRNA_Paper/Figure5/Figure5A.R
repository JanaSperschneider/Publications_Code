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
#--------------------------------------------------------------
#--------------------------------------------------------------
dinucleotides = c('CA', 'CC', 'CG', 'CT')
#---------------------------------
infected_leaves_freqs = c(0.8665, 
                          1.8750,
                          19.1904,
                          1.0066)

germinated_spores_freqs = c(0.7112,
                            1.5754,
                            17.0293,
                            0.8461)

df <- data.frame(infected_leaves_freqs, germinated_spores_freqs,
                 dinucleotides=c(dinucleotides))
names(df) <- c("Infected leaves", "Germinated spores", "Dinucleotide")

head(df)
mdf <- melt(df)
head(mdf)
#---------------------------------
g <- ggplot(mdf, aes(x = Dinucleotide, y = value, fill=variable)) +
  geom_bar(stat='identity', position = "dodge") +
  xlab("Dinucleotide") +   
  ylab("% of dinucleotides\nthat are methylated") +   
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.major.y = element_blank()) +  
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill = guide_legend(override.aes = list(size = NULL))) +
  labs(title = "5mC methylation in centromeres") 

g1_centromeres <- g + scale_fill_Publication() + scale_colour_Publication()
g1_centromeres <- g1_centromeres + guides(size = FALSE) + scale_y_continuous(limits = c(0, 20))
g1_centromeres
#---------------------------------
#---------------------------------
infected_leaves_freqs = c(0.7400, 
                          1.3869,
                          7.0923,
                          0.7723)

germinated_spores_freqs = c(0.6821,
                            1.2829,
                            7.1090,
                            0.7178)

df <- data.frame(infected_leaves_freqs, germinated_spores_freqs,
                 dinucleotides=c(dinucleotides))
names(df) <- c("Infected leaves", "Germinated spores", "Dinucleotide")

head(df)
mdf <- melt(df)
head(mdf)
#---------------------------------
g <- ggplot(mdf, aes(x = Dinucleotide, y = value, fill=variable)) +
  geom_bar(stat='identity', position = "dodge") +
  xlab("Dinucleotide") +   
  ylab("% of dinucleotides\nthat are methylated") +   
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.major.y = element_blank()) +  
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill = guide_legend(override.aes = list(size = NULL))) +
  labs(title = "5mC methylation outside the centromeres") 

g1_not_centromeres <- g + scale_fill_Publication() + scale_colour_Publication()
g1_not_centromeres <- g1_not_centromeres + guides(size = FALSE) + scale_y_continuous(limits = c(0, 20))
g1_not_centromeres
#---------------------------------
require(gridExtra)
png("5mC_DinucleotideFrequencies.png", height = 5, width = 15, units = 'in', res = 300)
grid.arrange(g1_centromeres, g1_not_centromeres, ncol=2)
dev.off()
#---------------------------------
#---------------------------------
#---------------------------------
dinucleotides = c('AA', 'AC', 'AG', 'AT')

infected_leaves_freqs = c(1.0860,
                          1.3178,
                          1.7659,
                          1.1700)

germinated_spores_freqs = c(0.8366,
                            1.0764,
                            1.3996,
                            0.9223)

df <- data.frame(infected_leaves_freqs, germinated_spores_freqs,
                 dinucleotides=c(dinucleotides))
names(df) <- c("Infected leaves", "Germinated spores", "Dinucleotide")

head(df)
mdf <- melt(df)
head(mdf)
#---------------------------------
g <- ggplot(mdf, aes(x = Dinucleotide, y = value, fill=variable)) +
  geom_bar(stat='identity', position = "dodge") +
  xlab("Dinucleotide") +   
  ylab("% of dinucleotides\nthat are methylated") +   
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.major.y = element_blank()) +  
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill = guide_legend(override.aes = list(size = NULL))) +
  labs(title = "6mA methylation in centromeres") 

g2_centromeres <- g + scale_fill_Publication() + scale_colour_Publication()
g2_centromeres <- g2_centromeres + guides(size = FALSE) + scale_y_continuous(limits = c(0, 20))
g2_centromeres
#---------------------------------
dinucleotides = c('AA', 'AC', 'AG', 'AT')

infected_leaves_freqs = c(1.1316,
                          1.1966,
                          1.5530,
                          1.0806)

germinated_spores_freqs = c(1.0315,
                            1.1402,
                            1.4225,
                            0.9947)

df <- data.frame(infected_leaves_freqs, germinated_spores_freqs,
                 dinucleotides=c(dinucleotides))
names(df) <- c("Infected leaves", "Germinated spores", "Dinucleotide")

head(df)
mdf <- melt(df)
head(mdf)
#---------------------------------
g <- ggplot(mdf, aes(x = Dinucleotide, y = value, fill=variable)) +
  geom_bar(stat='identity', position = "dodge") +
  xlab("Dinucleotide") +   
  ylab("% of dinucleotides\nthat are methylated") +   
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.major.y = element_blank()) +  
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill = guide_legend(override.aes = list(size = NULL))) +
  labs(title = "6mA methylation outside the centromeres") 

g2_not_centromeres <- g + scale_fill_Publication() + scale_colour_Publication()
g2_not_centromeres <- g2_not_centromeres + guides(size = FALSE) + scale_y_continuous(limits = c(0, 20))
g2_not_centromeres
#---------------------------------
require(gridExtra)
png("6mA_DinucleotideFrequencies.png", height = 5, width = 15, units = 'in', res = 300)
grid.arrange(g2_centromeres, g2_not_centromeres, ncol=2)
dev.off()
#---------------------------------
#---------------------------------
#---------------------------------
#---------------------------------
# CHG (H = A, C, or T)
# CHH (H = A, C, or T)
trinucleotides = c('CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CTA', 'CTC', 'CTG', 'CTT')

infected_leaves_freqs = c(0.7346,
                          1.0176,
                          0.6365,
                          1.1497,
                          0.7905,
                          1.4750,
                          8.1340,
                          1.0482,
                          0.8073,
                          2.0090,
                          0.3846,
                          0.8528)

germinated_spores_freqs = c(0.5974,
                            0.8204,
                            0.5058,
                            0.9776,
                            0.6686,
                            1.1463,
                            7.2052,
                            0.8054,
                            0.7293,
                            1.7022,
                            0.2993,
                            0.6997)

df <- data.frame(infected_leaves_freqs, germinated_spores_freqs,
                 trinucleotides=c(trinucleotides))
names(df) <- c("Infected leaves", "Germinated spores", "Trinucleotide")

head(df)
mdf <- melt(df)
head(mdf)
#---------------------------------
g <- ggplot(mdf, aes(x = Trinucleotide, y = value, fill=variable)) +
  geom_bar(stat='identity', position = "dodge") +
  xlab("Trinucleotide") +   
  ylab("% of trinucleotides\nthat are methylated") +   
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.major.y = element_blank()) +  
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill = guide_legend(override.aes = list(size = NULL))) +
  labs(title = "5mC methylation in centromeres") 

g3_centromeres <- g + scale_fill_Publication() + scale_y_continuous(limits = c(0, 20))
g3_centromeres <- g3_centromeres + guides(size = FALSE) 
g3_centromeres
#---------------------------------
#---------------------------------
# CHG (H = A, C, or T)
# CHH (H = A, C, or T)
trinucleotides = c('CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CTA', 'CTC', 'CTG', 'CTT')

infected_leaves_freqs = c(0.6468,
                          0.9722,
                          0.4951,
                          0.9177,
                          0.7457,
                          1.5259,
                          3.3746,
                          0.9205,
                          0.7127,
                          1.1218,
                          0.3527,
                          0.8428
)

germinated_spores_freqs = c(0.5904,
                            0.8630,
                            0.4641,
                            0.8642,
                            0.6680,
                            1.3198,
                            3.4319,
                            0.7844,
                            0.6928,
                            1.0369,
                            0.3196,
                            0.7792)

df <- data.frame(infected_leaves_freqs, germinated_spores_freqs,
                 trinucleotides=c(trinucleotides))
names(df) <- c("Infected leaves", "Germinated spores", "Trinucleotide")

head(df)
mdf <- melt(df)
head(mdf)
#---------------------------------
g <- ggplot(mdf, aes(x = Trinucleotide, y = value, fill=variable)) +
  geom_bar(stat='identity', position = "dodge") +
  xlab("Trinucleotide") +   
  ylab("% of trinucleotides\nthat are methylated") +   
  theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.major = element_blank(),panel.grid.major.y = element_blank()) +  
  theme(legend.title=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") +
  guides(fill = guide_legend(override.aes = list(colour = NULL))) +
  guides(fill = guide_legend(override.aes = list(size = NULL))) +
  labs(title = "5mC methylation outside the centromeres") 

g3_not_centromeres <- g + scale_fill_Publication() + scale_y_continuous(limits = c(0, 20))
g3_not_centromeres <- g3_not_centromeres + guides(size = FALSE) 
g3_not_centromeres
#---------------------------------
require(gridExtra)
png("5mC_TrinucleotideFrequencies.png", height = 5, width = 15, units = 'in', res = 300)
grid.arrange(g3_centromeres, g3_not_centromeres, ncol=2)
dev.off()
