theme_Publication <- function(base_size=14, base_family="helvetica") {
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
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.text = element_text(size = rel(1.2)),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#a6611a","#018571","#018571","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#a6611a","#018571","#018571","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

library(ggplot2)
library(scales)
library(gridExtra)
library(reshape2)
library(directlabels)
#---------------------------------
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Centromere_Figure/")
#---------------------------------
df <- read.delim('Centromere_Repeats.txt', sep="\t")
head(df)
df$Chromosome <- factor(df$Chromosome, levels=unique(df$Chromosome))
#---------------------------------
#---------------------------------
#---------------------------------
library(ggrepel)

df_scatter <- df
head(df_scatter)

g1 <- ggplot(df_scatter, aes(x = LTR.GYPSY, y = LTR.Copia, color=Region)) +
  geom_text_repel(data = subset(df_scatter, Region=="Centromere"), aes(label = Chromosome), color="black") +
  geom_point(size = 3) +
  xlab("% of sequence occupied\nby LTR/Gypsy elements") +   
  ylab("% of sequence occupied\nby LTR/Copia elements") +
  theme_Publication()


#g1 <- g1 +  stat_ellipse()#geom_density_2d()
g1 <- g1 + scale_fill_Publication() + scale_colour_Publication() + 
  scale_color_manual(name="",
                     labels=c("Centromeric region","Non-centromeric region"),
                     values=c("#a6611a","#018571"))

g1
#---------------------------------
g2 <- ggplot(df_scatter, aes(x = LTR.GYPSY, y = DNA.MULE.MuDR, color=Region)) +
  geom_text_repel(data = subset(df_scatter, Region=="Centromere"), aes(label = Chromosome), color="black") +
  geom_point(size = 3) +
  xlab("% of sequence occupied\nby LTR/Gypsy elements") +   
  ylab("% of sequence occupied\nby DNA/MULE-MuDR elements") +     
  theme_Publication()


#g2 <- g2 +  stat_ellipse()#geom_density_2d()
g2 <- g2 + scale_fill_Publication() + scale_colour_Publication() + 
  scale_color_manual(name="",
                     labels=c("Centromeric region","Non-centromeric region"),
                     values=c("#a6611a","#018571"))


g2
#---------------------------------
png("Gypsy_Copia.png", height = 7, width = 10, units = 'in', res = 300)
g1
dev.off()
#---------------------------------
png("Gypsy_DNA.png", height = 7, width = 10, units = 'in', res = 300)
g2
dev.off()
#---------------------------------
#---------------------------------
#---------------------------------
data <- data.frame(df$Region, df$Covered.by.repeats)
head(data)

names(data)[names(data) == "df.Region"] <- "Sequence"
names(data)[names(data) == "df.Covered.by.repeats"] <- "% repetitive sequence"

mdf <- melt(data)
head(mdf)

g <- ggplot(mdf, aes(x = Sequence, y = value, fill=Sequence)) +
  geom_boxplot() +
  xlab("") +   
  ylab("% of repetitive sequence") +   
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(legend.title=element_blank(),legend.position="none") 

g <- g + annotate("text", x = "Centromere", y = 38, label = "chr8A")
g <- g + annotate("text", x = "Centromere", y = 41.5, label = "chr3A")
g <- g + annotate("text", x = "NonCentromere", y = 26, label = "chr17A")

g <- g + scale_fill_Publication() + scale_colour_Publication() + scale_x_discrete(breaks=c("Centromere","NonCentromere"),
                                                                                  labels=c("Centromeric\nregions", "Non-centromeric\nregions"))
g

png("RepeatContent.png", height = 7, width = 10, units = 'in', res = 300)
g
dev.off()
