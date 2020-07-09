#---------------------------------
theme_Publication <- function(base_size=14, base_family="helvetica") {
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
            #legend.position = "right",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#a6611a","#dfc27d","#80cdc1","#018571","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#a6611a","#dfc27d","#80cdc1","#018571","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
#---------------------------------
library(ggplot2)
library(scales)
library(gridExtra)
library(reshape2)
library(directlabels)
library(ggpubr)
#---------------------------------
#---------------------------------
samples = c('20 nt','21 nt','22 nt', 'Other lengths')
samples <- as.character(samples)

lengths_noDE <-     c(22, 46.4, 25.9, 100-sum(c(22, 46.4, 25.9)))
lengths_upGS <-     c(29.7, 46.2, 20.4, 100-sum(c(29.7, 46.2, 20.4)))
lengths_up_early <- c(22, 44, 25.9, 100-sum(c(22, 44, 25.9)))
lengths_up_late <-  c(24.8, 24.9, 40.8, 100-sum(c(24.8, 24.9, 40.8)))

df <- data.frame(samples, lengths_upGS, lengths_up_early, lengths_up_late, lengths_noDE)
df

names(df)[names(df) == "lengths_upGS"] <- "Germinated\nspores"
names(df)[names(df) == "lengths_up_early"] <- "Early\ninfection"
names(df)[names(df) == "lengths_up_late"] <- "Late\ninfection"
names(df)[names(df) == "lengths_noDE"] <- "No differential\nexpression"

mdf <- melt(df)
head(mdf)

g <- ggplot(mdf, aes(x = variable, y = value, fill=samples)) +
  geom_bar(stat='identity', width=0.5) +
  xlab("") +   
  ylab("% of sequences") +   
  theme_bw(base_size = 24, base_family = "Helvetica") +
  theme(legend.title=element_blank(),legend.position="bottom") 

g <- g + scale_fill_Publication() + scale_colour_Publication()
g
#---------------------------------
#---------------------------------
png("UpRegulated_smRNACluster_Lengths.png", height = 6, width = 8, units = 'in', res = 300)
g
dev.off()
#---------------------------------
#---------------------------------
L20_noDE <- c(25.7, 8.8, 0, 65.4)
L21_noDE <- c(11.5, 5.9, 0, 82.5)
L22_noDE <- c(16.9, 5, 0, 78.1)

L20_GS   <- c(19.9, 12.5, 2.0, 65.6)
L21_GS   <- c(18.7, 4.9, 1.3, 75)
L22_GS   <- c(18, 2.5, 0.7, 78.7)

L20_early <- c(14.3, 9.8, 2.7, 73.2)
L21_early <- c(15.2, 2.7, 0.4, 81.7)
L22_early <- c(15.2, 2.3, 1.5, 81.1)

L20_late <- c(11.4, 2.2, 0.3, 86)
L21_late <- c(20.4, 2.2, 0.1, 77.3)
L22_late <- c(70.5, 0.4, 0.4, 28.8)
#---------------------------------
bases = c('%A', '%C', '%G', '%U')
bases <- as.character(bases)

df <- data.frame(bases, L22_GS, L22_early, L22_late, L22_noDE)
head(df)
names(df)[names(df)=="L22_GS"] <- "Germinated\nspores"
names(df)[names(df)=="L22_early"] <- "Early\ninfection"
names(df)[names(df)=="L22_late"] <- "Late\ninfection"
names(df)[names(df)=="L22_noDE"] <- "No differential\nexpression"

mdf <- melt(df)
mdf

g <- ggplot(mdf, aes(x = variable, y = value,fill=bases)) +
  geom_bar(stat='identity', width=0.5) +
  xlab("") +   
  ylab("% of sequences") +   
  theme_bw(base_size = 20, base_family = "Helvetica") +
  labs(title = "22nt sRNAs") +
  theme(legend.title=element_blank(),legend.position="bottom") 
g22 <- g + scale_fill_Publication() + scale_colour_Publication()
g22
#---------------------------------
#---------------------------------
df <- data.frame(bases, L21_GS, L21_early, L21_late, L21_noDE)
head(df)
names(df)[names(df)=="L21_GS"] <- "Germinated\nspores"
names(df)[names(df)=="L21_early"] <- "Early\ninfection"
names(df)[names(df)=="L21_late"] <- "Late\ninfection"
names(df)[names(df)=="L21_noDE"] <- "No differential\nexpression"

mdf <- melt(df)
mdf

g <- ggplot(mdf, aes(x = variable, y = value,fill=bases)) +
  geom_bar(stat='identity', width=0.5) +
  xlab("") +   
  ylab("% of sequences") +   
  theme_bw(base_size = 20, base_family = "Helvetica") +
  labs(title = "21nt sRNAs") +
  theme(legend.title=element_blank(),legend.position="bottom") 
g21 <- g + scale_fill_Publication() + scale_colour_Publication()
g21
#---------------------------------
#---------------------------------
df <- data.frame(bases, L20_GS, L20_early, L20_late, L20_noDE)
head(df)
names(df)[names(df)=="L20_GS"] <- "Germinated\nspores"
names(df)[names(df)=="L20_early"] <- "Early\ninfection"
names(df)[names(df)=="L20_late"] <- "Late\ninfection"
names(df)[names(df)=="L20_noDE"] <- "No differential\nexpression"

mdf <- melt(df)
mdf

g <- ggplot(mdf, aes(x = variable, y = value,fill=bases)) +
  geom_bar(stat='identity', width=0.5) +
  xlab("") +   
  ylab("% of sequences") +   
  theme_bw(base_size = 20, base_family = "Helvetica") +
  labs(title = "20nt sRNAs") +
  theme(legend.title=element_blank(),legend.position="bottom") 
g20 <- g + scale_fill_Publication() + scale_colour_Publication()
g20
#---------------------------------
#---------------------------------
png("UpRegulated_5PrimeBases.png", height = 6, width = 18, units = 'in', res = 300)
grid.arrange(g20,g21,g22,ncol=3)
dev.off()
#---------------------------------


