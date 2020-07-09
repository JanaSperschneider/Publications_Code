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
library(ggsignif)
library(reshape2)
#---------------------------------
#---------------------------------
#---------------------------------
#---------------------------------
df1 <- read.delim('methylation.freq.centromeres.dpi7.txt', sep="\t", header=FALSE)
df2 <- read.delim('methylation.freq.not_centromeres.dpi7.txt', sep="\t", header=FALSE)
df3 <- read.delim('methylation.freq.centromeres.GS.txt', sep="\t", header=FALSE)
df4 <- read.delim('methylation.freq.not_centromeres.GS.txt', sep="\t", header=FALSE)

df1 <- data.frame(val=df1,type="Centromeres\n(Infected leaves)")
df2 <- data.frame(val=df2,type="Non-centromeres\n(Infected leaves)")
df3 <- data.frame(val=df3,type="Centromeres\n(Germinated spores)")
df4 <- data.frame(val=df4,type="Non-centromeres\n(Germinated spores)")

data <- rbind(df1,df2,df3,df4)

g <- ggplot(data, aes(x=type, y=V1, fill=type)) + geom_boxplot() +
  labs(title = "CG methylation frequencies") +
  xlab("") + ylab("Methylation frequency") + 
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_fill_manual(values=c("#8c510a", "#8c510a", "#80cdc1", "#80cdc1")) + coord_flip()

g

g <- g  + theme_bw() +theme(legend.position="none") 
g

png("Methylation_Frequencies.png", height = 5, width = 5, units = 'in', res = 300)
g
dev.off()
