#################################################################################################################################################################
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
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
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
########################################################################################
library(edgeR)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
########################################################################################
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figures_Revision/Figure8")
##########################################################
##########################################################
countsTable = read.delim(file = "Counts_Filtered.txt", 
                         sep="\t", header=TRUE, 
                         row.names=1)
head(countsTable)
countsTable <- countsTable[,-1] 
countsTable <- countsTable[,-1] 
names(countsTable)
# Delete U1 and Uninfected plant samples
countsTable <- countsTable[,-5] 
countsTable <- countsTable[,-5] 
countsTable <- countsTable[,-5] 
countsTable <- countsTable[,-4] 

head(countsTable)
names(countsTable)
########################################################################################
########################################################################################
########################################################################################
# Analyze the rust samples
########################################################################################
########################################################################################
########################################################################################
samples <- c("Spores", "Spores", "Spores", 
             "dpi3", "dpi3", "dpi3", 
             "dpi5", "dpi5", "dpi5", 
             "dpi7", "dpi7", "dpi7")

designFull = model.matrix(~0+factor(samples))
colnames(designFull) <- levels(factor(samples))
designFull

d <- DGEList(counts=countsTable,group=factor(samples))
d <- calcNormFactors(d)
d$samples

plotMDS(d, col=rep(1:2, each=3))

pcaData <- plotMDS(d, col=rep(1:2, each=3))
x <- pcaData$x
y <- pcaData$y

pcaData <- data.frame(samples,x,y)

g1 <- ggplot(pcaData, aes(x=x, y=y)) + 
  geom_point(size=3)+
  xlab("Leading logFC dim1") + 
  ylab("Leading logFC dim2") + 
  coord_fixed() + 
  scale_shape_manual(values=1:4) + 
  geom_point(aes(fill=samples),colour="black",pch=21,size=5)+
  scale_fill_brewer("Sample",palette="Spectral")

g1 <- g1 + scale_fill_Publication() + theme_Publication() + scale_fill_discrete(name = "")

g1

png("MDS_Plot_smRNA.png", height=6, width=6, units="in", res=300)
g1
dev.off()

########################################################################################
d <- estimateDisp(d,designFull,robust=TRUE)
fit <- glmQLFit(d,designFull,robust=TRUE)
########################################################################################
my.contrasts <- makeContrasts(dpi3vsSpores = dpi3-Spores, 
                              dpi5vsSpores = dpi5-Spores, 
                              dpi7vsSpores = dpi7-Spores,
                              dpi7vs3dpi = dpi7-dpi3,
                              dpi7vs5dpi = dpi7-dpi5,
                              dpi3vs5dpi = dpi3-dpi5,
                              dpi3vs7dpi = dpi3-dpi7,
                              dpi5vs3dpi = dpi5-dpi3,
                              dpi5vs7dpi = dpi5-dpi7,
                              levels=designFull)

# Compare 3dpi versus spores
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi3vsSpores"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter, file="EdgeR_3dpivsSpores.csv")

# Compare 5dpi versus spores
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi5vsSpores"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter, file="EdgeR_5dpivsSpores.csv")

# Compare 7dpi versus spores
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi7vsSpores"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter, file="EdgeR_7dpivsSpores.csv")
########################################################################################
# Compare 7dpi versus 3dpi
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi7vs3dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter, file="EdgeR_7dpivs3dpi.csv")

# Compare 7dpi versus 5dpi
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi7vs5dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter, file="EdgeR_7dpivs5dpi.csv")

# Compare 3dpi versus 5dpi
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi3vs5dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter, file="EdgeR_3dpivs5dpi.csv")

# Compare 3dpi versus 7dpi
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi3vs7dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter_dpi3vs7dpi = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter_dpi3vs7dpi, file="EdgeR_3dpivsdpi7.csv")

# Compare 5dpi versus 3dpi
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi5vs3dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter, file="EdgeR_5dpivsdpi3.csv")

# Compare 5dpi versus 7dpi
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi5vs7dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05, ]
table(DEtag_all$table[, "FDR"] < 0.05) 
write.csv(DEtag_filter, file="EdgeR_5dpivsdpi7.csv")
########################################################################################
########################################################################################
lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi7vsSpores"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter_dpi7vsSpores_up = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] > 1, ]
DEtag_filter_dpi7vsSpores_down = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] < -1, ]

lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi7vs3dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter_dpi7vsdpi3_up = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] > 1, ]
DEtag_filter_dpi7vsdpi3_down = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] < -1, ]

lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi5vs3dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter_dpi5vsdpi3_up = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] > 1, ]
DEtag_filter_dpi5vsdpi3_down = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] < -1, ]

lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi5vs7dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter_dpi5vsdpi7_up = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] > 1, ]
DEtag_filter_dpi5vsdpi7_down = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] < -1, ]

lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi7vs5dpi"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter_dpi7vsdpi5_up = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] > 1, ]
DEtag_filter_dpi7vsdpi5_down = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] < -1, ]

lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi3vsSpores"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter_dpi3vsSpores_up = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] > 1, ]
DEtag_filter_dpi3vsSpores_down = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] < -1, ]

lrt <- glmQLFTest(fit,contrast=my.contrasts[,"dpi5vsSpores"])
topTags(lrt)
DEtag_all = topTags(lrt, n = dim(d)[1])  
DEtag_filter_dpi5vsSpores_up = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] > 1, ]
DEtag_filter_dpi5vsSpores_down = DEtag_all$table[DEtag_all$table[, "FDR"] < 0.05 & DEtag_all$table[, "logFC"] < -1, ]

########################################################################################
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list("Up-regulated in 7dpi vs. germinated spores" = rownames(DEtag_filter_dpi7vsSpores_up),
           "Up-regulated in 7dpi vs. 5dpi" = rownames(DEtag_filter_dpi7vsdpi5_up),
           "Up-regulated in 7dpi vs. 3dpi" = rownames(DEtag_filter_dpi7vsdpi3_up)
  ),
  category.names = c("Up-regulated in 7dpi\nvs. germinated spores" , "Up-regulated in 7dpi\nvs. 5dpi" , "Up-regulated in 7dpi vs. 3dpi"),
  filename = 'Up7dpi_venn_diagramm.png',
  output=TRUE,
  # Output features
  imagetype="png", height = 1000 , width = 1000 , resolution = 300, compression = "lzw",
  # Circles
  lwd = 0.2, lty = 'blank', fill = myCol,
  # Numbers
  cex = .8, fontface = "bold", fontfamily = "sans",
  # Set names
  cat.cex = 0.8, cat.fontface = "bold", cat.default.pos = "outer",
  cat.pos = c(-18, 18, 180))

venn.diagram(
  x = list("Up-regulated in germinated spores vs. 7dpi" = rownames(DEtag_filter_dpi7vsSpores_down),
           "Up-regulated in germinated spores vs. 5dpi" = rownames(DEtag_filter_dpi5vsSpores_down),
           "Up-regulated in germinated spores vs. 3dpi" = rownames(DEtag_filter_dpi3vsSpores_down)
  ),
  category.names = c("Up-regulated in germinated\nspores vs. 7dpi" , "Up-regulated in germinated\nspores vs. 5dpi" , 
                     "Up-regulated in germinated\nspores vs. 3dpi"),
  filename = 'UpGS_venn_diagramm.png',
  output=TRUE,
  # Output features
  imagetype="png", height = 1000 , width = 1000 , resolution = 300, compression = "lzw",
  # Circles
  lwd = 0.2, lty = 'blank', fill = myCol,
  # Numbers
  cex = .8, fontface = "bold", fontfamily = "sans",
  # Set names
  cat.cex = 0.8, cat.fontface = "bold", cat.default.pos = "outer",
  cat.pos = c(-18, 18, 180))

########################################################################################

