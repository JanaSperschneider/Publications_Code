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
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
########################################################################################
setwd("H:/ANU/InProgress/Project_Pgt_smRNA_V2/Figure1")
########################################################################################
########################################################################################
tpms = read.delim(file = "Tx2_Abundance.txt", 
                  sep="\t", header=TRUE, row.names = 1)
head(tpms)
########################################################################################
gene_names <- c("PGT21_021399 (Argonaute A - chr14A)","PGT21_022388 (Argonaute A - chr14B)",
                "PGT21_001976 (Argonaute B - chr13A)","PGT21_002123 (Argonaute B - chr13B)",
                "PGT21_033256 (Dicer A - chr4A)","PGT21_033881 (Dicer A - chr4B)",
                "PGT21_033709 (Dicer B - chr4B)","PGT21_033021 (Dicer B - chr4A)",       
                "PGT21_028061 (Dicer C - chr6B)","PGT21_029367 (Dicer C - chr6A)",
                "PGT21_001684 (RdRP A - chr10B)","PGT21_002642 (RdRP A - chr10A)",        
                "PGT21_009102 (RdRP B - chr15B)","PGT21_009430 (RdRP B - chr15A)",
                "PGT21_011158 (RdRP C - chr14B)","PGT21_009651 (RdRP C - chr14A)",
                "PGT21_032301 (RdRP D - chr4B)","PGT21_031631 (RdRP D - chr4A",
                "PGT21_035256 (RdRP E - chr16B","PGT21_031875 (RdRP E - chr8A")
########################################################################################                   
include_list <- c("PGT21_021399-T1", "PGT21_022388-T1", "PGT21_001976-T1", "PGT21_002123-T1",
                  "PGT21_033256-T2", "PGT21_033881-T2", "PGT21_033709-T1", "PGT21_033021-T1",
                  "PGT21_028061-T1", "PGT21_029367-T1", "PGT21_001684-T1", "PGT21_002642-T1",
                  "PGT21_009102-T1", "PGT21_009430-T1", "PGT21_011158-T1", "PGT21_009651-T1",
                  "PGT21_032301-T1", "PGT21_031631-T1", "PGT21_035256-T1", "PGT21_031875-T1")
newdata <- tpms[include_list, ]
########################################################################################
dat_average <- data.frame(log1p(rowMeans(subset(newdata, select = c(V1, V2, V3)), na.rm = TRUE)),
                          log1p(rowMeans(subset(newdata, select = c(V4, V5, V6)), na.rm = TRUE)),
                          log1p(rowMeans(subset(newdata, select = c(V7, V8, V9)), na.rm = TRUE)),
                          log1p(rowMeans(subset(newdata, select = c(V10, V11, V12)), na.rm = TRUE)),
                          log1p(rowMeans(subset(newdata, select = c(V13, V14, V15)), na.rm = TRUE)),
                          log1p(rowMeans(subset(newdata, select = c(V16, V17, V18)), na.rm = TRUE)),
                          log1p(rowMeans(subset(newdata, select = c(V19, V20, V21)), na.rm = TRUE)),
                          log1p(rowMeans(subset(newdata, select = c(V22, V23, V24)), na.rm = TRUE)))
rownames(dat_average) <- gene_names
head(dat_average)

breaksList = seq(0, 5.9, by = 0.5)
colors = c("#91BFDB","#FEF0D9","#FDE0B8","#FDD19B","#FDC48D",
           "#FCB27C","#FC9964","#F88254","#F16C4B","#E65139",
           "#D93422","#BA1A10","#990000")

annot <- c("Haustoria", "Germinated spores", "2dpi", "3dpi", "4dpi", "5dpi", "6dpi", "7dpi")

annot_stages <- data.frame(c("Haustorial tissue", "Germinated spores", "Haustoria development", "Haustoria development", "Haustoria development", 
                             "Haustoria development", "Start of sporulation", "Start of sporulation"))
rownames(annot_stages) = annot
names(annot_stages) <- c("Infection stage")
annot_stages
names(dat_average) <- c("Haustoria", "Germinated spores", "2dpi", "3dpi", "4dpi", "5dpi", "6dpi", "7dpi")

my_colour = c("#a6611a", "#80cdc1", "#018571", "#018571", "#018571", "#dfc27d", "#dfc27d", "#dfc27d")


pheatmap(dat_average, cluster_rows = TRUE, cluster_cols=FALSE, gaps_col = c(1,2,3,4,5,6,7), 
         annotation_col = annot_stages,
         fontsize = 24,
         color = colors, 
         breaks = breaksList)


png("RNAi_Expression_TMPs.png", height = 8, width = 15, units = 'in', res = 300)
pheatmap(dat_average, cluster_rows = TRUE, cluster_cols=FALSE, gaps_col = c(1,2,3,4,5,6,7), 
         annotation_col = annot_stages,
         fontsize = 18,
         color = colors, 
         breaks = breaksList)
dev.off()

########################################################################################
# Simplified for slides
########################################################################################
gene_names <- c("Argonaute A - chr14A","Argonaute A - chr14B",
                "Argonaute B - chr13A","Argonaute B - chr13B",
                "Dicer A - chr4A","Dicer A - chr4B",
                "Dicer B - chr4B","Dicer B - chr4A",       
                "Dicer C - chr6B","Dicer C - chr6A",
                "RdRP A - chr10B","RdRP A - chr10A",        
                "RdRP B - chr15B","RdRP B - chr15A",
                "RdRP C - chr14B","RdRP C - chr14A",
                "RdRP D - chr4B","RdRP D - chr4A",
                "RdRP E - chr16B","RdRP E - chr8A")

rownames(dat_average) <- gene_names

png("RNAi_Expression_TMPs_Simple.png", height = 8, width = 15, units = 'in', res = 300)
pheatmap(dat_average, cluster_rows = TRUE, cluster_cols=FALSE, gaps_col = c(1,2,3,4,5,6,7), 
         annotation_col = annot_stages,
         fontsize = 18,
         color = colors, 
         breaks = breaksList)
dev.off()
