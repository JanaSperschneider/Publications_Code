#################################################################################################################################################################
library(ggplot2)
library(drawProteins)
#################################################################################################################################################################
prot_data <- data.frame(
  type <- c('CHAIN', 'DOMAIN', 'DOMAIN', 'DOMAIN', 'REPEAT',
            'CHAIN', 'DOMAIN', 'DOMAIN', 'DOMAIN', 'REPEAT',
            'CHAIN', 'DOMAIN', 'DOMAIN', 'DOMAIN',
            'CHAIN', 'DOMAIN', 'DOMAIN', 'DOMAIN'),
  description <- c('', 'N-terminal\n domain of\n argonaute' , 'Argonaute\n linker 2\ndomain', 'Piwi domain', 'DEDX',
                   '', 'N-terminal\n domain of\n argonaute' , 'Argonaute\n linker 2\ndomain', 'Piwi domain', 'DEDX',
                   '', 'N-terminal\n domain of\n argonaute' , 'Argonaute\n linker 1\ndomain', 'Piwi domain',
                   '', 'N-terminal\n domain of\n argonaute' , 'Argonaute\n linker 1\ndomain', 'Piwi domain'),
  begin <- c(1, 13, 378, 528, 652,
             1, 18, 383, 533, 657,
             1, 16, 170, 523,
             1, 16, 170, 523),
  end <- c(810, 162, 425, 809, 656,
           815, 167, 430, 814, 660,
           882, 160, 225, 831,
           882, 160, 225, 831),
  length <- c(810, 162-13+1, 425-378+1, 809-528+1, 4, 
              815, 167-19+1, 430-383+1, 814-533+1, 4,
              882, 160-16+1, 225-170+1, 831-523+1,
              882, 160-16+1, 225-170+1, 831-523+1 ),
  accession <- c('AGO1A', 'AGO1A', 'AGO1A', 'AGO1A', 'AGO1A',
                 'AGO1B', 'AGO1B', 'AGO1B', 'AGO1B', 'AGO1B',
                 'AGO2A', 'AGO2A', 'AGO2A', 'AGO2A',
                 'AGO2B', 'AGO2B', 'AGO2B', 'AGO2B'),
  entryName <- c('AGO1A', 'AGO1A', 'AGO1A', 'AGO1A', 'AGO1A',
                 'AGO1B', 'AGO1B', 'AGO1B', 'AGO1B', 'AGO1B',
                 'AGO2A', 'AGO2A', 'AGO2A', 'AGO2A',
                 'AGO2B', 'AGO2B', 'AGO2B', 'AGO2B'),
  taxid <- c(321614, 321614, 321614, 321614, 321614,
             321614, 321614, 321614, 321614, 321614,
             321614, 321614, 321614, 321614,
             321614, 321614, 321614, 321614),
  order <- c(1, 1, 1, 1, 1,
             2, 2, 2, 2, 2,
             3, 3, 3, 3,
             4, 4, 4, 4)
)

names(prot_data) <- c('type', 'description', 'begin', 'end', 'length',
                      'accession', 'entryName', 'taxid', 'order')

prot_data

p <- draw_canvas(prot_data)
p <- draw_chains(p, prot_data, label_size=6, fill="lightgray")
p <- draw_repeat(p, prot_data, label_size=6, fill="blue", outline="black", label_repeats=FALSE) #+xlim(-120,1600) + ylim(-5, 37) 
p <- draw_domains(p, prot_data, label_size=6, label_domains=TRUE)
p

# background and y-axis
p <- p + theme_bw(base_size = 20) + # white backgnd & change text size
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank(), title = element_blank()) 
p

# move legend to top
p <- p + theme(legend.position="none") + labs(fill="") + 
  annotate("text", x = 654, y = 1.3, label = "DEDG", size = 6) +
  annotate("text", x = 658, y = 2.3, label = "DEDG", size = 6) 

 
p

#################################################################################################################################################################
png("Figure12B.png", height = 10, width = 20, units = 'in', res = 300)
p
dev.off()
#################################################################################################################################################################



