library(gggenes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)


#import prophages/defence islands regions
t <-read.csv(file = 'merged_table.csv',quote ="",sep = ",",
             row.names = NULL, 
             stringsAsFactors =FALSE,header=TRUE)
# reorder the columns
GenesDF2 <- t[, c(1, 3, 4, 5,6,7,8)]
##rename columns
names(GenesDF2)[names(GenesDF2) == "strand"] <- "orientation"



#write out to incorporate VFDB annotations
write.table(finaldf, file="temprnew.txt", quote=F, row.names=F,sep="\t")

newfinaldf <-read.csv(file = 'temprnew.txt',quote ="",sep = "\t",
                      row.names = NULL, 
                      stringsAsFactors =FALSE,header=TRUE)
#plot
g <- ggplot(GenesDF2, aes(xmin = start, xmax = end, y = molecule , fill=system,forward = orientation)) +
  geom_gene_arrow(arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + 
  scale_fill_manual(values = c("SoFic" = "darkgreen",
                               "PDC-M11" = "thistle1",
                               "PD-Lambda-5" = "brown",
                               "RM_type_II"="midnightblue",
                               "PDC-S20" ="darkorange",
                               "pycsar_effector" = "yellow",
                               "PDC-M53"= "mediumspringgreen",
                               "phage late control D family protein "="red",
                               "HTH"="cyan",
                               "AntiRM"= "peachpuff1"))+
  theme_genes() %+replace% theme(axis.title.y = element_blank(),
                                 axis.ticks.x = element_blank(),
                                 legend.direction="vertical",
                                 legend.position="right",
                                 panel.grid.major.y =element_blank())+scale_y_discrete(expand = expand_scale(mult = c(1, 6)))

ggsave("genes.pdf", width = 30, height = 30, units = "cm",limitsize = FALSE)
