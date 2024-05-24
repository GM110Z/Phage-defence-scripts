library(gggenes)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)


#import prophages/defence islands regions
t <-read.csv(file = 'gggenes-input.csv',quote ="",sep = ",",
               row.names = NULL, 
               stringsAsFactors =FALSE,header=TRUE)
# reorder the columns
GenesDF2 <- t[, c(1,2, 3, 4, 5,6)]
##rename columns
names(GenesDF2)[names(GenesDF2) == "add_column"] <- "molecule"
names(GenesDF2)[names(GenesDF2) == "strand"] <- "orientation"

#process padloc
a <- read.csv(file = 'padloc.txt',quote = "",sep = "\t",
                              stringsAsFactors = FALSE,header=TRUE)

#add padloc Ids to dataset
names(a)[names(a) == "target.name"] <- "Protein_ID"
names(a)[names(a) == "seqid"] <- "molecule"

#join the dataframes
d <-merge(GenesDF2, a, by=c("Protein_ID","molecule"), all.x=T, all.y=T)

#Select the bloody columns now of new dataframe
finaldf <- d[, c(1,2, 3, 4, 5,6,9)]
#write out to incorporate VFDB annotations
write.table(finaldf, file="temprnew.txt", quote=F, row.names=F,sep="\t")

newfinaldf <-read.csv(file = 'temprnew.txt',quote ="",sep = "\t",
                      row.names = NULL, 
                      stringsAsFactors =FALSE,header=TRUE)
#plot
g <- ggplot(newfinaldf, aes(xmin = start.x, xmax = end.x, y = molecule , fill=system,forward = orientation)) +
  geom_gene_arrow(arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) + 
  scale_fill_manual(values = c("Colicin" = "darkgreen",
                               "Colicin immunity" = "darkgreen",
                                     "RM_type_I" = "thistle1",
                                     "Esterase" = "brown",
                                     "LysR regulator" = "#000000",
                                   "PrrC"="midnightblue",
                                   "Integrase"="seashell3",
                                   "Adherence"="#009292",
                                   'cbass_type_IIs'="khaki3",
                                   "Exotoxin"="slategray2",
                                    "TA" ="darkorange",
                               "Pseudogene"="snow",
                               "DRT_class_II"="mistyrose",
                               "RhuM virulence factor" = "yellow",
                               "brex_type_I"= "mediumspringgreen",
                               "BrxU" = "red3",
                               "DRT_class_III" = "plum3",
                               "qatABCD" = "lightcyan",
                               "septu_type_I" = "red",
                               "Transposase"="tomato",
                               "Integrase" = "lightcyan4",
                               "Menshen" = "lightgoldenrod1",
                               "Hydrolase-TM" = "mediumorchid",
                               "BrxR" = "blue",
                               "hachiman_type_I" = "palevioletred1",
                               "RM-like"="peachpuff1",
                               "PD-T4-3" = "cyan",
                               "Pyocin" = "limegreen",
                               "tmn" = "steelblue",
                               "GAO_19" = "violetred4",
                               "retron_XII" = "sienna2"))+
  theme_genes() %+replace% theme(axis.title.y = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  legend.direction="vertical",
                                  legend.position="right",
                                  panel.grid.major.y =element_blank())+scale_y_discrete(expand = expand_scale(mult = c(1, 6)))

ggsave("genes.pdf", width = 30, height = 30, units = "cm",limitsize = FALSE)
