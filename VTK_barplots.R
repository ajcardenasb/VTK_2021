library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(gridExtra)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK_2021/")

#############################################################
#####Taxonomic profiles of the 10 most abundant families#####
#############################################################
#stats=read.table("Input_files/VTK2021_ASV_stats.txt", header = TRUE, row.names = 1, sep = " ")
## load data in
asv=read.table("Input_files/VTK2021_ASV_table.txt", header = TRUE, row.names = 1, sep = " ")
map=read.table("Input_files/VTK_metadata.txt", header = T, row.names = 1, sep = "\t")
## fix sample names 
colnames(asv)=gsub("\\.", "-", colnames(asv)) # this replaces dots by "-"

## aggregate 
names(asv) # identify my samples and the family column
asv.tax.ag=aggregate(asv[, 1:20], by = list(asv[, 26]), FUN =  sum) #define sample range and group factor
topFamilies=asv.tax.ag[order(rowSums(asv.tax.ag[,2:ncol(asv.tax.ag) ]),decreasing = TRUE),][1:10,1]
asv.tax.ag$Group.1=ifelse(asv.tax.ag$Group.1 %in% topFamilies, as.character(asv.tax.ag$Group.1), "zOthers")
asv.gg=aggregate(asv.tax.ag[, 2:ncol(asv.tax.ag)], by = list(asv.tax.ag$Group.1), FUN =  sum)
all.l=reshape2::melt(asv.gg, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")

## Add sample information
all.l$Tempearature=map$Temperature[match(all.l$Sample, rownames(map))]
all.l$Genotype=map$Genotype[match(all.l$Sample, rownames(map))]

## Plot
P10=c("#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#ffcc00", "#ADADAD")

ggplot() + geom_bar(aes(y = Abundance, x = Genotype, fill = Family), data = all.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme_classic() +theme(axis.text.x=element_text(angle=90,hjust=1), legend.position = "right",  plot.title = element_text(hjust = 0.5), panel.spacing = unit(1.5, "lines")) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + scale_fill_manual(values=P10) + facet_grid(Temperature~.) 

### Assignment 1: plot the 10 most abundant genera

### Assignment 2: plot the 10 most abundant ASVs. Add taxonomy (ASV - Genus - Family) to the labels
