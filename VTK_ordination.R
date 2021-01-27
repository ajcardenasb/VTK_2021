setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK_2021/")
library(phyloseq)
library(ggplot2)
library(microbiome)
library(vegan)

asv=read.table("Input_files/VTK2021_ASV_table.txt", header = TRUE, row.names = 1)[,1:20]
map=read.table("Input_files/VTK_metadata.txt", header = T, row.names = 1, sep = "\t")
map$Temperature2=paste(map$Temperature, "ºC", sep = "")
tax=read.table("Input_files/VTK2021_ASV_table.txt", header = TRUE, row.names = 1)[,22:27]
colnames(asv)=gsub("\\.", "-", colnames(asv))

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

##transform data and subset phyloseq objects
phy.t=microbiome::transform(phy.all, transform = "compositional", target = "OTU", shift = 0, scale = 1) # try clr and log10 transformations

P4=c("#2E33D1", "#FFEE32","#D37D47", "#F43535") #27, 29, 32, 34
ordi = ordinate(phy.t, method = "PCoA", distance = "bray") #distance : try bray and euclidean. method: try PCoA, RDA, NMDS
plot_ordination(phy.t,ordi, color = "Temperature2")  + geom_point(size = 3, alpha = 1) + theme_bw()   + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4) + scale_shape_manual(values=c(0, 15)) + labs(title = "Compositional + Bray + PCoA")

## Assigment 1: Use the function grid.arrange to make a comparative panel with
# 1. compositional + PCoA + Bray
# 2. clr + RDA + euclidean
# 3. log10 + RDA + Euclidean
# 4. log10 + RDA + bray
# 5. compositional + NMDS + Bray
# 6. clr + NMDS + euclidean

#### Permanovas
source("VTK_pairwiseAdonis_function.R")
asv.n=as.data.frame(t(sweep(asv[, 1:20],2,colSums(asv[, 1:20]),`/`)))
asv.n$Temperature=map$Temperature2[match(rownames(asv.n), rownames(map))]

##overall 
adonis(asv.n[,1:1000] ~ asv.n$Temperature ) 

##pairwise comparisons
pairwise.adonis(asv.n[,1:1000], asv.n$Temperature, p.adjust.m = "fdr")

## Assigment 2: run permanovas to compare treatmens ONLY using ASVs classified as Endozoicomonadaceae
