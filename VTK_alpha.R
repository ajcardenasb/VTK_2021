library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(GUniFrac)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK_2021/")

asv=read.table("Input_files/VTK2021_ASV_table.txt", header = TRUE, row.names = 1)[,1:20]
colnames(asv)=gsub("\\.", "-", colnames(asv))
map=read.table("Input_files/VTK_metadata.txt", header = T, row.names = 1, sep = "\t")
map$Temperature2=paste(map$Temperature, "ºC", sep = "")

###rarefying
cnts=t(asv[, 1:20])
min(rowSums(cnts)) # determine sample with lowest counts
asv.rar=Rarefy(cnts, 11180)$otu.tab.rff

## alpha diversity
alpha=as.data.frame(t(estimateR(asv.rar,  smallsample = TRUE)))
alpha$Shannon=diversity(asv.rar, index = "shannon")#$shannon
alpha$Evenness=(alpha$Shannon/log10(alpha$S.obs)) ### adding eveness.  J = H'/ln(S) where H' is Shannon Weiner diversity and S is the total number of species 
alpha$Temperature=map$Temperature2[match(rownames(alpha),rownames(map))]


##################################################
##################### Stats ######################
##################################################

shapiro.test(alpha$S.obs) # p-value > 0.05 implying we can assume the normality.

#ANOVAs
summary(aov(alpha$S.obs ~ alpha$Temperature))

#Pairwise t-test
pairwise.t.test(alpha$S.obs,alpha$Temperature, p.adj = "fdr")
#Calculate paiwise for ACE, Shannon and Evenness

## Plot with significance groups
#pdf("./outputs/ASVs_skeleton16S_Shannon.pdf", width=6.5,height=3, pointsize = 12)
P4=c("#2E33D1", "#FFEE32","#D37D47", "#F43535")
ggplot(alpha, aes(x=Temperature, y=S.obs, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "Observed number of ASVs", x="", title = "Observed") + annotate(geom="text", x=1, y=150, label= "A") + annotate(geom="text", x=2, y=185, label= "B")+ annotate(geom="text", x=3, y=295, label= "C") + annotate(geom="text", x=4, y=245, label= "BC") +theme(legend.position = "none")      
#dev.off()
#grid.arrange(a,b,c,d, ncol=2)


