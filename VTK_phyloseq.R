############################################################
##################### phyloseq #############################
############################################################
library("phyloseq")
library(ggplot2)
library(reshape2)
library(microbiome)

setwd("/Users/anny/Documents/GitHub/RSS_16Shttps://github.com/ajcardenasb/VTK_2021/")

otu.2=read.table("./inputFiles/RSS_OTU_table.txt", header = TRUE)
colnames(otu.2)=gsub('\\.', '-', colnames(otu.2))
map=read.table("./inputFiles/RSS_metadata.txt", header = TRUE, row.names = 1, sep ='\t')
sha=read.table("shape_mapping.txt", header = FALSE,  sep ='\t')
map$shape1=paste(map$Time,"-",map$Genotype, sep = "")

tax = read.csv("./inputFiles/RSS.final.taxonomy", header = TRUE, sep = "\t", row.names = 1)
tax$Taxonomy = gsub("[(0-9)]", "", tax$Taxonomy) 
tax$Taxonomy = gsub("[a-z]__", "", tax$Taxonomy) 
tax[c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')] <- colsplit(tax$Taxonomy,';',c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'))
tax=tax[, -c(3)]

otu.t= otu_table(otu.2, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

##transform data and subset phyloseq objects
phy.t=microbiome::transform(phy.all, transform = "clr", target = "OTU", shift = 0, scale = 1)
cb=subset_samples(phy.t, Experiment=="CB")
rss=subset_samples(phy.t, Experiment=="RSS")

P4=c("#2E33D1", "#FFEE32","#D37D47", "#F43535") #27, 29, 32, 34
cb_PCOA_br = ordinate(cb, method = "RDA", distance = "euclidean")
c.p=plot_ordination(cb,cb_PCOA_br, color = "Temperature", shape = "shape1")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4) + scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6))

#plot_ordination(cb,cb_PCOA_br, color = "Genotype", shape = "Temperature")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS") + theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6))


rss_PCOA_br = ordinate(rss, method = "PCoA", distance = "euclidean")
r.p=plot_ordination(rss,rss_PCOA_br, color = "Temperature", shape = "shape1")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("RSS") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4) + scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6))

pdf("CBASS_vs_RSS_ordination.pdf", onefile = TRUE, width=15,height=15)
gridExtra::grid.arrange(c.p, r.p,  ncol=2)
dev.off()

write.table(cb_PCOA_br$vectors[,c(1,2)], "PCOA_bact_CBASS.txt", quote= FALSE, row.names = TRUE, sep = "\t" )
write.table(rss_PCOA_br$vectors[,c(1,2)], "PCOA_bact_RSS.txt", quote= FALSE, row.names = TRUE, sep = "\t" )

##### 4 ordination plots

cb.t1=subset_samples(phy.all, Experiment=="CB" & Time =="T1")
cb.t3=subset_samples(phy.all, Experiment=="CB" & Time =="T3")
rss.t1=subset_samples(phy.all, Experiment=="RSS" & Time =="T1")
rss.t3=subset_samples(phy.all, Experiment=="RSS" & Time =="T3")


cb1_PCOA_br = ordinate(cb.t1, method = "PCoA", distance = "bray")
plot_ordination(cb.t1,cb1_PCOA_br, color = "Temperature", shape = "shape1")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS T1") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4) + scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6)) #values=c(18, 5, 25, 6, 15, 0, 16, 1, 17, 2)

cb3_PCOA_br = ordinate(cb.t3, method = "PCoA", distance = "bray")
b=plot_ordination(cb.t3,cb3_PCOA_br, color = "Temperature", shape = "Genotype")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS T3") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4)

rss1_PCOA_br = ordinate(rss.t1, method = "PCoA", distance = "bray")
c=plot_ordination(rss.t1,rss1_PCOA_br, color = "Temperature", shape = "Genotype")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("RSS T1") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4)

rss3_PCOA_br = ordinate(rss.t3, method = "PCoA", distance = "bray")
d=plot_ordination(rss.t3,rss3_PCOA_br, color = "Temperature", shape = "Genotype")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("RSS T3") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4)

gridExtra::grid.arrange(a, b,c,d,  ncol=2)

#### No outliers and Edoz

#original
no.f=subset_samples(phy.all, Experiment=="CB" |   Experiment=="RSS" )
no.f_PCOA_br = ordinate(no.f, method = "PCoA", distance = "bray")
plot_ordination(no.f,no.f_PCOA_br , color = "Temperature", shape = "shape1")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS") + theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6)) + facet_grid(~Experiment) + scale_colour_manual(values=P4)

#no outliers
no.o=subset_samples(no.f, !Genotype=="G13A" &   !Genotype=="G8A" )
PCOA_br = ordinate(no.o, method = "PCoA", distance = "bray")
plot_ordination(no.o,PCOA_br , color = "Temperature", shape = "shape1")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS") + theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6)) + facet_grid(~Experiment) + scale_colour_manual(values=P4)


#only endos
endos=subset_taxa(no.f, Family=="Endozoicimonaceae")
endos_PCOA_br = ordinate(endos, method = "PCoA", distance = "bray")
plot_ordination(endos,endos_PCOA_br , color = "Temperature", shape = "shape1")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS") + theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6)) + facet_grid(~Experiment) + scale_colour_manual(values=P4)

#only endos no outliers
endos_noOut=subset_taxa(no.o, Family=="Endozoicimonaceae")
endos_noOut_PCOA_br = ordinate(endos_noOut, method = "PCoA", distance = "bray")
plot_ordination(endos_noOut,endos_noOut_PCOA_br  , color = "Temperature", shape = "shape1")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS") + theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6)) + facet_grid(~Experiment) + scale_colour_manual(values=P4)

# Treatment effect on T1 vs T2
library(plyr)

val=as.data.frame(no.f_PCOA_br$vectors[,1])
colnames(val) = c("EigenValues")
val$Temperature=map$Temperature[match(rownames(val), rownames(map))]
val$Experiment=map$Experiment[match(rownames(val), rownames(map))]
val$Time=map$Time[match(rownames(val), rownames(map))]
val$hybrid=paste(val$Temperature,val$Time, sep = "-")
sum=ddply(val, c("Experiment", "hybrid"), summarise, N = length(EigenValues), mean = mean(EigenValues), sd= sd(EigenValues), se= sd / sqrt(N))
sum$lower=sum$mean-sum$sd
sum$upper=sum$mean+sum$sd
sum$color=gsub("[0-9][0-9]-", "", sum$hybrid)

ggplot(sum, aes(x=hybrid, weight=mean, ymin=lower, ymax=upper, fill=color)) + geom_boxplot(position=position_dodge(), aes(y=mean), stat="identity", alpha=0.7) + geom_errorbar(position=position_dodge(width=0.5), colour="black")  + labs( y= "PCoA1 (38.7%)", x="", title= "")+ theme(legend.title = element_blank()) + scale_fill_manual(values=c("#508CA4", "#FF5714")) + facet_grid(~Experiment)

ggplot(val, aes(x=hybrid, y=mean, fill=Time)) + geom_boxplot()   + labs( y= "PCoA1 (38.7%)", x="", title= "") + theme(legend.title = element_blank()) + scale_fill_manual(values=c("#508CA4", "#FF5714")) + facet_grid(~Experiment)

ggplot(val, aes(x=hybrid, y=EigenValues, fill=Time)) + geom_boxplot() + facet_grid(~Experiment) +  scale_fill_manual(values=c("#E01811", "#43C62B"))

############################################################
##################### Alpha-diversity ######################
############################################################
library(plyr)

alpha=estimate_richness(phy.all, split = TRUE, measures = c("Observed", "Chao1"))
alpha.2=alpha[,-3]
message("Average OTU per sample: ", mean(subset(alpha, !rownames(alpha) %like% "field")$Observed))
#write.table(alpha.2, "alpha_div.txt", row.names = TRUE, quote=FALSE, sep= "\t")

# rownames(alpha)=gsub("\\.", "-", rownames(alpha))
# alpha.met=merge(alpha, met, by = "row.names")
# alpha.sum=ddply(alpha.met, c("Experiment", "Time", "Temperature"), summarise, N = length(Observed), mean = mean(Observed), sd= sd(Observed), se= sd / sqrt(N))
