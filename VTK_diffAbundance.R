library(ANCOMBC)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(compositions)
library(pheatmap)
library(dplyr)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/VTK_2021/")

asv=read.table("Input_files/VTK2021_ASV_table.txt", header = TRUE, row.names = 1)[,1:20]
map=read.table("Input_files/VTK_metadata.txt", header = T, row.names = 1, sep = "\t")
map$Temperature2=paste(map$Temperature, "ºC", sep = "")
tax=read.table("Input_files/VTK2021_ASV_table.txt", header = TRUE, row.names = 1)[,22:27]
colnames(asv)=gsub("\\.", "-", colnames(asv))

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

#27 vs 29
C1=subset_samples(phy.all, Temperature %in% c("27", "29"))
res1=ancombc(phyloseq=C1,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_df=data.frame(Beta=res1[["res"]][["beta"]], Beta_se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])
colnames(res1_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_sig=subset(res1_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_sig$Diff_more_abundant=ifelse( res1_sig$W < 0 , "27", "29")
res1_sig$Taxa=paste(rownames(tax), " (",tax$Genus,";", tax$Family, ";",tax$Order, ") " , sep = "" )[match(rownames(res1_sig),rownames(tax))]
res1_sig$Comparison="Comparison 1: 27 vs 29"
message("Number of DA ASVs: ", nrow(res1_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res1_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 29: ", nrow(subset(res1_sig, Diff_more_abundant == "29" )))

#27 vs 32
C2=subset_samples(phy.all, Temperature %in% c("27", "32"))
res2=ancombc(phyloseq=C2,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res2_df=data.frame(Beta=res2[["res"]][["beta"]], Beta_se=res2[["res"]][["se"]], W=res2[["res"]][["W"]],pval=res2[["res"]][["p_val"]], qval=res2[["res"]][["q_val"]],DA=res2[["res"]][["diff_abn"]])
colnames(res2_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res2_sig=subset(res2_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res2_sig$Diff_more_abundant=ifelse( res2_sig$W < 0 , "27", "32")
res2_sig$Taxa=paste(rownames(tax), " (",tax$Genus,";", tax$Family, ";",tax$Order, ") " , sep = "" )[match(rownames(res2_sig),rownames(tax))]
res2_sig$Comparison="Comparison 2: 27 vs 32"
message("Number of DA ASVs: ", nrow(res2_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res2_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 32: ", nrow(subset(res2_sig, Diff_more_abundant == "32" )))

#27 vs 34
C3=subset_samples(phy.all, Temperature %in% c("27", "34"))
res3=ancombc(phyloseq=C3,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res3_df=data.frame(Beta=res3[["res"]][["beta"]], Beta_se=res3[["res"]][["se"]], W=res3[["res"]][["W"]],pval=res3[["res"]][["p_val"]], qval=res3[["res"]][["q_val"]],DA=res3[["res"]][["diff_abn"]])
colnames(res3_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res3_sig=subset(res3_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res3_sig$Diff_more_abundant=ifelse( res3_sig$W < 0 , "27", "34")
res3_sig$Taxa=paste(rownames(tax), " (",tax$Genus,";", tax$Family, ";",tax$Order, ") " , sep = "" )[match(rownames(res3_sig),rownames(tax))]
res3_sig$Comparison="Comparison 3: 27 vs 34"
message("Number of DA ASVs: ", nrow(res3_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res3_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 34: ", nrow(subset(res3_sig, Diff_more_abundant == "34" )))

#29 vs 32
C4=subset_samples(phy.all, Temperature %in% c("29", "32"))
res4=ancombc(phyloseq=C4,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res4_df=data.frame(Beta=res4[["res"]][["beta"]], Beta_se=res4[["res"]][["se"]], W=res4[["res"]][["W"]],pval=res4[["res"]][["p_val"]], qval=res4[["res"]][["q_val"]],DA=res4[["res"]][["diff_abn"]])
colnames(res4_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res4_sig=subset(res4_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res4_sig$Diff_more_abundant=ifelse( res4_sig$W < 0 , "29", "32")
res4_sig$Taxa=paste(rownames(tax), " (",tax$Genus,";", tax$Family, ";",tax$Order, ") " , sep = "" )[match(rownames(res4_sig),rownames(tax))]
res4_sig$Comparison="Comparison 4: 29 vs 32"
message("Number of DA ASVs: ", nrow(res4_sig), "\nNumber of DA ASVs enriched in 29: ", nrow(subset(res4_sig, Diff_more_abundant == "29" )), "\nNumber of DA ASVs enriched in 32: ", nrow(subset(res4_sig, Diff_more_abundant == "32" )))

#29 vs 34
C5=subset_samples(phy.all, Temperature %in% c("29", "34"))
res5=ancombc(phyloseq=C5,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res5_df=data.frame(Beta=res5[["res"]][["beta"]], Beta_se=res5[["res"]][["se"]], W=res5[["res"]][["W"]],pval=res5[["res"]][["p_val"]], qval=res5[["res"]][["q_val"]],DA=res5[["res"]][["diff_abn"]])
colnames(res5_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res5_sig=subset(res5_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res5_sig$Diff_more_abundant=ifelse( res5_sig$W < 0 , "29", "34")
res5_sig$Taxa=paste(rownames(tax), " (",tax$Genus,";", tax$Family, ";",tax$Order, ") " , sep = "" )[match(rownames(res5_sig),rownames(tax))]
res5_sig$Comparison="Comparison 5: 29 vs 34"
message("Number of DA ASVs: ", nrow(res5_sig), "\nNumber of DA ASVs enriched in 29: ", nrow(subset(res5_sig, Diff_more_abundant == "29" )), "\nNumber of DA ASVs enriched in 34: ", nrow(subset(res5_sig, Diff_more_abundant == "34" )))

#32 vs 34
C6=subset_samples(phy.all, Temperature %in% c("32", "34"))
res6=ancombc(phyloseq=C6,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res6_df=data.frame(Beta=res6[["res"]][["beta"]], Beta_se=res6[["res"]][["se"]], W=res6[["res"]][["W"]],pval=res6[["res"]][["p_val"]], qval=res6[["res"]][["q_val"]],DA=res6[["res"]][["diff_abn"]])
colnames(res6_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res6_sig=subset(res6_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res6_sig$Diff_more_abundant=ifelse( res6_sig$W < 0 , "32", "34")
res6_sig$Taxa=paste(rownames(tax), " (",tax$Genus,";", tax$Family, ";",tax$Order, ") " , sep = "" )[match(rownames(res6_sig),rownames(tax))]
res6_sig$Comparison="Comparison 6: 32 vs 34"
message("Number of DA ASVs: ", nrow(res6_sig), "\nNumber of DA ASVs enriched in 32: ", nrow(subset(res6_sig, Diff_more_abundant == "32" )), "\nNumber of DA ASVs enriched in 34: ", nrow(subset(res6_sig, Diff_more_abundant == "34" )))

## plots

ANCOMresults=rbind(res1_sig,res2_sig,res3_sig ,res4_sig,res5_sig,res6_sig)
#write.table(ANCOMresults,  "outputs/ANCOMBC_ASVs_results.txt", sep = "\t", quote = F, row.names = T )

ANCOMresults_plot=ANCOMresults %>% group_by(Diff_more_abundant, Comparison) %>% tally()
ANCOMresults_plot$n=ifelse(ANCOMresults_plot$Comparison == "Comparison 1: 27 vs 29" & ANCOMresults_plot$Diff_more_abundant == "27", ANCOMresults_plot$n*-1, ANCOMresults_plot$n*1)
ANCOMresults_plot$n=ifelse(ANCOMresults_plot$Comparison == "Comparison 2: 27 vs 32" & ANCOMresults_plot$Diff_more_abundant == "27", ANCOMresults_plot$n*-1, ANCOMresults_plot$n*1)
ANCOMresults_plot$n=ifelse(ANCOMresults_plot$Comparison == "Comparison 3: 27 vs 34" & ANCOMresults_plot$Diff_more_abundant == "27", ANCOMresults_plot$n*-1, ANCOMresults_plot$n*1)
ANCOMresults_plot$n=ifelse(ANCOMresults_plot$Comparison == "Comparison 4: 29 vs 32" & ANCOMresults_plot$Diff_more_abundant == "29", ANCOMresults_plot$n*-1, ANCOMresults_plot$n*1)
ANCOMresults_plot$n=ifelse(ANCOMresults_plot$Comparison == "Comparison 5: 29 vs 34" & ANCOMresults_plot$Diff_more_abundant == "29", ANCOMresults_plot$n*-1, ANCOMresults_plot$n*1)
ANCOMresults_plot$n=ifelse(ANCOMresults_plot$Comparison == "Comparison 6: 32 vs 34" & ANCOMresults_plot$Diff_more_abundant == "32", ANCOMresults_plot$n*-1, ANCOMresults_plot$n*1)

#pdf("./outputs/ANCOMBC_DA_barplots.pdf", width=6,height=4, pointsize = 12)
ggplot(data=ANCOMresults_plot, aes(x=Comparison, y=n)) + geom_bar(stat="identity", position = "dodge")  + geom_text(aes(label=n), vjust=0.5, color="white", position = position_dodge(1), size=3) +  theme_classic() + theme(axis.text.x=element_text(angle=90,hjust=1)) 
#dev.off()

#heat maps
ANCOMresults_subset=subset(ANCOMresults, abs(W) > 5)
asv_clr=apply(asv,2,clr)
asv_clr_subset=subset(asv_clr, rownames(asv_clr) %in% rownames(ANCOMresults_subset))
colnames(asv_clr_subset)=paste(map$Temperature, map$Genotype, sep = "_")[match(colnames(asv_clr_subset), rownames(map))]
asv_clr_sorted=asv_clr_subset[,c(1:5,16:20,6:10,11:15)]
labels=paste(rownames(tax), tax$Genus,tax$Family,tax$Order, sep = " | ")[match(rownames(asv_clr_sorted),rownames(tax))]
#pdf("outputs/ANCOM_heatmap.pdf", width = 7, height = 10, pointsize = 12)
pheatmap(asv_clr_sorted, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),  cellwidth = 7, labels_row = labels, cellheight =  5, fontsize_row = 5, fontsize_col= 8, legend = T, gaps_col = c(5,10,15),  cluster_rows = T,cluster_col = F, scale = "row", clustering_method = "centroid")#ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#dev.off()

