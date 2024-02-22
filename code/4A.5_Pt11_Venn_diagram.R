library(devtools)    
install_github("guokai8/VennDetail")
library(VennDetail)
library(readxl)
library(tidyverse)

#Load the included dataset
earlyGMP<-read_csv("results/Pt11/DEA/DEA_pt11_early GMP_padj.csv")
LMPP<-read_csv("results/Pt11/DEA/DEA_pt11_LMPP_padj.csv")
MPP<-read_csv("results/Pt11/DEA/DEA_pt11_MPP_padj.csv")
Progenitors<-read_csv("results/Pt11/DEA/DEA_pt11_Progenitors_padj.csv")

#Selecting terms of interest
l1<-earlyGMP[, c("gene", "log2FoldChange")]
l2<-LMPP[, c("gene", "log2FoldChange")]
l3<-MPP[, c("gene", "log2FoldChange")]
l4<-Progenitors[, c("gene", "log2FoldChange")]


#Make a very basic Venn diagramm
ven<-venndetail(list(GMP = l1$gene, LMPP = l2$gene, MPP = l3$gene, Progenitors = l4$gene))

plot(ven)

#Get the gene list
detail(ven)
write_csv(result(ven), "results/Pt11_shared_genes.csv")
