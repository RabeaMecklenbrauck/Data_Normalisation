library(devtools)    
install_github("guokai8/VennDetail")
library(VennDetail)
library(readxl)
library(tidyverse)
library(utils)

#Load the included dataset
HSC<-read_csv("results/Pt4/DEA/DEA_pt4_HSC_padj.csv")
LMPP<-read_csv("results/Pt4/DEA/DEA_pt4_LMPP_padj.csv")
MPP<-read_csv("results/Pt4/DEA/DEA_pt4_MPP_padj.csv")
MEP<-read_csv("results/Pt4/DEA/DEA_pt4_MEP_padj.csv")

#Selecting terms of interest
l1<-HSC[, c("gene", "log2FoldChange")]
l2<-LMPP[, c("gene", "log2FoldChange")]
l3<-MPP[, c("gene", "log2FoldChange")]
l4<-MEP[, c("gene", "log2FoldChange")]


#Make a very basic Venn diagramm
ven<-venndetail(list(HSC = l1$gene, LMPP = l2$gene, MPP = l3$gene, MEP = l4$gene))

plot(ven)

#Get the gene list
detail(ven)
write_csv(result(ven), "results/Pt4/DEA/Pt4_shared_genes.csv")
