library(devtools)    
install_github("guokai8/VennDetail")
library(VennDetail)
library(readxl)
library(tidyverse)
library(utils)

#Load the included dataset
LMPP<-read_csv("results/Pt14/DEA/DEA_pt14_LMPP_padj.csv")
MPP<-read_csv("results/Pt14/DEA/DEA_pt14_MPP_padj.csv")

#Selecting terms of interest
l1<-LMPP[, c("gene", "log2FoldChange")]
l2<-MPP[, c("gene", "log2FoldChange")]


#Make a very basic Venn diagramm
ven<-venndetail(list(LMPP = l1$gene, MPP = l2$gene))

plot(ven)

#Get the gene list
detail(ven)
write_csv(result(ven), "results/Pt14/DEA/Pt14_shared_genes.csv")
