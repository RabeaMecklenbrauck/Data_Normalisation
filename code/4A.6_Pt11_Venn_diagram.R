library(devtools)    
install_github("guokai8/VennDetail")
library(VennDetail)
library(readxl)
library(tidyverse)
library(utils)

#Load the included dataset
earlyGMP<-read_csv("results/pt11/DEA/DEA_pt11_early GMP_padj.csv")
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

#Clean environment
#Compare that to the marker genes found in Find Markers
library(readxl)
library(WriteXLS)
GMP<-read_xlsx("results/Pt11/DEA/Pt11_FindMarkers_GMP_padj.xlsx")
LMPP<-read_xls("results/Pt11/DEA/Pt11_FindMarkers_LMPP_padj.xls")
MPP<-read_xlsx("results/Pt11/DEA/Pt11_FindMarkers_MPP_padj.xlsx")
Progenitors<-read_xlsx("results/Pt11/DEA/Pt11_FindMarkers_Progenitors_padj.xlsx")

#Selecting terms of interest
l1<-GMP[, c("...1", "avg_log2FC")]
l2<-LMPP[, c("...1", "avg_log2FC")]
l3<-MPP[, c("...1", "avg_log2FC")]
l4<-Progenitors[, c("...1", "avg_log2FC")]

ven<-venndetail(list(GMP = l1$...1, LMPP = l2$...1, MPP = l3$...1, Progenitors = l4$...1))
plot(ven)

#Get the gene list
detail(ven)
WriteXLS(result(ven), "results/Pt11_shared_genes_FindMarkers.csv")
