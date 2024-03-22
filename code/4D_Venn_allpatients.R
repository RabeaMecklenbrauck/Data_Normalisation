#Let's plot a Venn diagram of all shared genes
library(devtools)    
install_github("guokai8/VennDetail")
library(VennDetail)
library(readxl)
library(tidyverse)
library(utils)
library(readr)

#Load the included dataset
pt11<-read.csv("results/Pt11/DEA/Pt11_shared_genes.csv")
pt11_shared = subset(pt11, pt11$Subset == "Shared")
pt9<-read.csv("results/Pt9/Pt9_shared_genes.csv")
pt9_shared = subset (pt9, pt9$Subset=="Shared")
pt4<-read.csv("results/Pt4/DEA/Pt4_shared_genes.csv")
pt4_shared = subset (pt4, pt4$Subset=="Shared")
pt14<-read.csv("results/Pt14/DEA/Pt14_shared_genes.csv")
pt14_shared = subset (pt14, pt14$Subset=="Shared")


#Selecting terms of interest
l1<-earlyGMP[, c("gene", "log2FoldChange")]
l2<-LMPP[, c("gene", "log2FoldChange")]
l3<-MPP[, c("gene", "log2FoldChange")]
l4<-Progenitors[, c("gene", "log2FoldChange")]


#Make a very basic Venn diagramm
ven<-venndetail(list(pt11=pt11_shared$Detail, pt9=pt9_shared$Detail, pt4=pt4_shared$Detail, pt14=pt14_shared$Detail))

plot(ven)

#Get the gene list
detail(ven)
write_csv(result(ven), "results/allpt_shared_genes.csv")

#Now look at shared pathways
#KEGG
#Load the necessary files
Pt14_shared_KEGG<-read_csv("results/Pt14/Pt14_shared_KEGG_sig.csv")
Pt4_shared_KEGG<-read_csv("results/Pt4/Pt4_shared_KEGG_p<0.1.csv")
Pt9_shared_KEGG<-read_csv("results/Pt9/Pt9_shared_KEGG_p<0.1.csv")
Pt11_shared_KEGG<-read_csv("results/Pt11/fgsea/Pt11_shared_KEGG_sig.csv")

venn_KEGG<-venndetail(list(pt4 = Pt4_shared_KEGG$Detail, pt9 = Pt9_shared_KEGG$Detail, pt11 = Pt11_shared_KEGG$Detail, pt14 = Pt14_shared_KEGG$Detail))
plot(venn_KEGG)
write_csv(result(venn_KEGG), "results/ allpt_shared_pathways_KEGG.csv")

#GO
Pt14_shared_GO<-read_csv("results/Pt14/Pt14_shared_GO_sig.csv")
Pt4_shared_GO<-read_csv("results/Pt4/Pt4_shared_GO_sig.csv")
Pt9_shared_GO<-read_csv("results/Pt9/Pt9_shared_GO_sig.csv")
Pt11_shared_GO<-read_csv("results/Pt11/fgsea/Pt11_shared_GO_sig.csv")

venn_GO<-venndetail(list(pt4 = Pt4_shared_GO$Detail, pt9 = Pt9_shared_GO$Detail, pt11 = Pt11_shared_GO$Detail, pt14 = Pt14_shared_GO$Detail))
plot(venn_GO)
write_csv(result(venn_GO), "results/ allpt_shared_pathways_GO.csv")

#Hallmark
Pt14_shared_HALLMARK<-read_csv("results/Pt14/Pt14_shared_HALLMARK_sig.csv")
Pt4_shared_HALLMARK<-read_csv("results/Pt4/Pt4_shared_HALLMARK_sig.csv")
Pt9_shared_HALLMARK<-read_csv("results/Pt9/Pt9_shared_HALLMARK_sig.csv")
Pt11_shared_HALLMARK<-read_csv("results/Pt11/fgsea/Pt11_shared_HALLMARK_sig.csv")

venn_HALLMARK<-venndetail(list(pt4 = Pt4_shared_HALLMARK$Detail, pt9 = Pt9_shared_HALLMARK$Detail, pt11 = Pt11_shared_HALLMARK$Detail, pt14 = Pt14_shared_HALLMARK$Detail.1))
plot(venn_HALLMARK)
write_csv(result(venn_HALLMARK), "results/ allpt_shared_pathways_HALLMARK.csv")


