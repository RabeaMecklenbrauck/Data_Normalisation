#Venn diagramms of shared pathways
library(devtools)    
install_github("guokai8/VennDetail")
library(VennDetail)
library(readxl)
library(tidyverse)
library(utils)


#Shared genes Patient 9 and 11
pt9<-read_csv("results/Pt9/DEA_pt9_LMPP_MPP_padj.csv")
pt9_up= subset(pt9, pt9$log2FoldChange>0)
pt9_down=subset(pt9, pt9$log2FoldChange<0)
pt11<-read_csv("results/Pt11/DEA/DEA_pt11_LMPP_padj.csv")
pt11_up= subset(pt11, pt11$log2FoldChange>0)
pt11_down=subset(pt11, pt11$log2FoldChange<0)
#Select the terms of interest
pt9<-pt9[,c("gene", "log2FoldChange")]
pt11<-pt11[,c("gene", "log2FoldChange")]

venn_down<-venndetail(list(pt9_down=pt9_down$gene, pt11_down = pt11_down$gene))
plot(venn_down)
write_csv(result(venn_down),"results/Pt9and11_sharedgenes_down.csv")

venn_up<-venndetail(list(pt9_up=pt9_up$gene, pt11_up = pt11_up$gene))
plot(venn_up)
write_csv(result(venn_up),"results/Pt9and11_sharedgenes_up.csv")

#Look at shared pathways
Pt9_PW<-read_tsv("results/Pt9/GSEA_pt9_BL_REL_LMPP_MPP_all.tsv")
Pt11_PW<-read_tsv("results/GSEA_pt11_BL_REL_LMPP_all.tsv")
#Subset for significant genes
Pt9_PW_sig = subset(Pt9_PW, Pt9_PW$padj<0.05)
Pt11_PW_sig = subset(Pt11_PW, Pt11_PW$padj<0.05)

venn_PW<-venndetail(list(Pt9 = Pt9_PW_sig$pathway, Pt11=Pt11_PW_sig$pathway))
plot(venn_PW)
write_csv(result(venn_PW), "results/Pt9_11_shared_pathways.csv")
