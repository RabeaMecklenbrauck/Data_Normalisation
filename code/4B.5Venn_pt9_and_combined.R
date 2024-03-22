library(devtools)    
install_github("guokai8/VennDetail")
library(VennDetail)
library(readxl)
library(tidyverse)
library(utils)

#Load the included dataset
MPP<-read_csv("results/DEA_pt9_MPP_padj.csv")
LMPP<-read_csv("results/DEA_pt9_LMPP_padj.csv")
HSC<-read_csv("results/DEA_pt9_HSC_padj.csv")


#Selecting terms of interest
l1<-HSC[, c("gene", "log2FoldChange")]
l2<-LMPP[, c("gene", "log2FoldChange")]
l3<-MPP[, c("gene", "log2FoldChange")]



#Make a very basic Venn diagramm
ven<-venndetail(list(HSC= l1$gene, LMPP = l2$gene, MPP = l3$gene))

plot(ven)

#Get the gene list
detail(ven)
write_csv(result(ven), "results/Pt9_shared_genes.csv")

#Bonus: Create a Ven diagramm with all FLT3 mutant patients
#add gebne lists for pt11
#Load the included dataset
MPP11<-read_csv("results/DEA_pt11_MPP_padj.csv")
LMPP11<-read_csv("results/DEA_pt11_LMPP_padj.csv")
earlyGMP<-read_csv("results/DEA_pt11_early GMP_padj.csv")
Progenitors<-read_csv("results/DEA_pt11_Progenitors_padj.csv")

#Selecting terms of interest
l4<-MPP11[, c("gene", "log2FoldChange")]
l5<-LMPP11[, c("gene", "log2FoldChange")]
l6<-earlyGMP[, c("gene", "log2FoldChange")]
l7<-Progenitors[, c("gene", "log2FoldChange")]


venn_comb<-venndetail(list(HSC9= l1$gene, LMPP9 = l2$gene, MPP9 = l3$gene, MPP11 = l4$gene, LMPP11 = l5$gene, earlyGMP11 = l6$gene, progenitors = l7$gene))
plot(venn_comb)


if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library(ggVennDiagram)                         
ggVennDiagram(list_all, label_alpha=0, set_color = c("blue", "red", "yellow", "green", "orange", "purple", "pink"))+scale_fill_gradient(low = "white", high = "lightblue")

list_all<-list(HSC9=l1$gene,
               LMPP9=l2$gene,
               MPP9=l3$gene,
               MPP11=l4$gene,
               LMPP11=l5$gene,
               GMP11=l6$gene,
               Progenitors11=l7$gene
)
ven<-ggVennDiagram(list_all, label_alpha=0)
ven_result<-result(ven)
v=Venn(list_all)
intersect(l1$gene, l2$gene, l3$gene, l4$gene, l5$gene, l6$gene, l7$gene)
d=process_data(v)
list_shared_genes<-venn_region(d)|>as.data.frame()
WriteXLS::WriteXLS(list_shared_genes, "results/List_shared_genes_pt11_pt9.xls")
#So that looks like a pretty flower, not very useful
#Let's try something else
library(UpSetR)
library(ComplexHeatmap)
list_to_matrix(list_all)
m1= make_comb_mat(list_all) #uses default mode "distinct"
m1
m2=make_comb_mat(list_all, mode="intersect")
m3=make_comb_mat(list_all, mode="union")

UpSet(m1)
UpSet(m2)
UpSet(m3)


