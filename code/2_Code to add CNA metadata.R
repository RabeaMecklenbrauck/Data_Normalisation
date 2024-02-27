#Load packages
library("Seurat")
library("SeuratObject")
library("readxl")
library("dplyr")
library("tidyverse")
library(tibble)

#Load Seurat-object
obj<-readRDS("data/2_lognorm_seurat.rds")

#Load the files with correct clone annotation
pt4_annotations<-read.csv("data/Excel files with final clone annotations/pt4_FINAL_CNA_SNV_clones.csv")
pt5_annotations<-read.csv("data/Excel files with final clone annotations/pt5_FINAL_CNA_SNV_clones.csv")
pt_9_annotations<-read.csv("data/Excel files with final clone annotations/pt9_FINAL_CNA_SNV_clones.csv")
pt10_annotations<-read.csv("data/Excel files with final clone annotations/pt10_clones_population.csv")
pt11_annotations<-read.csv("data/Excel files with final clone annotations/pt11_clones_population.csv")
pt14_annotations<-read.csv("data/Excel files with final clone annotations/pt14_FINAL_CNA_SNV_clones.csv")
pt15_annotations<-read.csv("data/Excel files with final clone annotations/pt15_FINAL_CNA_SNV_clones.csv")
pt18_annotations<-read.csv("data/Excel files with final clone annotations/pt18_clones_population.csv")
pt20_annotations<-read.csv("data/Excel files with final clone annotations/pt20_clones_population.csv")

#Merge the annotation df
#Make sure that in each DF the clone column is called the same
names(pt_9_annotations)[names(pt_9_annotations) == 'CNA_SNV_clone'] <- 'clone'
names(pt_9_annotations)[names(pt_9_annotations) == 'sampleID'] <- 'sampleID'
head(pt14_annotations)
names(pt14_annotations)[names(pt14_annotations) == 'CNA_SNV_clone'] <- 'clone'
names(pt14_annotations)[names(pt14_annotations) == 'sampleID'] <- 'sampleID'
head(pt14_annotations)
head(pt15_annotations)
names(pt15_annotations)[names(pt15_annotations) == 'CNA_SNV_clone'] <- 'clone'
head(pt15_annotations)
head(pt5_annotations)
names(pt5_annotations)[names(pt5_annotations) == 'CNA_SNV_clone'] <- 'clone'
names(pt5_annotations)[names(pt5_annotations) == 'Sample'] <- 'sampleID'
head(pt5_annotations)
head(pt4_annotations)
names(pt4_annotations)[names(pt4_annotations) == 'CNA_SNV_clone'] <- 'clone'
names(pt4_annotations)[names(pt4_annotations) == 'sampleID'] <- 'sampleID'
head(pt4_annotations)
head(pt20_annotations)
names(pt20_annotations)[names(pt20_annotations) == 'Sample.x'] <- 'sampleID'
names(pt20_annotations)[names(pt20_annotations) == 'id'] <- 'cell_id'
head(pt20_annotations)
head(pt18_annotations)
names(pt18_annotations)[names(pt18_annotations) == 'Sample'] <- 'sampleID'
names(pt18_annotations)[names(pt18_annotations) == 'id'] <- 'cell_id'
head(pt18_annotations)
head(pt11_annotations)
names(pt11_annotations)[names(pt11_annotations) == 'Sample'] <- 'sampleID'
head(pt11_annotations)
head(pt10_annotations)
names(pt10_annotations)[names(pt10_annotations) == 'Sample'] <- 'sampleID'
names(pt10_annotations)[names(pt10_annotations) == 'id'] <- 'cell_id'
pt4_clone<-pt4_annotations[ , c("cell_id", "clone", "Population", "sampleID")] 
pt5_clone<-pt5_annotations[ , c("cell_id", "clone", "Population", "sampleID")]
pt9_clone<-pt_9_annotations[ , c("cell_id", "clone", "Population", "sampleID")]
pt10_clone<-pt10_annotations[ , c("cell_id", "clone", "Population", "sampleID")]
pt11_clone<-pt11_annotations[ , c("cell_id", "clone", "Population", "sampleID")]
pt14_clone<-pt14_annotations[ , c("cell_id", "clone", "Population", "sampleID")]
pt15_clone<-pt15_annotations[ , c("cell_id", "clone", "Population", "sampleID")]
pt18_clone<-pt18_annotations[ , c("cell_id", "clone", "Population", "sampleID")]
pt20_clone<-pt20_annotations[ , c("cell_id", "clone", "Population", "sampleID")]

combined<-rbind(pt4_clone, pt5_clone, pt9_clone, pt10_clone, pt11_clone, pt14_clone, pt15_clone, pt18_clone, pt20_clone)
table(combined$sampleID, combined$clone)
write.csv( combined,file="results/new_clones_combined.csv")
combined_new_clones<-read.csv("results/new_clones_combined.csv")
rownames(combined)<-combined[,1]
combined_clones<-combined[,c('clone','sampleID')]
new_clones<-tibble::rownames_to_column(combined_clones, "cellID")
write.csv(new_clones, file="results/new_clones.csv")
write.csv(obj@meta.data, file="results/Metadata_transcriptome.csv")
rownames(combined_new_clones)<-combined_new_clones[,2]
combined_new_clones<-combined_new_clones[, c('clone')]
object<-obj@meta.data
object <- tibble::rownames_to_column(object, "cellID")

write.csv(object, "results/Seurat_metadata.csv")
combined_metadata<-left_join(object, new_clones, by="cellID")
write.csv(combined_metadata, "results/Metadata_combined.csv")
combined_metadata<-read.csv("results/Metadata_combined.csv")
sum(is.na(combined_metadata$clone.y))
new_clones_meta<-combined_metadata[c("clone.y", "cellID")]


#Add metadata with CNA annotations
combined_metadata<-read.csv("results/Metadata_combined.csv")
sum(is.na(combined_metadata$clone.y))
new_clones_meta<-combined_metadata[c("clone.y", "cellID")]
library(tibble)
new_clones_meta <- new_clones_meta %>% remove_rownames %>% column_to_rownames(var="cellID") 
obj<-AddMetaData(object=obj, new_clones_meta, col.name = 'clone.y')
metadata<-obj@meta.data



saveRDS(obj, "data/3_Seurat_CloneswithCNA.rds")
