#Use FindMarkers on patient 11
df<-readRDS("data/4_Seurat_obj_CNA_referenceannotations.rds")
library(tidyverse)
library(WriteXLS)

library(Seurat)
library(SeuratObject)

#Create a subset for patient 11 BL and REL
pt11<-subset(x=df, subset=Patient=='pt11')
#Extract metadata
data11<-pt11@meta.data
#Create variable marking the dominant clone at BL and REL
data11<-mutate(data11, dominant = if_else(
  Sample == 'pt11-BL'&clone.y %in% c("NDI", "NDIN")|
    Sample == 'pt11-REL' & clone.y == 'NDINF', "1","0"))
#Check whether it worked
table(data11$dominant, data11$clone.y, data11$Sample)
#Add metadata back to Seurat object
pt11<-AddMetaData(object = pt11, metadata = data11, col.name = 'dominant')
pt11_BL_REL<-subset(x=pt11, Sample %in% c("pt11-BL","pt11-REL"))

#Chose one population 
##Create a subset with only dominant population which is early GMPs in this case
pt11_BL_REL_GMP<-subset(x=pt11_BL_REL, subset=predicted_CellType=='Early GMP')
##From this subset select the dominant clones
pt11_BL_REL_GMP <- subset(x=pt11_BL_REL_GMP, dominant == "1")
#Check whether the subsetting worked
table(pt11_BL_REL_GMP$predicted_CellType)
table(pt11_BL_REL_GMP$dominant)
Idents(pt11_BL_REL_GMP)<-"Sample"
relapse.markers<-FindMarkers(pt11_BL_REL_GMP, ident.1 = "pt11-BL", ident.2 = "pt11-REL")
head(relapse.markers)
write.csv(relapse.markers, "results/Pt11_FindMarkers_GMP.csv")

##Create a subset with only dominant population which is LMPPs in this case
pt11_BL_REL_LMPP<-subset(x=pt11_BL_REL, subset=predicted_CellType=='LMPP')
##From this subset select the dominant clones
pt11_BL_REL_LMPP <- subset(x=pt11_BL_REL, dominant == "1")
#Check whether the subsetting worked
table(pt11_BL_REL$predicted_CellType)
table(pt11_BL_REL$dominant)
Idents(pt11_BL_REL_LMPP)<-"Sample"
relapse.markers.LMPP<-FindMarkers(pt11_BL_REL_LMPP, ident.1 = "pt11-BL", ident.2 = "pt11-REL")
head(relapse.markers.LMPP)
write.csv(relapse.markers.LMPP, "results/Pt11_FindMarkers_LMPP.csv")

#Create a subset with only dominant population which is MPPs in this case
pt11_BL_REL_MPP<-subset(x=pt11_BL_REL, subset=predicted_CellType %in% c("MPP-MkEry", "MPP-MyLy"))
##From this subset select the dominant clones
pt11_BL_REL_MPP <- subset(x=pt11_BL_REL, dominant == "1")

Idents(pt11_BL_REL_MPP)<-"Sample"
relapse.markers.MPP<-FindMarkers(pt11_BL_REL_MPP, ident.1 = "pt11-BL", ident.2 = "pt11-REL")
head(relapse.markers.MPP)
write.csv(relapse.markers.MPP, "results/Pt11_FindMarkers_MPP.csv")

#Create a subset with only dominant population which is Cycling Progenitors in this case
pt11_BL_REL_Progenitors<-subset(x=pt11_BL_REL, subset=predicted_CellType == "Cycling Progenitor")
##From this subset select the dominant clones
pt11_BL_REL_Progenitors <- subset(x=pt11_BL_REL_Progenitors, dominant == "1")

Idents(pt11_BL_REL_Progenitors)<-"Sample"
relapse.markers.Progenitors<-FindMarkers(pt11_BL_REL_Progenitors, ident.1 = "pt11-BL", ident.2 = "pt11-REL")
head(relapse.markers.Progenitors)
write.csv(relapse.markers.Progenitors, "results/Pt11_FindMarkers_Progenitors.csv")

#Let's plot a UMAP of the whole patient 11 and try to locate the DE genes
pt11<-FindVariableFeatures(pt11, selection.method = "vst",nfeatures =2000)
all.genes <- rownames(pt11)
pt11<-ScaleData(pt11, features = all.genes)
pt11<-RunPCA(pt11, features = VariableFeatures(object = pt11))
ElbowPlot(pt11)
pt11<-FindNeighbors(pt11, dims = 1:20)
pt11<-FindClusters(pt11, resolution = 0.5)
DimPlot(pt11, reduction = "umap") #The UMAP now looks like the reference map??? Why is that? I'm just using the annotations...
#I'll just check this with the whole df
df<-FindVariableFeatures(df, selection.method = "vst",nfeatures =2000)
all.genes <- rownames(df)
df<-ScaleData(df, features = all.genes)
df<-RunPCA(df, features = VariableFeatures(object = df))
ElbowPlot(df)
df<-FindNeighbors(df, dims = 1:20)
df<-FindClusters(df, resolution = 0.5)
DimPlot(df, reduction = "umaps")
#This also maps to the reference now. Why is that?
FeaturePlot(pt11, features = c("PDE4D", "MSRB3", "INPP4B"))
#Lets try to find the markers in the UMAP
library(scCustomize)
DimPlot_scCustom(pt11, reduction = "umap", group.by = c("Sample"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
DimPlot_scCustom(pt11, reduction = "umap", group.by = c("clone.y"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
#These Plots show quite nicely that at REL the cells map much more to LMPP/MLP and cycling progenitors when at BL they map mostly to early GMPs
#So what happens if we take the dominant clones only
DimPlot_scCustom(pt11_BL_REL, reduction = "umap", group.by = c("Sample"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
DimPlot_scCustom(pt11_BL_REL, reduction = "umap", group.by = c("clone.y"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
#The cell numbers are too low to mak nice UMAPS...