#Use FindMarkers on patient 11
df<-readRDS("data/4_Seurat_obj_CNA_referenceannotations.rds")
library(tidyverse)

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
pt11_BL_REL<-subset(x=pt11_BL_REL, subset=predicted_CellType=='Early GMP')
##From this subset select the dominant clones
pt11_BL_REL <- subset(x=pt11_BL_REL, dominant == "1")
#Check whether the subsetting worked
table(pt11_BL_REL$predicted_CellType)
table(pt11_BL_REL$dominant)

library(Seurat)
library(SeuratObject)
Idents(pt11_BL_REL)<-"Sample"
relapse.markers<-FindMarkers(pt11_BL_REL, ident.1 = "pt11-BL", ident.2 = "pt11-REL")
head(relapse.markers)

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

