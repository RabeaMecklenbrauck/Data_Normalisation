#Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratObject)
library(cowplot)

#Read file
obj <- readRDS("data/4_Seurat_obj_CNA_referenceannotations.rds")
View(obj)

#Check QC
VlnPlot(obj, features = c("nFeature_gene", "nCount_gene"), ncol=2) #mtDNA high cells already excluded

plot1 <-FeatureScatter(obj, feature1 = "nCount_gene", feature2 = "nFeature_gene")
plot1

#Normalizing the data -> already done

#Identification of highly variable features
obj<-FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000) #change the number of variable features?
top10<-head(VariableFeatures(obj_new), 10)
plot1<-VariableFeaturePlot(obj_new)
plot2<-LabelPoints (plot = plot1, points = top10, repel = TRUE)
plot2

#Scale the data
all.genes<-rownames (obj)
obj<-ScaleData(obj, features = all.genes)

#Linear dimensional reduction
obj<-RunPCA(obj, features=VariableFeatures(object = obj))
#Visualize the PCA
print(obj[["pca"]], dims = 1:5, nfeatures =5)
VizDimLoadings(obj, dims = 1:5, reduction = 'pca') &
  theme (axis.text = element_text(size=5))
DimPlot(obj, reduction = 'pca')+NoLegend()
DimHeatmap(obj, dims = 1:20, cells = 500, balanced = TRUE)

#Dimensionality of the dataset
?ElbowPlot
ElbowPlot(obj,ndims=50)

#Cluster the cells
obj<-FindNeighbors(obj, dims = 1:40)
obj<-FindClusters(obj, resolution = 1)

#Run Umap
obj<-RunUMAP(obj, dims = 1:40)
DimPlot(obj, reduction = "umap")

#Annotate the Umap per patient
DimPlot(obj,label=T, reduction = "umap", group.by="predicted_CellType")

##############################################################################
##Create Umap for every patient
#Pt4
#Identification of highly variable features
pt4<-FindVariableFeatures(pt4, selection.method = "vst", nfeatures = 2000) #change the number of variable features?
top10<-head(VariableFeatures(pt4), 10)
plot1<-VariableFeaturePlot(pt4)
plot2<-LabelPoints (plot = plot1, points = top10, repel = TRUE)
plot2

#Scale the data
all.genes<-rownames (pt4)
pt4<-ScaleData(pt4, features = all.genes)

#Linear dimensional reduction
pt4<-RunPCA(pt4, features=VariableFeatures(object = pt4))
#Visualize the PCA
print(pt4[["pca"]], dims = 1:5, nfeatures =5)
VizDimLoadings(pt4, dims = 1:5, reduction = 'pca') &
  theme (axis.text = element_text(size=5))
DimPlot(pt4, reduction = 'pca')+NoLegend()
DimHeatmap(pt4, dims = 1:20, cells = 500, balanced = TRUE)

#Dimensionality of the dataset
?ElbowPlot
ElbowPlot(pt4,ndims=50)

#Cluster the cells
pt4<-FindNeighbors(pt4, dims = 1:25)
pt4<-FindClusters(pt4, resolution = 1)

#Run Umap
pt4<-RunUMAP(pt4, dims = 1:25)
DimPlot(pt4, reduction = "umap")

#Annotate the Umap per Timepoint
DimPlot(pt4,label=T, reduction = "umap", group.by="Sample")

#Pt9
#Identification of highly variable features
pt9<-FindVariableFeatures(pt9, selection.method = "vst", nfeatures = 2000) #change the number of variable features?
top10<-head(VariableFeatures(pt9), 10)
plot1<-VariableFeaturePlot(pt9)
plot2<-LabelPoints (plot = plot1, points = top10, repel = TRUE)
plot2

#Scale the data
all.genes<-rownames (pt9)
pt9<-ScaleData(pt9, features = all.genes)

#Linear dimensional reduction
pt9<-RunPCA(pt9, features=VariableFeatures(object = pt9))
#Visualize the PCA
print(pt9[["pca"]], dims = 1:5, nfeatures =5)
VizDimLoadings(pt9, dims = 1:5, reduction = 'pca') &
  theme (axis.text = element_text(size=5))
DimPlot(pt9, reduction = 'pca')+NoLegend()
DimHeatmap(pt9, dims = 1:20, cells = 500, balanced = TRUE)

#Dimensionality of the dataset
?ElbowPlot
ElbowPlot(pt9,ndims=50)

#Cluster the cells
pt9<-FindNeighbors(pt9, dims = 1:25)
pt9<-FindClusters(pt9, resolution = 1)

#Run Umap
pt9<-RunUMAP(pt9, dims = 1:25)
DimPlot(pt9, reduction = "umap")

#Annotate the Umap per Timepoint
DimPlot(pt9,label=T, reduction = "umap", group.by="Sample")
DimPlot(pt9, label = T, reduction = "umap", group.by = "clone.y")

#Pt11
#Identification of highly variable features
pt11<-FindVariableFeatures(pt11, selection.method = "vst", nfeatures = 2000) #change the number of variable features?
top10<-head(VariableFeatures(pt11), 10)
plot1<-VariableFeaturePlot(pt11)
plot2<-LabelPoints (plot = plot1, points = top10, repel = TRUE)
plot2

#Scale the data
all.genes<-rownames (pt11)
pt11<-ScaleData(pt11, features = all.genes)

#Linear dimensional reduction
pt11<-RunPCA(pt11, features=VariableFeatures(object = pt11))
#Visualize the PCA
print(pt11[["pca"]], dims = 1:5, nfeatures =5)
VizDimLoadings(pt11, dims = 1:5, reduction = 'pca') &
  theme (axis.text = element_text(size=5))
DimPlot(pt11, reduction = 'pca')+NoLegend()
DimHeatmap(pt11, dims = 1:20, cells = 500, balanced = TRUE)

#Dimensionality of the dataset
?ElbowPlot
ElbowPlot(pt11,ndims=50)

#Cluster the cells
pt11<-FindNeighbors(pt11, dims = 1:25)
pt11<-FindClusters(pt11, resolution = 1)

#Run Umap
pt11<-RunUMAP(pt11, dims = 1:25)
DimPlot(pt11, reduction = "umap")

#Annotate the Umap per Timepoint
DimPlot(pt11,label=T, reduction = "umap", group.by="Sample")
DimPlot(pt11, label=T, reduction = "umap", group.by="clone.y")


#Pt14
#Identification of highly variable features
pt14<-FindVariableFeatures(pt14, selection.method = "vst", nfeatures = 2000) #change the number of variable features?
top10<-head(VariableFeatures(pt14), 10)
plot1<-VariableFeaturePlot(pt14)
plot2<-LabelPoints (plot = plot1, points = top10, repel = TRUE)
plot2

#Scale the data
all.genes<-rownames (pt14)
pt14<-ScaleData(pt14, features = all.genes)

#Linear dimensional reduction
pt14<-RunPCA(pt14, features=VariableFeatures(object = pt14))
#Visualize the PCA
print(pt14[["pca"]], dims = 1:5, nfeatures =5)
VizDimLoadings(pt14, dims = 1:5, reduction = 'pca') &
  theme (axis.text = element_text(size=5))
DimPlot(pt14, reduction = 'pca')+NoLegend()
DimHeatmap(pt14, dims = 1:20, cells = 500, balanced = TRUE)

#Dimensionality of the dataset
ElbowPlot(pt14,ndims=50)

#Cluster the cells
pt14<-FindNeighbors(pt14, dims = 1:25)
pt14<-FindClusters(pt14, resolution = 1)

#Run Umap
pt14<-RunUMAP(pt14, dims = 1:25)
DimPlot(pt14, reduction = "umap")

#Annotate the Umap per Timepoint
DimPlot(pt14,label=T, reduction = "umap", group.by="Sample")

#Pt15
table(pt15$Sample)
