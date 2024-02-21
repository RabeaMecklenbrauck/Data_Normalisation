#Data Integration using Harmony
install.packages("Seurat")
library(Seurat)
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
devtools::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
if (!requireNamespace('remotes', quietly = TRUE) {
  install.packages('remotes')
}
remotes::install_github('satijalab/azimuth', ref = 'master')
library(Azimuth)
install.packages("ggplot2")
library(ggplot2)
install.packages("patchwork")
library(patchwork)
options(future.globals.maxSize = 1e9)

#Load annotated Seurat object
df<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240110_Mapping_to_reference/Seurat_obj_CNA_referenceannotations.rds")

#Split layers
df[["RNA"]]<-split(df[["RNA"]], f = df$Patient)

#Run Seurat workflow
df <- FindVariableFeatures(df)
df<-ScaleData(df)
df<-RunPCA(df)
df<-FindNeighbors(df, dims = 1:30, reduction = "pca")
df<-FindClusters(df, resolution = 2, cluster.name = "unintegrated")
df<-RunUMAP (df, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(df, reduction = "umap.unintegrated", group.by = c ("Patient", "predicted_CellType"))

#Perform HARMONY Integration
df<-IntegrateLayers(object = df, method = HarmonyIntegration, new.reduction = "harmony", verbose = FALSE)
#Cluster Cells after integration
df<-FindNeighbors(df, reduction = "harmony", dims = 1:30 )
df<-FindClusters(df, reduction = 2, cluster.name = "harmony")
df<-RunUMAP(df, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(df, reduction = "harmony", group.by = c("Patient", "predicted_CellType", "seurat_clusters"), combine = FALSE, label.size =2)

#Try again and split by Sample
obj<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240110_Mapping_to_reference/Seurat_obj_CNA_referenceannotations.rds")

#Split layers
obj[["RNA"]]<-split(obj[["RNA"]], f =obj$Sample)

#Run Seurat workflow
obj <- FindVariableFeatures(obj)
obj<-ScaleData(obj)
obj<-RunPCA(obj)
obj<-FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj<-FindClusters(obj, resolution = 2, cluster.name = "unintegrated")
obj<-RunUMAP (obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(obj, reduction = "umap.unintegrated", group.by = c ("Patient", "predicted_CellType"))
table(obj$predicted_CellType, obj$harmony_clusters)
#Perform HARMONY Integration
obj<-IntegrateLayers(object = obj, method = HarmonyIntegration, new.reduction = "harmony", verbose = FALSE)
#Cluster Cells after integration
obj<-FindNeighbors(obj, reduction = "harmony", dims = 1:10 )
obj<-FindClusters(obj, reduction = 2, cluster.name = "harmony_clusters")
obj<-RunUMAP(obj, reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony")
DimPlot_scCustom(obj, reduction = "umap.harmony", group.by = c("Patient", "harmony_clusters", "Sample", "clone.y"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
DimPlot(obj, reduction = "umap.harmony", group.by = c("Patient", "predicted_CellType", "harmony_clusters"), combine = FALSE, label.size =2)
saveRDS(obj, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240110_Mapping_to_reference/Seurat_obj_HARMONY.rds")
table(obj$predicted_CellType, obj$harmony_clusters)
#Run Harmony on a per patient basis
pt4<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/pt4_all_genes.rds")
#Split layers
pt4[["RNA"]]<-split(pt4[["RNA"]], f = pt4$Sample)

#Run Seurat workflow
pt4 <- FindVariableFeatures(pt4)
pt4<-ScaleData(pt4)
pt4<-RunPCA(pt4)
pt4<-FindNeighbors(pt4, dims = 1:30, reduction = "pca")
pt4<-FindClusters(pt4, resolution = 2, cluster.name = "unintegrated")
pt4<-RunUMAP (pt4, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(pt4, reduction = "umap.unintegrated", group.by = c ("Sample", "predicted_CellType"))
DimPlot_scCustom()
#Perform HARMONY Integration
pt4<-IntegrateLayers(object = pt4, method = HarmonyIntegration, new.reduction = "harmony", verbose = FALSE)
#Cluster Cells after integration
pt4<-FindNeighbors(pt4, reduction = "harmony", dims = 1:10 )
pt4<-FindClusters(pt4, reduction = 2, cluster.name = "harmony_clusters")
pt4<-RunUMAP(pt4, reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony")


DimPlot_scCustom(pt4, reduction = "umap.harmony", group.by = c("clone.y"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
table(pt4$clone.y, pt4$harmony_clusters)
table(pt4$predicted_CellType, pt4$harmony_clusters)
DimPlot(pt4, reduction = "umap.harmony", group.by = c("Sample", "predicted_CellType", "harmony_clusters"), combine = FALSE, label.size =2)
DimPlot_scCustom(pt4, reduction = "umap.harmony", group.by = c("clone.y", "harmony_clusters"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))

saveRDS(pt4, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/pt4_integrated.rds")
 #Change the colors of the clusters
umap4<-as.data.frame(pt4@meta.data)
coordinates<-pt4[["umap.harmony"]]@cell.embeddings
umap4<-merge(umap4, coordinates, by = "row.names", all = TRUE)
coordinates<-as.data.frame(coordinates)
clusters<-as.character(umap4$harmony_clusters)
gg_color_hue<-function(n){
  hues = seq (30,1000, length = n +1)
  hcl(h = hues, l = 65, c = 100) [1:n]
}
pal<-polychrome(alphabet)
pal<-gg_color_hue(length(unique(clusters)))
palopal<-Palo(coordinates, clusters, pal, color_blind_fun = 'deutan')
med<-aggregate(coordinates, list(clusters), mean)
install.packages("pals")
library(pals)
colnames(med) [1] <- "Cluster"              
ggplot() + geom_point(data = umap4, aes (x = umapharmony_1, y = umapharmony_2, col = clusters), size = 1, alpha = 1) +
  theme_void () +
  scale_color_manual(values = palopal)
table(pt4$predicted_CellType, pt4$harmony_clusters)
#pt9
#Run Harmony on a per patient basis
pt9<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/pt9_all_genes.rds")
#Split layers
pt9[["RNA"]]<-split(pt9[["RNA"]], f = pt9$Sample)

#Run Seurat workflow
pt9 <- FindVariableFeatures(pt9)
pt9<-ScaleData(pt9)
pt9<-RunPCA(pt9)
pt9<-FindNeighbors(pt9, dims = 1:30, reduction = "pca")
pt9<-FindClusters(pt9, resolution = 2, cluster.name = "unintegrated")
pt9<-RunUMAP (pt9, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(pt9, reduction = "umap.unintegrated", group.by = c ("Sample", "predicted_CellType"))

#Perform HARMONY Integration
pt9<-IntegrateLayers(object = pt9, method = HarmonyIntegration, new.reduction = "harmony", verbose = FALSE)
#Cluster Cells after integration
pt9<-FindNeighbors(pt9, reduction = "harmony", dims = 1:10 )
pt9<-FindClusters(pt9, reduction = 2, cluster.name = "harmony_clusters")
pt9<-RunUMAP(pt9, reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony")

DimPlot(pt9, reduction = "umap.harmony", group.by = c("Sample", "predicted_CellType", "harmony_clusters"), combine = FALSE, label.size =2)
DimPlot_scCustom(pt9, reduction = "umap.harmony", group.by = c("clone.y", "harmony_clusters"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
table(pt9$predicted_CellType, pt9$harmony_clusters)
#Change the colors of the clusters
umap9<-as.data.frame(pt9@meta.data)
coordinates<-pt9[["umap.harmony"]]@cell.embeddings
umap9<-merge(umap9, coordinates, by = "row.names", all = TRUE)
coordinates<-as.data.frame(coordinates)
clusters<-as.character(umap9$harmony_clusters)
gg_color_hue<-function(n){
  hues = seq (30,1000, length = n +1)
  hcl(h = hues, l = 65, c = 100) [1:n]
}
pal<-polychrome(alphabet)
pal<-gg_color_hue(length(unique(clusters)))
palopal<-Palo(coordinates, clusters, pal, color_blind_fun = 'deutan')
med<-aggregate(coordinates, list(clusters), mean)
colnames(med) [1] <- "Cluster"              
ggplot() + geom_point(data = umap9, aes (x = umapharmony_1, y = umapharmony_2, col = clusters), size = 1, alpha = 1) +
  theme_void () +
  scale_color_manual(values = palopal)
table(pt9$predicted_CellType, pt9$harmony_clusters)
saveRDS(pt9, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/pt9_integrated.rds")


#Pt11
pt11<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/pt11_all_genes.rds")
#Split layers
pt11[["RNA"]]<-split(pt11[["RNA"]], f = pt11$Sample)

#Run Seurat workflow
pt11 <- FindVariableFeatures(pt11)
pt11<-ScaleData(pt11)
pt11<-RunPCA(pt11)
pt11<-FindNeighbors(pt11, dims = 1:30, reduction = "pca")
pt11<-FindClusters(pt11, resolution = 2, cluster.name = "unintegrated")
pt11<-RunUMAP (pt11, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(pt11, reduction = "umap.unintegrated", group.by = c ("Sample", "predicted_CellType"))

#Perform HARMONY Integration
pt11<-IntegrateLayers(object = pt11, method = HarmonyIntegration, new.reduction = "harmony", verbose = FALSE)
#Cluster Cells after integration
pt11<-FindNeighbors(pt11, reduction = "harmony", dims = 1:10 )
pt11<-FindClusters(pt11, reduction = 2, cluster.name = "harmony_clusters")
pt11<-RunUMAP(pt11, reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony")
polychrome_pal<-DiscretePalette_scCustomize(num_colors =46, palette = "polychrome")
DimPlot_scCustom(pt11, reduction = "umap.harmony", group.by = c("Sample", "harmony_clusters"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
DimPlot_scCustom(pt11, reduction = "umap.harmony", group.by = c("clone.y", "harmony_clusters"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
table(pt11$predicted_CellType, pt11$harmony_clusters)
saveRDS(pt11, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/pt11_integrated.rds")
#Patient 14
pt14<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/pt14_all_genes.rds")
#Split layers
pt14[["RNA"]]<-split(pt14[["RNA"]], f = pt14$Sample)

#Run Seurat workflow
pt14 <- FindVariableFeatures(pt14)
pt14<-ScaleData(pt14)
pt14<-RunPCA(pt14)
pt14<-FindNeighbors(pt14, dims = 1:30, reduction = "pca")
pt14<-FindClusters(pt14, resolution = 2, cluster.name = "unintegrated")
pt14<-RunUMAP (pt14, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(pt14, reduction = "umap.unintegrated", group.by = c ("Sample", "predicted_CellType"))

#Perform HARMONY Integration
pt14<-IntegrateLayers(object = pt14, method = HarmonyIntegration, new.reduction = "harmony", verbose = FALSE)
#Cluster Cells after integration
pt14<-FindNeighbors(pt14, reduction = "harmony", dims = 1:10 )
pt14<-FindClusters(pt14, reduction = 2, cluster.name = "harmony_clusters")
pt14<-RunUMAP(pt14, reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony")
saveRDS(pt14,"/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/pt14_integrated.rds")
polychrome_pal<-DiscretePalette_scCustomize(num_colors =46, palette = "polychrome")
DimPlot_scCustom(pt14, reduction = "umap.harmony", group.by = c("Sample", "harmony_clusters"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
table(pt14$predicted_CellType, pt14$harmony_clusters)
#Control cells
ctrl<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Subsets for analysis/control_all_genes.rds")
#Split layers
ctrl[["RNA"]]<-split(ctrl[["RNA"]], f = ctrl$Plate)

#Run Seurat workflow
ctrl <- FindVariableFeatures(ctrl)
ctrl<-ScaleData(ctrl)
ctrl<-RunPCA(ctrl)
ctrl<-FindNeighbors(ctrl, dims = 1:30, reduction = "pca")
ctrl<-FindClusters(ctrl, resolution = 2, cluster.name = "unintegrated")
ctrl<-RunUMAP (ctrl, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ctrl, reduction = "umap.unintegrated", group.by = c ("predicted_CellType"))
DimPlot_scCustom(ctrl, reduction = "umap.unintegrated", group.by = c("predicted_CellType", "umap.unintegrated"), combine = FALSE, pt.size = 1, label.size =2, colors_use = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
#Perform HARMONY Integration
ctrl<-IntegrateLayers(object = ctrl, method = HarmonyIntegration, new.reduction = "harmony", verbose = FALSE)
#Cluster Cells after integration
ctrl<-FindNeighbors(ctrl, reduction = "harmony", dims = 1:10 )
ctrl<-FindClusters(ctrl, reduction = 2, cluster.name = "harmony_clusters")
ctrl<-RunUMAP(ctrl, reduction = "umap.harmony", dims = 1:10, reduction.name = "umap.harmony")

DimPlot(ctrl, reduction = "umap.harmony", group.by = c( "predicted_CellType", "harmony_clusters"), combine = FALSE, label.size =2)
Cluster_Highlight_Plot(obj, "control", highlight_color = "darkblue", background_color = "grey")
install.packages("SCpubr")
obj<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240110_Mapping_to_reference/Seurat_obj_HARMONY.rds")


SCpubr::do_DimPlot(obj, split.by = "Patient", integration = "umap.harmony", ncol = 6)
#Install different color palette
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("Winnie09/Palo")

library(Palo)
library(ggplot2)
obj1<-as.data.frame(obj@meta.data)
coordinates<-obj[["umap.harmony"]]@cell.embeddings
obj2<-merge(obj1, coordinates, by = "row.names", all = TRUE)
coordinates<-as.data.frame(coordinates)
clusters<-as.character(obj2$harmony_clusters)
gg_color_hue<-function(n){
  hues = seq (30,1000, length = n +1)
  hcl(h = hues, l = 65, c = 100) [1:n]
}
pal<-gg_color_hue(length(unique(clusters)))
palopal<-Palo(coordinates, clusters, pal)
med<-aggregate(coordinates, list(clusters), mean)
colnames(med) [1] <- "Cluster"              
ggplot() + geom_point(data = obj2, aes (x = umapharmony_1, y = umapharmony_2, col = clusters), size = 1, alpha = 1) +
  theme_void () +
  scale_color_manual(values = palopal)
table(obj1$predicted_CellType, obj1$harmony_clusters)  

