#Try to project the whole sample set to a reference map
#Load libraries
library(Seurat)
library(tidyverse)
library(symphony)
library(ggpubr)
library(patchwork)
#Load other packages not available on CRAN
if(!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AUCell", "doMC"))
if(!require(devtools, quietly = TRUE)) install.packages("devtools")
devtools::install_github("jaredhuling/jcolors")

#Load the BoneMarrowMap package
devtools::install_github('andygxzeng/BoneMarrowMap', force = TRUE)
library(BoneMarrowMap)

#Download the reference 

# Set directory to store projection reference files
projection_path = '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/'
#Extract metadata'

# Download Bone Marrow Reference - 344 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_SymphonyRef.rds', 
                    destfile = paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
# Download uwot model file - 221 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot', 
                    destfile = paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot'))

# Load Symphony reference
ref <- readRDS(paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
# Set uwot path for UMAP projection
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot')

#Visualize the Reference BM
ReferenceSeuratObj <- create_ReferenceObject(ref)
DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', 
        raster=FALSE, label=TRUE, label.size = 4)

#Optional: Visiualze different existing annotations within the reference
p1 <- DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CyclePhase', raster=FALSE)
p2 <- FeaturePlot(ReferenceSeuratObj, reduction = 'umap', features = 'Pseudotime', raster=FALSE) 

p1 + p2

#Load the Seurat object which you want to map to the reference
obj_new <- readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Seurat_CloneswithCNA.rds")

#Correct for a batch effect
#in this case I would correct for Patient
batchvar<-'Patient'

#Map the data, mapping done using Symphony (Kang et al, 2021)
mapped <- map_Query(
  exp_query = obj_new@assays$gene@counts, 
  metadata_query = obj_new@meta.data,
  ref_obj = ref,
  vars = batchvar
)

#Exclude cells with a high mapping error
mapped<-mapped%>% calculate_MappingError(., reference=ref, MAD_threshold = 2.5)

#Visualize QCs
QC_plots <- plot_MappingErrorQC(mapped)
patchwork::wrap_plots (QC_plots, ncol= 4, widths = c(0.8, 0.3, 0.8, 0.3))

#Remove cells with a high mapping error (optional)
mappedQC <- subset( mapped, mapping_error_QC == 'Pass')

#Cell Type assignment using 30 K-Nearest Neighbours
# Predict Hematopoietic Cell Types by KNN classification
mappedQC <- predict_CellTypes(
  query_obj = mappedQC, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

#Save mapped Seurat Object
saveRDS(mappedQC, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240129_Data_mapped_to reference")
#Show Umap
DimPlot(mappedQC, reduction = 'umap', group.by = c('predicted_CellType'), raster=FALSE, label=TRUE, label.size = 4)

#Predict Pseudotime values by k-nearest neighbours
mappedQC <- predict_Pseudotime(
  query_obj = mappedQC,
  ref_obj = ref,
  initial_label = 'initial Pseudotime',
  final_label = 'predicted_Pseudotime')
FeaturePlot(mappedQC, features = c ('predicted_Pseudotime'))

#Save an Excel file with the mapped annotations for each cell
save_ProjectionResults(
  query_obj = mappedQC,
  celltype_label = 'predicted_CellType',
  celltype_KNNprob_label = 'predicted_CellType_prob',
  pseudotime_label = 'predicted_Pseudotime',
  file_name = '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240129IvoVen_data_mapped_labelled.csv'
)

#Visualize projection density
#Set each condition to be visualized individually
batch_key <- 'Patient'
projection_plots <- plot_Projection_byDonor(
  query_obj = mappedQC, 
  batch_key = batch_key, 
  ref_obj = ref, 
  Hierarchy_only = FALSE, # Whether to exclude T/NK/Plasma/Stromal cells 
  downsample_reference = TRUE, 
  downsample_frac = 0.25,   # down-sample reference cells to 25%; reduces figure file size
  query_point_size = 0.2,   # adjust size of query cells based on # of cells
  saveplot = TRUE, 
  save_folder = '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis'
)

#Show plots together
patchwork::wrap_plots (projection_plots, ncol = 6)

#Generate a table with the number of each cell type in each sample
#with row per patient and cell type
query_composition <- get_Composition(
  query_obj = mappedQC, 
  donor_key = 'Patient', 
  celltype_label = 'predicted_CellType', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = NULL, 
  return_type = 'long')

query_composition
write.csv(query_composition, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240129_Ivo_Ven_all_cells_mapped.csv", row.names = FALSE)
#with one row per patient
query_composition2 <- get_Composition(
  query_obj = mappedQC, 
  donor_key = 'Patient', 
  celltype_label = 'predicted_CellType', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = NULL, 
  return_type = 'count')

query_composition2
write.csv(query_composition2, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240129_Ivo_Ven_all_cells_mapped_persample.csv", row.names = FALSE)
#with one row per patient and proportion of cells
query_composition3 <- get_Composition(
  query_obj = mappedQC, 
  donor_key = 'Patient', 
  celltype_label = 'predicted_CellType', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = NULL, 
  return_type = 'proportion')

query_composition3 
write.csv(query_composition2, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240129_Ivo_Ven_all_cells_mapped_percellproprotion.csv", row.names = FALSE)

table1<-table(mappedQC$Population, mappedQC$predicted_CellType, mappedQC$Sample)
write.csv(table1, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240129_Ivo_Ven_comparison_FACS_mapping.csv")


#Save the annotated Seurat object
saveRDS(mappedQC,"/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Seurat_obj_CNA_referenceannotations.rds" )




#Visualize projection density
#Set each condition to be visualized individually
batch_key <- 'Sample'
projection_plots <- plot_Projection_byDonor(
  query_obj = mappedQC, 
  batch_key = batch_key, 
  ref_obj = ref, 
  Hierarchy_only = FALSE, # Whether to exclude T/NK/Plasma/Stromal cells 
  downsample_reference = TRUE, 
  downsample_frac = 0.25,   # down-sample reference cells to 25%; reduces figure file size
  query_point_size = 0.2,   # adjust size of query cells based on # of cells
  saveplot = TRUE, 
  save_folder = '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis'
)

#Show plots together
patchwork::wrap_plots (projection_plots, ncol = 16)

#Generate a table with the number of each cell type in each sample
#with row per patient and cell type
query_composition_sample <- get_Composition(
  query_obj = mappedSampleQC, 
  donor_key = 'Sample', 
  celltype_label = 'predicted_CellType', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = NULL, 
  return_type = 'long')

query_composition_sample
write.csv(query_composition_sample, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Ivo_Ven_all_cells_mapped_sample.csv", row.names = FALSE)
#with one row per patient
query_composition_sample2 <- get_Composition(
  query_obj = mappedSampleQC, 
  donor_key = 'Sample', 
  celltype_label = 'predicted_CellType', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = NULL, 
  return_type = 'count')

query_composition_sample2 
write.csv(query_composition_sample2, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Ivo_Ven_all_cells_mapped_counts_sampleperrow.csv", row.names = FALSE)
#with one row per patient and proportion of cells
query_composition_sample3 <- get_Composition(
  query_obj = mappedSampleQC, 
  donor_key = 'Patient', 
  celltype_label = 'predicted_CellType', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = NULL, 
  return_type = 'proportion')

query_composition_sample3 #no cells from pt15 EOC3 and relapse are included?

#Check if those cells get eliminated by QC
mappedSample <- predict_CellTypes(
  query_obj = mappedSample, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

#Save mapped Seurat Object
saveRDS(mappedSampleQC, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240110_Data_mapped_to reference_batchcorrectionbysample")
#Show Umap
DimPlot(mappedSample, reduction = 'umap', group.by = c('predicted_CellType'), raster=FALSE, label=TRUE, label.size = 4)

#Predict Pseudotime values by k-nearest neighbours
batch_key='Sample'
mappedSample <- predict_Pseudotime(
  batch_key = batch_key,
  query_obj = mappedQC,
  ref_obj = ref,
  initial_label = 'initial Pseudotime',
  final_label = 'predicted_Pseudotime')
FeaturePlot(mappedSample, features = c ('predicted_Pseudotime'))

#Save an Excel file with the mapped annotations for each cell
save_ProjectionResults(
  query_obj = mappedSampleQC,
  celltype_label = 'predicted_CellType',
  celltype_KNNprob_label = 'predicted_CellType_prob',
  pseudotime_label = 'predicted_Pseudotime',
  file_name = '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/IvoVen_data_projected_labeled_batchcorrectionbysample.csv'
)

#Visualize projection density
#Set each condition to be visualized individually
batch_key <- 'Sample'
projection_plots <- plot_Projection_byDonor(
  query_obj = mapped, 
  batch_key = batch_key, 
  ref_obj = ref, 
  Hierarchy_only = FALSE, # Whether to exclude T/NK/Plasma/Stromal cells 
  downsample_reference = TRUE, 
  downsample_frac = 0.25,   # down-sample reference cells to 25%; reduces figure file size
  query_point_size = 0.2,   # adjust size of query cells based on # of cells
  saveplot = TRUE, 
  save_folder = '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis'
)

#Show plots together
patchwork::wrap_plots (projection_plots, ncol = 16)

#Project haematopoietic pseudotime
query <- predict_Pseudotime(
  query_obj = mappedQC, 
  ref_obj = ref, 
  initial_label = 'initial_Pseudotime',  # pseudotime assignments before filtering on mapping QC
  final_label = 'predicted_Pseudotime'   # pseudotime assignments with map QC failing cells assigned as NA
)

# Visualize Hematopoietic Pseudotime in query data
FeaturePlot(subset(mappedQC, mapping_error_QC == 'Pass'), features = c('predicted_Pseudotime'), split.by = 'Sample') & 
  scale_color_gradientn(colors = rev(brewer.pal(11, 'RdBu')))

# Simple heatmap to visualize composition of projected samples
p <- query_composition %>% 
  # show celltypes present in >1% of total cells
  select(Patient, colnames(query_composition)[-1][colSums(query_composition[-1]) > 0.01]) %>% 
  # convert to matrix and display heatmap
  column_to_rownames('Sample') %>% data.matrix() %>% ComplexHeatmap::Heatmap()
p

curl::curl_download('https://raw.githubusercontent.com/andygxzeng/BoneMarrowMap_Extras/main/AML_CellType_Genesets.gmt', 
                    destfile = paste0(projection_path, 'AML_CellType_Genesets.gmt'))
AMLgenesets <- load_Genesets_gmt('/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/AML_CellType_Genesets.gmt')
AMLgenesets %>% summary()
query <- score_Genesets_AUCell(mappedQC, genesets = AMLgenesets, nbatches = 1, ncores = 10, output = 'metadata')
query@meta.data
FeaturePlot(subset(query, mapping_error_QC == 'Pass'), 
            features = c('LSPC_Quiescent_AUC', 'LSC104_Ng2016_UP_AUC'), 
            split.by = 'Sample', max.cutoff = 'q99', min.cutoff = 'q5') & 
  scale_color_gradientn(colors = rev(brewer.pal(11, 'RdBu')))
