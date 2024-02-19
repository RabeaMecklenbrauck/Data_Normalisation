install.packages("Seurat")
install.packages("SeuratObject")
library(Seurat)
library(SeuratObject)
library(tidyverse)

df<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240109_QC/sce_filtered.rds")

assay(df, "counts")  
assays(df)
colData(df) 
?as.Seurat
obj<-CreateSeuratObject(counts = counts(df), meta.data = as.data.frame(colData(df)))
obj@meta.data

#Now we have an object which is not normalized
obj_norm<-NormalizeData(obj)
DefaultAssay(obj_norm) <- "RNA"
obj_norm$nCount_RNA
#Add metadata with CNA annotations
combined_metadata<-read.csv("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Metadata_combined.csv")
sum(is.na(combined_metadata$clone.y))
new_clones_meta<-combined_metadata[c("clone.y", "cellID")]
install.packages("tibble")
library(tibble)
new_clones_meta <- new_clones_meta %>% remove_rownames %>% column_to_rownames(var="cellID") 


obj_norm<-AddMetaData(object=obj_norm, new_clones_meta, col.name = 'clone.y')

#Map this object to reference
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
#Load reference
ref<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240110_Mapping_to_reference/BoneMarrow_RefMap_SymphonyRef.rds")
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

#Optional: Visualize different existing annotations within the reference
p1 <- DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CyclePhase', raster=FALSE)
p2 <- FeaturePlot(ReferenceSeuratObj, reduction = 'umap', features = 'Pseudotime', raster=FALSE) 

p1 + p2
#Correct for a batch effect
#in this case I would correct for Patient
batchvar<-'Sample'

#Map the data, mapping done using Symphony (Kang et al, 2021)
mapped <- map_Query(
  exp_query = obj_norm@assays$RNA@layers$counts,  #at this point I can't directly access counts
  metadata_query = obj_norm@meta.data,
  ref_obj = ref,
  vars = batchvar
)
#Can't map to reference, no overlap in features???

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

#Mapping to reference failed
#Try Deseq2 with annotations from flow
#Create a subset for patient 11 BL and REL
pt11<-subset(x=obj_norm, subset=Patient=='pt11')
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
##Create a subset with only dominant population which is Precursors in this case
pt11_BL_REL<-subset(x=pt11_BL_REL, subset=Population=='Precursor')
##From this subset select the dominant clones
pt11_BL_REL <- subset(x=pt11_BL_REL, dominant == "1")
#Check whether the subsetting worked
table(pt11_BL_REL$Population)
table(pt11_BL_REL$dominant)

########DE####################
#Load libraries
library(ExperimentHub)
library(Seurat)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# a0counts aggregate to sample level
View(pt11_BL_REL@meta.data)
pt11_BL_REL$samples <- paste0(pt11_BL_REL$Sample, pt11_BL_REL$Sample_well)

DefaultAssay(pt11_BL_REL)

#Create dataset with counts and genes
cts <- AggregateExpression(pt11_BL_REL, 
                           group.by = c("samples"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

View(cts)
cts <- cts$RNA
View(cts)

# transpose
cts.t <- t(cts)

# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split cell type and sample
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame according row cell types which have been generated above
#cts.split <- split.data.frame(cts.t,
#f = factor(splitRows)) #only necessary when looking at a specific cell type

# fix colnames and transpose
cts.t.modified <- t(cts.t)
cts.t.modified[(1:10), (1:10)]


# Create Count matrix for DE
#Analyse for all cells
# 1. generate sample level metadata
colData_all <- data.frame(samples = colnames(cts.t.modified))
View(colData_all)

#3. generate a variable for the condition to compare (e.g. relapse and BL)
colData_all <- colData_all %>%
  mutate(condition = ifelse(grepl('BL', samples), 'baseline', 'relapse')) 



# perform DESeq2 
# 1.Create DESeq2 object   
dds_all <- DESeqDataSetFromMatrix(countData = cts.t.modified,
                                  colData = colData_all,
                                  design = ~ condition)

# 2.filter -> optional to exclude genes with a low number of reads
keep <- rowSums(counts(dds_all)) >=10
dds_all <- dds_all[keep,]

# 3. run DESeq2
dds_all <- DESeq(dds_all)

# 4. Check the coefficients for the comparison
resultsNames(dds_all)

# 5. Generate results object
res_all <- results(dds_all)
res_all
res_all<-as.data.frame(res_all)

#Visualize the results
#Create a table
res_tbl_all <-res_all%>% 
  data.frame()%>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_tbl_all
#Save the table
write.csv(res_tbl_all,"/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/pt11_BL_REL_Precursors_DE.csv" )

#Filter for only significant Genes
#Set threshold
padj_cutoff <- 0.05
#Subset for signifcant results 
sig_res_all<- subset(res_tbl_all, res_tbl_all$padj<0.05)
write.csv(sig_res_all, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/DataNormalisation/pt11_BL_REL_DE_Precursors_library_norm.csv" )

sig_res_mt <- sig_res_all[grepl('MT-', sig_res_all$gene),]

#Volcano Plot
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
#Create a vector where TRUE values denot padj ,0.05 and fold change >1.5
res_table_thres_all<- res_tbl_all %>% mutate(threshold=padj <0.05 & abs(log2FoldChange)>0.58)
EnhancedVolcano(res_table_thres_all,
                lab = res_table_thres_all$gene,
                x = 'log2FoldChange',
                y = 'pvalue')

ggplot(res_table_thres_all)+
  geom_point(aes(x=log2FoldChange, y= -log10(padj), color = threshold))+
  ggtitle("Volcano plot Pt11 BL vs REL")+
  xlab("log2 fold Change")+
  scale_y_continuous (limits = c(0,100))+
  lab = rownames(res_table_thres_all)+
  theme(legend.position = "none",
        plot.title = element_text (size =rel(1.5), hjust = 0.5),
        axis.title=element_text(size = rel(1.25)))
View(res_table_thres_all)
write.csv(res_table_thres_all, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Pt4_all_BL_REL_clones_DE_genes_VolcanoPlot.csv")