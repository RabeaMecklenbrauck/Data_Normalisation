####Script performs library normalisation on not normalised sce, creates a Seurat object, annotated from mapping and CNAs and dominant clone, does DEA with DeSeq2 and GSEA with fgsea
#Load libraries
#Reinstall Matrix to make sure it is at least version 1.6.5, otherwise it is not compatible with SeuratObject
install.packages("Matrix")
library(Matrix)
library(Seurat)
library(SeuratObject)
library(tidyverse)

#Load the not normalised sce
df<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/sce_filtered_not_normalised.rds")

#Create Seurat Object
obj<-CreateSeuratObject(counts = counts(df), meta.data = as.data.frame(colData(df)))

#Now we have an object which is not normalized, normalise the Seurat Object
obj_norm<-NormalizeData(obj)
DefaultAssay(obj_norm) <- "RNA"
obj_norm$nCount_RNA


#Add metadata with CNA annotations
combined_metadata<-read.csv("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Metadata_combined.csv")
new_clones_meta<-combined_metadata[c("clone.y", "cellID")]
library(tibble)
new_clones_meta <- new_clones_meta %>% remove_rownames %>% column_to_rownames(var="cellID") 
obj_norm<-AddMetaData(object=obj_norm, new_clones_meta, col.name = 'clone.y')

#Map this object to reference and create annotations of predicted celltype
#Add more packages
library(symphony)
library(ggpubr)
library(patchwork)
#Load the BoneMarrowMap package
library(BoneMarrowMap)
#Load reference
ref<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240110_Mapping_to_reference/BoneMarrow_RefMap_SymphonyRef.rds")
projection_path = '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/'

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

#Correct for a batch effect
#in this case I would correct for Sample
batchvar<-'Sample'

#Map the data, mapping done using Symphony (Kang et al, 2021)
mapped <- map_Query(
  exp_query = obj_norm@assays$RNA$counts,  #at this point I can't directly access counts
  metadata_query = obj_norm@meta.data,
  ref_obj = ref,
  vars = batchvar
)


#Exclude cells with a high mapping error
mapped<-mapped%>% calculate_MappingError(., reference=ref, MAD_threshold = 2.5)

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
saveRDS(mappedQC, '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/Seurat_object_lib_norm_annotated.rds')

#####Annotate dominant clone in patient
#Create a subset for patient 11 BL and REL
pt11<-subset(x=mappedQC, subset=Patient=='pt11')
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
pt11_BL_REL<-subset(x=pt11_BL_REL, subset=predicted_CellType=='LMPP')
##From this subset select the dominant clones
pt11_BL_REL <- subset(x=pt11_BL_REL, dominant == "1")
#Check whether the subsetting worked
table(pt11_BL_REL$predicted_CellType)
table(pt11_BL_REL$dominant)

##You should end up with 517 cells LMPP and dominant clone

########DEA using DeSeq2####################
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
pt11_BL_REL$samples <- paste0(pt11_BL_REL$Sample, pt11_BL_REL$Sample_well)

DefaultAssay(pt11_BL_REL)

#Create dataset with counts and genes
cts <- AggregateExpression(pt11_BL_REL, 
                           group.by = c("samples"),
                           assays = 'RNA',
                           slot = "counts", #here you get the warning that the slot argument is no longer used, but it works anyway
                           return.seurat = FALSE)

cts <- cts$RNA

# transpose
cts.t <- t(cts)

# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split cell type and sample
splitRows <- gsub('_.*', '', rownames(cts.t))

# fix colnames and transpose
cts.t.modified <- t(cts.t)
cts.t.modified[(1:10), (1:10)]


# Create Count matrix for DE
# 1. generate sample level metadata
colData_all <- data.frame(samples = colnames(cts.t.modified))

#2. generate a variable for the condition to compare (e.g. relapse and BL)
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
res_all <- results(dds_all, contrast = c("condition", "relapse", "baseline")) #the second park of the arguments specifies the direction of the comparison
res_all<-as.data.frame(res_all)

#Visualize the results
#Create a table
res_tbl_all <-res_all%>% 
  data.frame()%>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_tbl_all
#Save the table
write.csv(res_tbl_all, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/240220_DEA_pt11_LMPP_lib_norm.csv")

#Filter for only significant Genes
#Set threshold
padj_cutoff <- 0.05
#Subset for signifcant results 
sig_res_all<- subset(res_tbl_all, res_tbl_all$padj<0.05)
write.csv(sig_res_all, '/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/240220_DEA_pt11_LMPP_libnorm_padj.csv')

#Volcano Plot
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
#Create a vector where TRUE values denot padj ,0.05 and fold change >1.5
res_table_thres_all<- res_tbl_all %>% mutate(threshold=padj <0.05 & abs(log2FoldChange)>0.58)
EnhancedVolcano(res_table_thres_all,
                lab = res_table_thres_all$gene,
                x = 'log2FoldChange',
                y = 'pvalue')

#GSEA analysis
# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Set seed
set.seed(123456)
# Set project library
.libPaths("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation")
# Loading relevant libraries 
install.packages("tidyverse")
install.packages("RColorBrewer")
BiocManager::install("fgsea")
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)
# Set relevant paths
list.files()
bg_path <- "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/GSEA_references/Reference gene sets/msigdb_v2023.2.Hs_GMTs/"
out_path <- "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation"
in_path<-"/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation"

# Functions ===================================================
## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}
# Analysis ====================================================
## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/240220_DEA_pt11_LMPP_libnorm_padj.csv", row.names = 1)
## 2. Prepare background genes #leave this out if already done -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df$gene
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [13], GO = [53], reactome = [17]
bg_genes <- prepare_gmt(gmt_files[17], my_genes, savefile = FALSE)

#Prepare ranking of your DF gene
rankings<- sign(df$log2FoldChange)*(-log10(df$pvalue)) #use signed p values as preferred
names(rankings)<-df$gene
#Check list
head(rankings)
#Sort genes by signed p value
rankings<-sort(rankings, decreasing= TRUE)
#Plot
plot(rankings)
#Check for infinite values
max(rankings)
min(rankings)
#if there are infinite values adjust them for very high numbers e.g. 10x the highest non-infinite value
#max_ranking <- max(rankings[is.finite(rankings)])
#min_ranking <- min(rankings[is.finite(rankings)])
#rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
#rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
#rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

#Plot first ranked genes
ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Run GSEA
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres)
write_csv(GSEAres, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/GSEA_LMPP_pt11_libnorm_REACTOME.csv")

#Order by pvalue
head(GSEAres[order(pval),])
sum(GSEAres[,padj<0.05])
sum(GSEAres[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)

#Save
## 5. Save the results -----------------------------------------------
saveRDS(GSEAres, file = "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/240220GSEA_pt11_BL_REL_LMPP_Reactome.rds")
data.table::fwrite(GSEAres, file ="/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/240220_GSEA_pt11_BL_REL_LMPP_lib_norm_REACTOME.tsv" , sep = "\t", sep2 = c("", " ", ""))


######GO pathways
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [13], GO = [53], reactome = [16]
bg_genes <- prepare_gmt(gmt_files[53], my_genes, savefile = FALSE)

#Run GSEA
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres)
write_csv(GSEAres,"/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/240220_GSEA_pt11_LMPP_libnorm_GO.csv")

#Order by pvalue
head(GSEAres[order(pval),])
sum(GSEAres[,padj<0.05])
sum(GSEAres[,pval<0.01])


#topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
#topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
#topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)

# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)

#Save
name_of_comparison <- 'BL_REL'
background_genes <- 'GO'
filename <- paste0(out_path, 'GSEA/', name_of_comparison, '_', background_genes) 
saveRDS(GSEAres, file = "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/GSEA_pt11_BL_REL_LMPP_libnorm_GO.rds")
data.table::fwrite(GSEAres, file ="/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Data_Normalisation/240220_Comparison_normalisation/GSEA_pt11_BL_REL_LMPP_libnorm_GO.tsv" , sep = "\t", sep2 = c("", " ", ""))

#Visualize
install.packages("readr")
library(readr)
res_df<-readr::read_tsv("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/GSEA_pt11_BL_REL_GO.tsv")
BiocManager::install("enrichplot")
a
mutate(GSEAres, qscore=-log(padj, base =10)) %>% barplot(x="qscore")
