# Felix A. Radtke
# 01 - creating of an SCE object, filtering, normalization, and generation of a filtered/normalized Seurat object

# Load libraries
library(tidyverse)
library(cowplot)
library(singleCellTK)
library(scuttle)
library(scater)

remotes::install_version("Matrix", version = "1.6-1.1") 



####START NEW SCE CREATION
transcriptome_indexes_merged <- readRDS("/data/0_transcriptome_indexes.rds")
matrix <- readRDS("data/0_matrix_reformatted.rds")

df <- data.frame(cell_id = colnames(matrix)) %>%
  mutate(n = row_number()) %>%
  left_join(transcriptome_indexes_merged[, c("cell_id", "Sort", "Plate", "Sample", "Patient", "Population", "clone", "OligodT_barcode_plate", "Sample_well")], by = "cell_id") %>%
  column_to_rownames("cell_id") # this generates the metadata that will be used for SCE object generation

# Make sce with reformatted matrix
sce <- SingleCellExperiment(assays = list(counts = matrix), colData = df)

# Feed ERCC features into AltExp within SCE (needs different processing)
is.ERCC <- grepl("ERCC-", rownames(sce))
sce <- splitAltExps(sce, ifelse(is.ERCC, "ERCC", "gene"))
altExpNames(sce)

if(sum(duplicated(rowData(sce)$feature_name)) > 0){
  stop("Problem with object assembly :(")
}

# Remove 
sce <- sce[, !is.na(sce$Sample_well)] # If only 192 barcodes are used then this will remove the unused barcodes

dim(sce) # Number of cells before filtering (after removal of unused barcodes) - since we use RAW input this should be all cells without any QC - in our case this is 14247

# Which features are mito
mito <- which(grepl("MT-", rownames(sce)))
# Add QC metrics to the SCE - this makes it easy to use them for inspection downstream, as well as regression
sce <- addPerCellQCMetrics(sce, subsets=list(Mt=mito))

metrics <- sce@colData %>% 
  as.data.frame() %>% 
  write_csv("results/IVO_VEN_REL_QC_metrics.csv")
??plotColData
plot <- plotColData(sce, x="Sample", y="detected") + ggtitle("Detected features") # way to plot directly for SCE object
plot + geom_hline(yintercept = 2000, linetype = "dashed", color = "red", size = 2)
print("results/IVO_VEN_REL_detected_per_sample.pdf", sep = ""), device = "pdf", width = 15)

# Summarize QC metrics

summary <- metrics %>% 
  group_by(Sample) %>% 
  summarize(Mean_genes_detected_per_cell = mean(detected),
            Mean_seq_depth_per_cell = mean(sum),
            Mean_mito_perc_per_cell = mean(subsets_Mt_percent),
            Mean_ERCC_perc_per_cell = mean(altexps_ERCC_percent)) 
%>% 
  write_csv("results/IVO_VEN_REL_QC_average_per_sample.csv")

summary <- metrics %>% 
  group_by(Patient) %>% 
  summarize(Mean_genes_detected_per_cell = mean(detected),
            Mean_seq_depth_per_cell = mean(sum),
            Mean_mito_perc_per_cell = mean(subsets_Mt_percent),
            Mean_ERCC_perc_per_cell = mean(altexps_ERCC_percent)) 
%>% 
  write_csv("results/IVO_VEN_REL_QC_average_per_patient.csv")
summary<-metrics %>% 
  group_by(Sample) %>%
  summarize (SD_genes_detected_per_cell = sd(detected),
             SD_seq_depth_per_cell = sd(sum),
             SD_mito_perc_per_cell = sd(subsets_Mt_percent),
             SD_ERCC_perc_per_cell = sd(altexps_ERCC_percent)) 
%>% 
  write_csv("results/IVO_VEN_REL_QC_sd_per_sample.csv")

# Plot QC metrics to get an idea

ggplot(metrics, aes(x=Sample, y=detected)) +
  geom_boxplot() +
  theme_cowplot()
ggsave(paste("results/IVO_VEN_REL_detected_per_sample_boxplots.pdf", sep = ""), device = "pdf", width = 15)

ggplot(metrics, aes(x=Patient, y=detected)) +
  geom_boxplot() +
  theme_cowplot()
ggsave(paste("results/IVO_VEN_REL_detected_per_patient_boxplots.pdf", sep = ""), device = "pdf", width = 10)

ggplot(metrics, aes(x=Sample, y=subsets_Mt_percent)) +
  geom_boxplot() +
  theme_cowplot() +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 2)

ggsave(paste("results/IVO_VEN_REL_MT_percent_per_sample_boxplots.pdf", sep = ""), device = "pdf", width = 15)

ggplot(metrics, aes(x=Sample, y = subsets_Mt_percent))+
  geom_violin()+
  theme_cowplot()+
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 2)

ggplot(metrics, aes(x=Sample, y=altexps_ERCC_percent)) +
  geom_boxplot() +
  theme_cowplot()

ggplot(metrics, aes(x=Sample, y=altexps_ERCC_percent)) +
  geom_violin() +
  theme_cowplot() +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", size = 2)
  
ggsave(paste("results/IVO_VEN_REL_ERCC_percent_per_sample_boxplots.pdf", sep = ""), device = "pdf", width = 15)

# Check if these are similar per batch - you might need to manually adjust thresholds per patient
qc.lib <- sce$sum < 25000  # Remove empty wells
qc.nexprs <- sce$detected < 2000 # Remove low complexity cells
qc.spike <- sce$altexps_ERCC_percent > 50 # Remove cells with > 50% ERCC (this is conservative, could potentially set much lower)
qc.mito <- sce$subsets_Mt_percent > 20 # Remove high mito cells
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito # Which cells are discarded

sce$discard <- discard # Add discard info to sce

# First get an idea how many cells would be discarded per sample if those metrics would be applied

test <- sce@colData %>% 
  as.data.frame() %>% 
  transform(discard = as.character(discard)) %>% 
  group_by(Sample,discard) %>% 
  summarise(n = n()) %>% 
  write_csv("results/IVO_VEN_REL_QC_pass_depth25000_det2000_ercc50_mt20.csv")

# How many cells are discarded, and why?
DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          SpikeProp=sum(qc.spike), MitoProp=sum(qc.mito), Total=sum(discard)) 



sce$discard <- discard # Add discard info to sce


plotColData(sce, x="Sample", y="sum", 
            colour_by="discard") + scale_y_log10() + ggtitle("Total count") + geom_hline(aes(yintercept = 25000), linetype = "longdash", color = "red") 
ggsave(paste("results/IVO_VEN_REL_depth_per_sample_discard.pdf", sep = ""), device = "pdf", width = 15)

plotColData(sce, x="Sample", y="detected", 
            colour_by="discard") + ggtitle("Detected features") + geom_hline(aes(yintercept = 2000), linetype = "longdash", color = "red")
ggsave(paste("results/IVO_VEN_REL_genes_per_sample_discard.pdf", sep = ""), device = "pdf", width = 15)

plotColData(sce, x="Sample", y="subsets_Mt_percent", 
            colour_by="discard") + ggtitle("Mito percent") + geom_hline(aes(yintercept = 15), linetype = "longdash", color = "red")
ggsave(paste("results/IVO_VEN_REL_mt_percent_per_sample_discard.pdf", sep = ""), device = "pdf", width = 15)

plotColData(sce, x="Sample", y="altexps_ERCC_percent", 
            colour_by="discard") + ggtitle("ERCC percent") + geom_hline(aes(yintercept = 50), linetype = "longdash", color = "red")
ggsave(paste("results/IVO_VEN_REL_ERCC_percent_per_sample_discard.pdf", sep = ""), device = "pdf", width = 15)

dim(sce)
sce_filtered <- sce[,discard == FALSE]
dim(sce_filtered)


saveRDS(sce_filtered, file = "/data/sce_filtered_not_normalised.rds") #this can be used if you want to do other form of normalisation
sce_filtered<-readRDS("/data/1_filtered.rds")
# SCRAN normalisation #
# Quick clustering
remove.packages("Matrix")
remove.packages("Seurat")
remotes::install_version("Matrix", version = "1.6-1.1") 
library(Matrix)
qclust<- scran::quickCluster(sce_filtered,)
# Size factors
sf.libsize <- librarySizeFactors(sce_filtered)
sf <- scran::calculateSumFactors(sce_filtered, clusters = qclust)
print("Size factors")
print(summary(sf))

plot(sf.libsize, sf, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(sce_filtered$Sample)))
abline(a=0, b=1, col="red")
install.packages("scuttle")
library(scuttle)
install.packages("Seurat")
sce_filtered <- logNormCounts(sce_filtered, size_factors=sf) # This is performing the log + normalization on the counts

sce_filtered <- swapAltExp(sce_filtered, "ERCC")

ercc.sf.libsize <- scuttle::librarySizeFactors(sce_filtered)

sce_filtered <- logNormCounts(sce_filtered, size_factors=ercc.sf.libsize)


# Add total ERCC count to metadata
sce_filtered$sf <- Matrix::colSums(counts(sce_filtered))

sce_filtered <- swapAltExp(sce_filtered, "gene")

library(Seurat)

obj <- as.Seurat(sce_filtered, counts = "counts") # converting the filtered and normalized sce to Seurat object

DefaultAssay(obj) <- "gene"

saveRDS(obj, file = "processed/2_lognorm_seurat.rds")

### Scatterplots

# Extract Seurat metadata for plotting in ggplot

metadata <- obj@meta.data

p1 <- ggplot(metadata, aes(nCount_gene, nFeature_gene, colour = Patient))+
  geom_point(alpha = 0.6, size = 0.5)+
  geom_vline(xintercept = 25000, color = "grey", linetype = "dashed")+
  geom_hline(yintercept = 2000, color = "grey", linetype = "dashed")+
  scale_x_log10()+ 
  labs(x= "Read count")+
  theme(legend.position = "none")+
  theme_cowplot()
p1
ggsave(paste("results/IVO_VEN_REL_depth_detected.pdf", sep = ""), device = "pdf", width = 9)

p2 <- ggplot(metadata, aes(altexps_ERCC_percent, subsets_Mt_percent, colour = Patient))+
  geom_point(size = 0.5)+
  geom_hline(yintercept = 15, color = "grey", linetype = "dashed")+
  geom_vline(xintercept = 50, color = "grey", linetype = "dashed")+
  labs()+
  theme_cowplot()
p2
ggsave(paste("results/IVO_VEN_REL_ERCC_MT.pdf", sep = ""), device = "pdf", width = 9)

p3 <- ggplot(metadata, aes(nFeature_gene, subsets_Mt_percent, colour = Patient))+
  geom_point(size = 0.5)+
  geom_hline(yintercept = 15, color = "grey", linetype = "dashed")+
  geom_vline(xintercept = 2000, color = "grey", linetype = "dashed")+
  labs()+
  theme_cowplot()
p3
ggsave(paste("results/IVO_VEN_REL_detected_MT.pdf", sep = ""), device = "pdf", width = 9)
