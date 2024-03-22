#This script uses the normalised data and performs DEA (DeSeq2) and fgsea
#Open the annotated Seurat object
df<-readRDS("data/4_Seurat_obj_CNA_referenceannotations.rds")
library(tidyverse)
library(Seurat)


#Create a subset for patient 11 BL and REL
pt4<-subset(x=df, subset=Patient=='pt4')
#Extract metadata
data4<-pt4@meta.data
#Create variable marking the dominant clone at BL and REL
#I chose SIAJ at BL and SIA_20_1, SIA_20_7 and SIAJ_1 at REL
data4<-mutate(data4, dominant = if_else(
  pt_status == 'pt4 BL'&clone.y %in% c("SIAJ")|
    pt_status == 'pt4 REL' & clone.y %in% c("SIA_20_1", "SIA_20_7", "SIAJ_1"), "1","0"))
#Check whether it worked
table(data4$dominant, data4$clone.y, data4$Sample)
#Add metadata back to Seurat object
pt4<-AddMetaData(object = pt4, metadata = data4, col.name = 'dominant')
pt4_BL_REL<-subset(x=pt4, pt_status %in% c("pt4 BL","pt4 REL"))

#According to this table the major clone at BL and at REL is mainly found in the continuum of HSC, MPP, LMPP and MEP
#We'lldo MPP and combine MPP-MkEry and MPP-MyLy
##Create a subset with only dominant population which is LMPP in this case
pt4_BL_REL<-subset(x=pt4_BL_REL, subset=predicted_CellType %in% c("MPP-MkEry","MPP-MyLy"))
##From this subset select the dominant clones
pt4_BL_REL <- subset(x=pt4_BL_REL, dominant == "1")
#Check whether the subsetting worked
table(pt4_BL_REL$predicted_CellType)
table(pt4_BL_REL$dominant)

#This leaves you with 418 MPP


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
View(pt4_BL_REL@meta.data)
pt4_BL_REL$samples <- paste0(pt4_BL_REL$Sample, pt4_BL_REL$Sample_well)

DefaultAssay(pt4_BL_REL)

#Create dataset with counts and genes
cts <- AggregateExpression(pt4_BL_REL, 
                           group.by = c("samples"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA
#Ignore warnings about Slot and layers
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
#Analyse for all cells
# 1. generate sample level metadata
colData_all <- data.frame(samples = colnames(cts.t.modified))
View(colData_all)

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
res_all <- results(dds_all, contrast = c("condition", 'relapse', 'baseline'))
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
write.csv(res_tbl_all,"results/Pt4/DEA/DEA_pt4_MPP.csv" )

#Filter for only significant Genes
#Set threshold
padj_cutoff <- 0.05
#Subset for signifcant results 
sig_res_all<- subset(res_tbl_all, res_tbl_all$padj<0.05)
write.csv(sig_res_all, "results/Pt4/DEA/DEA_pt4_MPP_padj.csv")

#Explore the data
#Normalise the data, don't use rlog, that does take ages
vst<-vst(dds_all)
assay(vst) [1:5, 1:5]
#Compare the normalised and not normalised data
plot(log2(1+counts(dds_all, normalized = TRUE) [,1:2]), col = "black", pch = 20, ces = 0.3)
plot(assay(vst)[,1:2], col = "black", pch = 20, cex = 0.3)

sampleDist<-dist(t(assay(vst)))
as.matrix(sampleDist)[1:5, 1:5]

sampleDistMatrix <- as.matrix( sampleDist )
rownames(sampleDistMatrix) <- paste( vst$condition, 
                                     vst$samples, sep="-" )
colnames(sampleDistMatrix) <- NULL   
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

#Plot the distance between conditions
plotPCA(vst, intgroup = c("condition"))

library(genefilter)

topVarGenes<-head(order(rowVars(assay(vst)), decreasing = TRUE), 100)
heatmap.2( assay(vst)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c(baseline="darkgreen", relapse = "darkred")[
             colData(vst)$condition ] )

#Volcano Plot
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
#Create a vector where TRUE values denote padj ,0.05 and fold change >1.5
res_table_thres_all<- res_tbl_all %>% mutate(threshold=padj <0.05 & abs(log2FoldChange)>0.58)
EnhancedVolcano(res_table_thres_all,
                lab = res_table_thres_all$gene,
                x = 'log2FoldChange',
                y = 'pvalue')#Explore the data
#Normalise the data, don't use rlog, that does take ages
vst<-vst(dds_all)
assay(vst) [1:5, 1:5]
#Compare the normalised and not normalised data
plot(log2(1+counts(dds_all, normalized = TRUE) [,1:2]), col = "black", pch = 20, ces = 0.3)
plot(assay(vst)[,1:2], col = "black", pch = 20, cex = 0.3)

sampleDist<-dist(t(assay(vst)))
as.matrix(sampleDist)[1:5, 1:5]

sampleDistMatrix <- as.matrix( sampleDist )
rownames(sampleDistMatrix) <- paste( vst$condition, 
                                     vst$samples, sep="-" )
colnames(sampleDistMatrix) <- NULL   
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

#Plot the distance between conditions
plotPCA(vst, intgroup = c("condition"))

library(genefilter)

topVarGenes<-head(order(rowVars(assay(vst)), decreasing = TRUE), 100)
heatmap.2( assay(vst)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c(baseline="darkgreen", relapse = "darkred")[
             colData(vst)$condition ] )

#Volcano Plot
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
#Create a vector where TRUE values denote padj ,0.05 and fold change >1.5
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
.libPaths("results/")
# Loading relevant libraries 
install.packages("tidyverse")
install.packages("RColorBrewer")
BiocManager::install("fgsea")
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)
# Set relevant paths
list.files()
bg_path <- "data/GSEA_references/Reference gene sets/msigdb_v2023.2.Hs_GMTs/"
out_path <- "results/Pt4/fgsea"
in_path<-"results/"

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
df_new <- read.csv("results/Pt4/DEA/DEA_pt4_MPP_padj.csv", row.names = 1)
## 2. Prepare background genes #leave this out if already done -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- df_new$gene
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
#GO pathways
bg_genes_GO <- prepare_gmt(gmt_files[53], my_genes, savefile = FALSE)

#Prepare ranking of your DF gene
rankings<- sign(df_new$log2FoldChange)*(-log10(df_new$padj)) #use signed p values as preferred
names(rankings)<-df_new$gene
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

#Plot first ranked genes of you want (optional)
ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Run GSEA
GSEAres_GO <- fgsea(pathways = bg_genes_GO, # List of gene sets to check
                    stats = rankings,
                    scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                    minSize = 10,
                    maxSize = 500,
                    nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres_GO)
write_csv(GSEAres_GO, "results/Pt4/DEA/GSEA_MPP_pt4_GO.csv")

#Order by pvalue
head(GSEAres_GO[order(pval),])
sum(GSEAres_GO[,padj<0.05])
sum(GSEAres_GO[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes_GO[[head(GSEAres_GO[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_GO[order(padj), ], 1)$pathway)

#Save
## 5. Save the results -----------------------------------------------
saveRDS(GSEAres_GO, file = "results/Pt4/fgsea/GSEA_pt4_BL_REL_MPP_GO.rds")
data.table::fwrite(GSEAres_GO, file ="results/Pt4/fgsea/GSEA_pt9_BL_REL_MPP_GO.tsv" , sep = "\t", sep2 = c("", " ", ""))

##6.Visualise
GSEAres_GO$adjPvalue <- ifelse(GSEAres_GO$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_GO_sig=subset(GSEAres_GO, GSEAres_GO$padj<=0.1)
ggplot(GSEAres_GO_sig, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO pathways Enrichment Score from GSEA (p<0.1)")


######KEGG pathways
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes_KEGG <- prepare_gmt(gmt_files[11], my_genes, savefile = FALSE)

#Run GSEA
GSEAres_KEGG <- fgsea(pathways = bg_genes_KEGG, # List of gene sets to check
                      stats = rankings,
                      scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                      minSize = 10,
                      maxSize = 500,
                      nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres_KEGG)
write_csv(GSEAres_KEGG,"results/GSEA_pt4_MPP_KEGG.csv")

#Order by pvalue
head(GSEAres_KEGG[order(pval),])
sum(GSEAres_KEGG[,padj<0.05])
sum(GSEAres_KEGG[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes_KEGG[[head(GSEAres_KEGG[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_KEGG[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres_KEGG, file = "results/GSEA_pt4_BL_REL_MPP_KEGG.rds")
data.table::fwrite(GSEAres_KEGG, file ="results/GSEA_pt4_BL_REL_MPP_KEGG.tsv" , sep = "\t", sep2 = c("", " ", ""))

#Visualise
GSEAres_KEGG$adjPvalue <- ifelse(GSEAres_KEGG$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_KEGG_sig=subset(GSEAres_KEGG, GSEAres_KEGG$padj<=0.1)
ggplot(GSEAres_KEGG , aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways Enrichment Score from GSEA")
######Hallmark pathways
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes_HALLMARK <- prepare_gmt(gmt_files[67], my_genes, savefile = FALSE)

#Run GSEA
GSEAres_HALLMARK <- fgsea(pathways = bg_genes_HALLMARK, # List of gene sets to check
                          stats = rankings,
                          scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                          minSize = 10,
                          maxSize = 500,
                          nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres_HALLMARK)
write_csv(GSEAres_HALLMARK,"results/Pt4/fgsea/GSEA_pt4_MPP_Hallmark.csv")

#Order by pvalue
head(GSEAres_HALLMARK[order(pval),])
sum(GSEAres_HALLMARK[,padj<0.05])
sum(GSEAres_HALLMARK[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes_HALLMARK[[head(GSEAres_HALLMARK[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_HALLMARK[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres_HALLMARK, file = "results/GSEA_pt4_BL_REL_MPP_HALLMARK.rds")
data.table::fwrite(GSEAres_HALLMARK, file ="results/GSEA_pt4_BL_REL_MPP_HALLMARK.tsv" , sep = "\t", sep2 = c("", " ", ""))

#Visualise
GSEAres_HALLMARK$adjPvalue <- ifelse(GSEAres_HALLMARK$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_HALLMARK_sig=subset(GSEAres_HALLMARK, GSEAres_HALLMARK$padj<=0.1)
ggplot(GSEAres_HALLMARK, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA (p<0.1)")

