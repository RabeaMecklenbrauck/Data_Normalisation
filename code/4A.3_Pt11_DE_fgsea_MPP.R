#This script uses the normalised data and performs DEA (DeSeq2) and fgsea
#Open the annotated Seurat object
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
##Create a subset with only dominant population which is MPP in this case, I decided to combine bith MPP populations
pt11_BL_REL<-subset(x=pt11_BL_REL, subset=predicted_CellType %in% c("MPP-MkEry","MPP-MyLy"))
##From this subset select the dominant clones
pt11_BL_REL <- subset(x=pt11_BL_REL, dominant == "1")
#Check whether the subsetting worked
table(pt11_BL_REL$predicted_CellType)
table(pt11_BL_REL$dominant)

#This leaves you with 143 MPP-MkEry and 54 MPP-MyLy

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
write.csv(res_tbl_all,"results/DEA_pt11_MPP.csv" )

#Filter for only significant Genes
#Set threshold
padj_cutoff <- 0.05
#Subset for signifcant results 
sig_res_all<- subset(res_tbl_all, res_tbl_all$padj<0.05)
write.csv(sig_res_all, "results/DEA_pt11_MPP_padj.csv")


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
out_path <- "results/"
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
df_new <- read.csv("results/DEA_pt11_MPP_padj.csv", row.names = 1)
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
bg_genes <- prepare_gmt(gmt_files[53], my_genes, savefile = FALSE)

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
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres)
write_csv(GSEAres, "results/GSEA_MPP_pt11_GO.csv")

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
saveRDS(GSEAres, file = "results/GSEA_pt11_BL_REL_MPP_GO.rds")
data.table::fwrite(GSEAres, file ="results/GSEA_pt11_BL_REL_LMPP_GO.tsv" , sep = "\t", sep2 = c("", " ", ""))

######wiki pathways
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes <- prepare_gmt(gmt_files[21], my_genes, savefile = FALSE)

#Run GSEA
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres)
write_csv(GSEAres,"results/GSEA_pt11_MPP_wiki.csv")

#Order by pvalue
head(GSEAres[order(pval),])
sum(GSEAres[,padj<0.05])
sum(GSEAres[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres, file = "results/GSEA_pt11_BL_REL_MPP_wiki.rds")
data.table::fwrite(GSEAres, file ="results/GSEA_pt11_BL_REL_MPP_wiki.tsv" , sep = "\t", sep2 = c("", " ", ""))

######KEGG pathways
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes <- prepare_gmt(gmt_files[11], my_genes, savefile = FALSE)

#Run GSEA
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres)
write_csv(GSEAres,"results/GSEA_pt11_MPP_KEGG.csv")

#Order by pvalue
head(GSEAres[order(pval),])
sum(GSEAres[,padj<0.05])
sum(GSEAres[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres, file = "results/GSEA_pt11_BL_REL_MPP_KEGG.rds")
data.table::fwrite(GSEAres, file ="results/GSEA_pt11_BL_REL_LMPP_KEGG.tsv" , sep = "\t", sep2 = c("", " ", ""))
