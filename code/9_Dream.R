#Extract the raw counts from the Seurat Object
df<-readRDS("data/4_Seurat_obj_CNA_referenceannotations.rds")
library(tidyverse)
library(Seurat)
#As I want to analyse only Pt11 BL and REL I make a new dataframe
#Create a subset for patient 11 BL and REL
pt11<-subset(x=df, subset=Patient=='pt11')
#Extract metadata
data11<-pt11@meta.data
#Create variable marking the dominant clone at BL and REL
data11<-mutate(data11, dominant = if_else(
  pt_status == 'pt11 BL'&clone.y %in% c("NDI", "NDIN")|
    pt_status == 'pt11 REL' & clone.y == 'NDINF', "1","0"))
#Check whether it worked
table(data11$dominant, data11$clone.y, data11$Sample)
#Add metadata back to Seurat object
pt11<-AddMetaData(object = pt11, metadata = data11, col.name = 'dominant')
pt11_BL_REL<-subset(x=pt11, pt_status %in% c("pt11 BL","pt11 REL"))
##Create a subset with only dominant population which is LMPP in this case
pt11_BL_REL<-subset(x=pt11_BL_REL, subset=predicted_CellType=='LMPP')
##From this subset select the dominant clones
pt11_BL_REL <- subset(x=pt11_BL_REL, dominant == "1")

#Get the metadata table
metadata<-pt11_BL_REL@meta.data
metadata$Cell<-rownames(metadata)
#I already have a conditions variable, so I can split right away
REL<-metadata[which(metadata$Sample == "REL"),]
BL<-metadata[which(metadata$Sample == "BL"),]
#Cells have already been filtered, so I guess I don't have to

#Next get the raw count matrix
#Get the count matrix
cts<-pt11_BL_REL@assays$RNA$counts
cts <- as.matrix(cts)
#Filter out ERCC
ERCC.genes<-grep("^ERCC-", rownames(cts), value = T)
#-> are not included anymore

cellsREL.m<-cts[,colnames(cts) %in% row.names(REL)]
cellsBL.m<-cts[,colnames(cts) %in% row.names(BL)]

MeanREL<-Matrix::rowMeans(cellsREL.m)
MeanBL<-Matrix::rowMeans(cellsBL.m)
Foldchange<-MeanREL/MeanBL

cellsREL.bi <- cellsREL.m
cellsREL.bi[cellsREL.bi >0] <- 1
REL.freq<-Matrix::rowSums(cellsREL.bi)

cellsBL.bi <- cellsBL.m
cellsBL.bi[cellsBL.bi >0] <- 1
BL.freq<-Matrix::rowSums(cellsBL.bi)

ExpFractionREL=REL.freq/ncol(cellsREL.m)
ExpFractionBL=BL.freq/ncol(cellsBL.m)

z0<-data.frame(Gene=rownames(cts),
               ExpREL=MeanREL,
               ExpBL=MeanBL,
               FoldChange=Foldchange,
               ExpFreqREL=REL.freq,
               ExpFreqBL=BL.freq,
               TotalREL=ncol(cellsREL.m),
               TotalBL=ncol(cellsBL.m),
               ExpFractionREL=ExpFractionREL,
               ExpFractionBL=ExpFractionBL)
rm(cellsREL.m, cellsBL.m, cellsREL.bi, cellsBL.bi)

# Filter out genes expressed in fewer cells than min.expFraction
z0<-subset(z0,ExpFractionREL >=1 | ExpFractionBL >=1)
cts<-cts[rownames(cts) %in% as.character(z0$Gene),]

# Normalize using scran size factors
#this is  absolutely necessary
# Get scran size factors
sf <- metadata$sizeFactor

nsf <- log(sf/Matrix::colSums(cts))
nsf <- exp(nsf - mean(nsf, na.rm=T))

dge <- edgeR::DGEList(counts = cts,
                      lib.size = Matrix::colSums(cts),
                      norm.factors = nsf,
                      remove.zeros = FALSE)

#Make a conditions dataframe
conditions.df <- data.frame(Cell = metadata$Cell,
                            condition = metadata$Sample,
                            covariate = metadata$Sort)

param = BiocParallel::SnowParam(2, "SOCK", progressbar=TRUE)

# estimate weights using linear mixed model of dream
print("Estimating weights using linear mixed model of dream")
formula <- ~ Sample
vobjDream = variancePartition::voomWithDreamWeights( dge, ~condition, conditions.df , BPPARAM=param )


# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
print("Fitting dream model on each gene")
fitmm = variancePartition::dream( vobjDream, ~condition, conditions.df , BPPARAM=param )
fitmm = variancePartition::eBayes(fitmm)

# Examine design matrix
print(head(fitmm$design, 3))

# Get results of hypothesis test on coefficients of interest
resT <- variancePartition::topTable( fitmm, coef=2, number=Inf )
resT <- tibble::rownames_to_column(resT, "Gene")

resT <- resT %>% dplyr::rename(p.value = P.Value, 
                               adjusted.pval = adj.P.Val, 
                               log2FC = logFC)
z0 <- dplyr::left_join(z0, resT,by = "Gene")

z0$ranking<- (-log10(z0$p.value))*sign(z0$log2FC)
z0<-z0 %>% dplyr::arrange(-ranking)

write_csv(z0, "results/Pt11/DREAM_LMPP.csv")

#Explore the data

#Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(z0,
                lab = z0$Gene,
                x = 'log2FC',
                y = 'adjusted.pval')

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
out_path <- "results/Pt11"
in_path<-"results/Pt11"

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
df_new <- read.csv("results/Pt11/DREAM_LMPP.csv", row.names = 1)
## 2. Prepare background genes #leave this out if already done -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- rownames(df_new)
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [12], GO = [54], reactome = [18], wikipathways [22], misigdb [70], this script uses KEGG, GO and Misigdb
#GO pathways
bg_genes_GO <- prepare_gmt(gmt_files[54], my_genes, savefile = FALSE)

#Prepare ranking of your DF gene
rankings<- sign(df_new$log2FC)*(-log10(df_new$adjusted.pval)) #use signed p values as preferred
names(rankings)<-rownames(df_new)
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
write_csv(GSEAres_GO, "results/Pt11?GSEA_LMPP_pt11_GO_Dream.csv")

#Order by pvalue
head(GSEAres_GO[order(pval),])
sum(GSEAres_GO[,padj<0.05])
sum(GSEAres_GO[,pval<0.01])


# plot the most significantly enriched pathway
library(ggplot2)
plotEnrichment(bg_genes_GO[[head(GSEAres_GO[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_GO[order(padj), ], 1)$pathway)
GSEAres_GO$adjPvalue <- ifelse(GSEAres_GO$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_GO_sig=subset(GSEAres_GO, GSEAres_GO$padj<=0.05)
ggplot(GSEAres_GO_sig, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO pathways Enrichment Score from GSEA (p<0.05)")

GSEAres_GO_sig=subset(GSEAres_GO, GSEAres_GO$padj<=0.1)
ggplot(GSEAres_GO_sig, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO pathways Enrichment Score from GSEA (p<0.1)")

#Save
## 5. Save the results -----------------------------------------------
saveRDS(GSEAres_GO, file = "results/Pt11/GSEA_pt11_BL_REL_LMPP_GO_dream.rds")
data.table::fwrite(GSEAres_GO, file ="results/GSEA_pt11_BL_REL_LMPP_GO_dream.tsv" , sep = "\t", sep2 = c("", " ", ""))


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
write_csv(GSEAres,"results/GSEA_pt11_LMPP_wiki.csv")

#Order by pvalue
head(GSEAres[order(pval),])
sum(GSEAres[,padj<0.05])
sum(GSEAres[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres, file = "results/GSEA_pt11_BL_REL_LMPP_wiki.rds")
data.table::fwrite(GSEAres, file ="results/GSEA_pt11_BL_REL_LMPP_wiki.tsv" , sep = "\t", sep2 = c("", " ", ""))

######KEGG pathways
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [12], GO = [54], reactome = [18], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes_KEGG <- prepare_gmt(gmt_files[12], my_genes, savefile = FALSE)

#Run GSEA
GSEAres_KEGG <- fgsea(pathways = bg_genes_KEGG, # List of gene sets to check
                      stats = rankings,
                      scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                      minSize = 10,
                      maxSize = 500,
                      nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres_KEGG)
write_csv(GSEAres_KEGG,"results/GSEA_pt11_LMPP_KEGG.csv")

#Order by pvalue
head(GSEAres_KEGG[order(pval),])
sum(GSEAres_KEGG[,padj<0.05])
sum(GSEAres_KEGG[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes_KEGG[[head(GSEAres_KEGG[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_KEGG[order(padj), ], 1)$pathway)

#Barplot
GSEAres_KEGG$adjPvalue <- ifelse(GSEAres_KEGG$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_KEGG_sig=subset(GSEAres_KEGG, GSEAres_KEGG$padj<=0.05)
ggplot(GSEAres_KEGG_sig , aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways Enrichment Score from GSEA (p<0.1)")
#Save
saveRDS(GSEAres, file = "results/GSEA_pt11_BL_REL_LMPP_KEGG.rds")
data.table::fwrite(GSEAres, file ="results/GSEA_pt11_BL_REL_LMPP_KEGG.tsv" , sep = "\t", sep2 = c("", " ", ""))

######HALLMARK
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes_HALLMARK <- prepare_gmt(gmt_files[68], my_genes, savefile = FALSE)

#Run GSEA
GSEAres_HALLMARK <- fgsea(pathways = bg_genes_HALLMARK, # List of gene sets to check
                          stats = rankings,
                          scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                          minSize = 10,
                          maxSize = 500,
                          nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres_HALLMARK)
write_csv(GSEAres_HALLMARK,"results/Pt11/GSEA_pt11_LMPP_HALLMARK_dream.csv")

#Order by pvalue
head(GSEAres_HALLMARK[order(pval),])
sum(GSEAres_HALLMARK[,padj<0.05])
sum(GSEAres_HALLMARK[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes_HALLMARK[[head(GSEAres_HALLMARK[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_HALLMARK[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres_HALLMARK, file = "results/Pt11/GSEA_pt11_BL_REL_LMPP_HALLMARK_dream.rds")
data.table::fwrite(GSEAres_HALLMARK, file ="results/GSEA_pt11_BL_REL_LMPP_HALLMARK_dream.tsv" , sep = "\t", sep2 = c("", " ", ""))
#Visualise
GSEAres_HALLMARK$adjPvalue <- ifelse(GSEAres_HALLMARK$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_HALLMARK_sig=subset(GSEAres_HALLMARK, GSEAres_HALLMARK$padj<0.05)
ggplot(GSEAres_HALLMARK_sig , aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA (p<0.05)")
GSEAres_HALLMARK_sig=subset(GSEAres_HALLMARK, GSEAres_HALLMARK$padj<0.1)
ggplot(GSEAres_HALLMARK , aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")

#All pathways
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes_all <- prepare_gmt(gmt_files[5], my_genes, savefile = FALSE)

#Run GSEA
GSEAres_all <- fgsea(pathways = bg_genes_all, # List of gene sets to check
                     stats = rankings,
                     scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                     minSize = 10,
                     maxSize = 500,
                     nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres_all)
write_csv(GSEAres_all,"results/Pt11_GSEA_pt11_LMPP_pathways_dream.csv")

#Order by pvalue
head(GSEAres_all[order(pval),])
sum(GSEAres_all[,padj<0.05])
sum(GSEAres_all[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes_all[[head(GSEAres_all[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_all[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres_all, file = "results/Pt11/GSEA_pt11_BL_REL_LMPP_all_dream.rds")
data.table::fwrite(GSEAres_all, file ="results/Pt11/GSEA_pt11_BL_REL_LMPP_all_dream.tsv" , sep = "\t", sep2 = c("", " ", ""))
#Visualise
GSEAres_all$padj <- ifelse(GSEAres_all$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_all_sig=subset(GSEAres_all, GSEAres_all$padj<=0.1)
ggplot(GSEAres_all , aes(reorder(pathway, NES), NES, fill = padj)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="All pathways Enrichment Score from GSEA")


############################################################################################################
#Setting the expression limit at 1 eliminates a huge number of cells, therefore I try with a limit of 0.5
#Extract the raw counts from the Seurat Object
df<-readRDS("data/4_Seurat_obj_CNA_referenceannotations.rds")
library(tidyverse)
library(Seurat)
#As I want to analyse only Pt11 BL and REL I make a new dataframe
#Create a subset for patient 11 BL and REL
pt11<-subset(x=df, subset=Patient=='pt11')
#Extract metadata
data11<-pt11@meta.data
#Create variable marking the dominant clone at BL and REL
data11<-mutate(data11, dominant = if_else(
  pt_status == 'pt11 BL'&clone.y %in% c("NDI", "NDIN")|
    pt_status == 'pt11 REL' & clone.y == 'NDINF', "1","0"))
#Check whether it worked
table(data11$dominant, data11$clone.y, data11$Sample)
#Add metadata back to Seurat object
pt11<-AddMetaData(object = pt11, metadata = data11, col.name = 'dominant')
pt11_BL_REL<-subset(x=pt11, pt_status %in% c("pt11 BL","pt11 REL"))
##Create a subset with only dominant population which is LMPP in this case
pt11_BL_REL<-subset(x=pt11_BL_REL, subset=predicted_CellType=='LMPP')
##From this subset select the dominant clones
pt11_BL_REL <- subset(x=pt11_BL_REL, dominant == "1")

#Get the metadata table
metadata<-pt11_BL_REL@meta.data
metadata$Cell<-rownames(metadata)
#I already have a conditions variable, so I can split right away
REL<-metadata[which(metadata$Sample == "REL"),]
BL<-metadata[which(metadata$Sample == "BL"),]
#Cells have already been filtered, so I guess I don't have to

#Next get the raw count matrix
#Get the count matrix
cts<-pt11_BL_REL@assays$RNA$counts
cts <- as.matrix(cts)
#Filter out ERCC
ERCC.genes<-grep("^ERCC-", rownames(cts), value = T)
#-> are not included anymore

cellsREL.m<-cts[,colnames(cts) %in% row.names(REL)]
cellsBL.m<-cts[,colnames(cts) %in% row.names(BL)]

MeanREL<-Matrix::rowMeans(cellsREL.m)
MeanBL<-Matrix::rowMeans(cellsBL.m)
Foldchange<-MeanREL/MeanBL

cellsREL.bi <- cellsREL.m
cellsREL.bi[cellsREL.bi >0] <- 1
REL.freq<-Matrix::rowSums(cellsREL.bi)

cellsBL.bi <- cellsBL.m
cellsBL.bi[cellsBL.bi >0] <- 1
BL.freq<-Matrix::rowSums(cellsBL.bi)

ExpFractionREL=REL.freq/ncol(cellsREL.m)
ExpFractionBL=BL.freq/ncol(cellsBL.m)

z0<-data.frame(Gene=rownames(cts),
               ExpREL=MeanREL,
               ExpBL=MeanBL,
               FoldChange=Foldchange,
               ExpFreqREL=REL.freq,
               ExpFreqBL=BL.freq,
               TotalREL=ncol(cellsREL.m),
               TotalBL=ncol(cellsBL.m),
               ExpFractionREL=ExpFractionREL,
               ExpFractionBL=ExpFractionBL)
rm(cellsREL.m, cellsBL.m, cellsREL.bi, cellsBL.bi)

# Filter out genes expressed in fewer cells than min.expFraction, here I change the cutoff to 0.5
z0<-subset(z0,ExpFractionREL >=0.5 | ExpFractionBL >=0.5)
cts<-cts[rownames(cts) %in% as.character(z0$Gene),]

# Normalize using scran size factors
#this is  absolutely necessary
# Get scran size factors
sf <- metadata$sizeFactor

nsf <- log(sf/Matrix::colSums(cts))
nsf <- exp(nsf - mean(nsf, na.rm=T))

dge <- edgeR::DGEList(counts = cts,
                      lib.size = Matrix::colSums(cts),
                      norm.factors = nsf,
                      remove.zeros = FALSE)

#Make a conditions dataframe
conditions.df <- data.frame(Cell = metadata$Cell,
                            condition = metadata$Sample,
                            covariate = metadata$Sort)

param = BiocParallel::SnowParam(2, "SOCK", progressbar=TRUE)

# estimate weights using linear mixed model of dream
vobjDream = variancePartition::voomWithDreamWeights( dge, ~condition, conditions.df , BPPARAM=param )


# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
print("Fitting dream model on each gene")
fitmm = variancePartition::dream( vobjDream, ~condition, conditions.df , BPPARAM=param )
fitmm = variancePartition::eBayes(fitmm)

# Examine design matrix
print(head(fitmm$design, 3))

# Get results of hypothesis test on coefficients of interest
resT <- variancePartition::topTable( fitmm, coef=2, number=Inf )
resT <- tibble::rownames_to_column(resT, "Gene")

resT <- resT %>% dplyr::rename(p.value = P.Value, 
                               adjusted.pval = adj.P.Val, 
                               log2FC = logFC)
z0 <- dplyr::left_join(z0, resT,by = "Gene")

z0$ranking<- (-log10(z0$p.value))*sign(z0$log2FC)
z0<-z0 %>% dplyr::arrange(-ranking)

write_csv(z0, "results/Pt11/Pt11_DREAM_LMPP_exp0.5.csv")

#Explore the data

#Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(z0,
                lab = z0$Gene,
                x = 'log2FC',
                y = 'adjusted.pval')

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
out_path <- "results/Pt11"
in_path<-"results/Pt11"

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
df_new <- read.csv("results/Pt11/Pt11_DREAM_LMPP_exp0.5.csv", row.names = 1)
## 2. Prepare background genes #leave this out if already done -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
my_genes <- rownames(df_new)
list.files(bg_path)
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [12], GO = [54], reactome = [18], wikipathways [22], misigdb [70], this script uses KEGG, GO and Misigdb
#GO pathways
bg_genes_GO <- prepare_gmt(gmt_files[54], my_genes, savefile = FALSE)

#Prepare ranking of your DF gene
rankings<- sign(df_new$log2FC)*(-log10(df_new$adjusted.pval)) #use signed p values as preferred
names(rankings)<-rownames(df_new)
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
write_csv(GSEAres_GO, "results/Pt11/GSEA_LMPP_pt11_GO_Dream_exp0.5.csv")

#Order by pvalue
head(GSEAres_GO[order(pval),])
sum(GSEAres_GO[,padj<0.05])
sum(GSEAres_GO[,pval<0.01])


# plot the most significantly enriched pathway
library(ggplot2)
plotEnrichment(bg_genes_GO[[head(GSEAres_GO[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_GO[order(padj), ], 1)$pathway)
GSEAres_GO$adjPvalue <- ifelse(GSEAres_GO$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_GO_sig=subset(GSEAres_GO, GSEAres_GO$padj<=0.05)
ggplot(GSEAres_GO_sig, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO pathways Enrichment Score from GSEA (p<0.05)")

GSEAres_GO_sig=subset(GSEAres_GO, GSEAres_GO$padj<=0.1)
ggplot(GSEAres_GO_sig, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO pathways Enrichment Score from GSEA (p<0.1)")

#Save
## 5. Save the results -----------------------------------------------
saveRDS(GSEAres_GO, file = "results/Pt11/GSEA_pt11_BL_REL_LMPP_GO_dream_exp0.5.rds")
data.table::fwrite(GSEAres_GO, file ="results/GSEA_pt11_BL_REL_LMPP_GO_dream_exp0.5.tsv" , sep = "\t", sep2 = c("", " ", ""))


######HALLMARK
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes_HALLMARK <- prepare_gmt(gmt_files[68], my_genes, savefile = FALSE)

#Run GSEA
GSEAres_HALLMARK <- fgsea(pathways = bg_genes_HALLMARK, # List of gene sets to check
                          stats = rankings,
                          scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                          minSize = 10,
                          maxSize = 500,
                          nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres_HALLMARK)
write_csv(GSEAres_HALLMARK,"results/Pt11/GSEA_pt11_LMPP_HALLMARK_dream_exp0.5.csv")

#Order by pvalue
head(GSEAres_HALLMARK[order(pval),])
sum(GSEAres_HALLMARK[,padj<0.05])
sum(GSEAres_HALLMARK[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes_HALLMARK[[head(GSEAres_HALLMARK[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_HALLMARK[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres_HALLMARK, file = "results/Pt11/GSEA_pt11_BL_REL_LMPP_HALLMARK_dream_exp0.5.rds")
data.table::fwrite(GSEAres_HALLMARK, file ="results/GSEA_pt11_BL_REL_LMPP_HALLMARK_dream_exp0.5.tsv" , sep = "\t", sep2 = c("", " ", ""))
#Visualise
GSEAres_HALLMARK$adjPvalue <- ifelse(GSEAres_HALLMARK$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_HALLMARK_sig=subset(GSEAres_HALLMARK, GSEAres_HALLMARK$padj<0.05)
ggplot(GSEAres_HALLMARK_sig , aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA (p<0.05)")
GSEAres_HALLMARK_sig=subset(GSEAres_HALLMARK, GSEAres_HALLMARK$padj<0.1)
ggplot(GSEAres_HALLMARK , aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")

#All pathways
gmt_files <- list.files(path = bg_path, pattern = '.gmt', full.names = TRUE)
gmt_files
#according to this list KEGG = [11], GO = [53], reactome = [17], wikipathways [21], misigdb [69], this script uses KEGG, GO and Misigdb
bg_genes_all <- prepare_gmt(gmt_files[5], my_genes, savefile = FALSE)

#Run GSEA
GSEAres_all <- fgsea(pathways = bg_genes_all, # List of gene sets to check
                     stats = rankings,
                     scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                     minSize = 10,
                     maxSize = 500,
                     nproc = 1) # for parallelisation
#Check whether it worked
head(GSEAres_all)
write_csv(GSEAres_all,"results/Pt11_GSEA_pt11_LMPP_pathways_dream_exp0.5.csv")

#Order by pvalue
head(GSEAres_all[order(pval),])
sum(GSEAres_all[,padj<0.05])
sum(GSEAres_all[,pval<0.01])


# plot the most significantly enriched pathway
plotEnrichment(bg_genes_all[[head(GSEAres_all[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres_all[order(padj), ], 1)$pathway)
#Save
saveRDS(GSEAres_all, file = "results/Pt11/GSEA_pt11_BL_REL_LMPP_all_dream_exp0.5.rds")
data.table::fwrite(GSEAres_all, file ="results/Pt11/GSEA_pt11_BL_REL_LMPP_all_dream_exp0.5.tsv" , sep = "\t", sep2 = c("", " ", ""))
#Visualise
GSEAres_all$adjPvalue <- ifelse(GSEAres_all$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
GSEAres_all_sig=subset(GSEAres_all, GSEAres_all$padj<=0.01)
ggplot(GSEAres_all_sig , aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="All pathways Enrichment Score from GSEA")

#Let's determine how much Deseq2 and dream overlap
library(VennDetail)
library(readxl)
library(tidyverse)
library(utils)

dream<-read_csv("results/Pt11/Pt11_DREAM_LMPP_exp0.5.csv")
deseq<-read_csv("results/Pt11/DEA/LMPP/DEA_pt11_LMPP_padj.csv")

venn<-venndetail(list(dream=dream$Gene, deseq= deseq$gene))
plot(venn)


