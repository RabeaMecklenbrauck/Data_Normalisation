#Prepare the count matrix to run in SCENIC
library(Seurat)
library(tidyverse)
library(magrittr)
library(SingleCellExperiment)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC@v1.1.2")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
df.loom<-as.loom(df)
df.loom
#Get expression matrix
library(SCopeLoomR)
exprMat<-get_dgem(df.loom)
SaveLoom(df, "data/5_Loom_obj_annotated.loom", overwrite = FALSE, verbose = TRUE)
?SaveLoom
