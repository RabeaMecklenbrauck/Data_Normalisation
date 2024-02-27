#This script should annotate the gene list Asger did to my data
#First read in the gmt.file 
#For that I need a package to read gmt files
ref<-cogena::gmt2list("data/human.hematopoiesis.signature.genes.curated.v6.gmt")

#I'm trying to follow the Star protocol 
#First there are a lot of packages to install, that takes a while...
if(!require(devtools)) {
  install.packages("devtools")
   }
if(!require(BiocManager)) {
   install.packages("BiocManager")
   }
  devtools::install_github('supatt-lab/SingCellaR',ref='master',
                           
                           repos=BiocManager::repositories())
  
  #Install the requored python modules
  library(reticulate)
  use_python("/miniconda3/envs/r-reticulate/bin/python")
  conda_create("r-reticulate", python_version = "3.8")
  py_install("fa2", envname="r-reticulate")
  py_install("networkx", envname = "r-reticulate")
 # py_install("Scrublet", envname = "r-reticulate") #that would be for doublet removal, I don't need that
#Install the required r packages, check with if require if the packages are already loaded
  if(!require(harmony)){
    install.packages("harmony")
  }
if(!require(AUCell)){
  BiocManager::install(("AUCell"))
}  
install.packages("doParallel")  
install.packages("doRNG")
install.packages("DAseq")
install.packages("destiny")
devtools::install_github('cole-trapnell-lab/monocle3',ref="develop")
install_github("https://github.com/theislab/destiny",build_vignettes=FALSE, dependencies=TRUE)
#install scanorama in th e pythin envirnoment using conda pip install scanorama
install.packages("Seurat")
install.packages ("rliger")
BiocManager::install("sva")
BiocManager::install("limma")

#Load the data frame 
df<-readRDS("data/2_lognorm_seurat.rds")

#Convert the data frame to single cell experiment
install.packages("Seurat")
BiocManager::install("scater")

library(SingCellaR)
library(Seurat)
FindVariableFeatures(df)
df.sce<-as.SingleCellExperiment(df)
df_SingCellaR<-new("SingCellaR", df.sce)
process_cells_annotation(df_SingCellaR, mito_genes_start_with = "MT-")
plot_cells_annotation(df_SingCellaR, type = "histogram") #this also includes UMIS which this data set does not have -> ignore
SingCellaR::runPCA(df_SingCellaR, use.components =  100, use.regressout.data =FALSE)
?SingCellaR::RunPCA
