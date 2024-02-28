# Felix A. Radtke and Sven Turkalj
#### Reformat counts matrix #####

# Load libraries
library(tidyverse)
library(stringr)
library(readxl)


# TODO SVEN: Please put in the exact code to construct IVO_VEN_REL_all_cells_clone_index_oligodT.csv from the genotyping meta here
# Generate the metadata sheet that will be used for reformatting of the IVO_VEN relapse patients matrix

# Merge all cell/clone metadata sheets - this metadata contains only cells that are successfully genotyped - this is done only once

pt4 <- read_csv("data/raw data per patient/pt4_clones_population.csv") %>% 
  select(cell_id, Sample, clone, Patient, Population)
pt9 <- read_csv("data/raw data per patient/pt9_clones_population.csv") %>% 
  select(cell_id, Sample, clone, Patient, Population)
pt11 <- read_csv("data/raw data per patient/pt11_clones_population.csv") %>% 
  select(cell_id, Sample, clone, Patient, Population)
pt14 <- read_csv("data/raw data per patient/pt14_clones_population.csv") %>% 
  select(cell_id, Sample, clone, Patient, Population)
pt15 <- read_csv("data/raw data per patient/pt15_clones_population.csv") %>% 
  select(cell_id, Sample, clone, Patient, Population)

plate_indexes <- rbind(pt4, pt9, pt11, pt14, pt15) %>% 
  select(cell_id, clone) # This is a list of all cells where genotyping did not completely fail (at least one locus was detected)

all_cells <- read_csv("data/raw data per patient/metadata_all_cells.csv") # This contains population information for all 14247 wells in which a cell was sorted


all_cells <- all_cells %>%
  separate(Sample, into = c("Patient", "Sample"), sep = "-") 
table(all_cells$Patient)

#I still need a column with pt and sample, so I'll create that now
all_cells$pt_status <- paste(all_cells$Patient,all_cells$Sample)
table(all_cells$pt_status)

all <- left_join(all_cells, plate_indexes, by = "cell_id") # Now we have clone + population information for all cells

all <- all %>% 
  mutate(new_clone = case_when(
    Sample == "control" ~ "WT",
    TRUE ~ clone
  )) # plate_indexes did not contain clone info about WT control cells so we add it here

all$clone <- NULL

all <- all %>% 
  dplyr::mutate(new_clone = replace_na(new_clone, "Undetermined")) %>% 
  dplyr::rename(clone = new_clone) %>% 
  write_csv("results/IVO_VEN_Relapse_all_cells_clone_index.csv")

# Now we have to add A1_A2 and A3_A4 info for each plate (based on sample sheets made by Bilyana)

barcodes <- read_csv("data/raw data per patient/IVO_VEN_REL_plate_oligodT_A1_to_A4.csv") # Contains A1_A2 or A3_A4 info for each of the 75 sorted plates

all$OligodT_barcode_plate <- barcodes$OligodT_barcode_plate[match(all$Plate, barcodes$Plate)]  #add oligodT info based on plate

# Moreover, we have to add the well name for each cell

all <- all %>% 
  mutate(Sample_well = sapply(strsplit(cell_id, "-"), `[`, 2)) # Takes the string between the first and second dash - this is the well

all$Sample <- gsub('_', '-', all$Sample) # Some samples had underscores - we convert them all to dashes

write_csv(all, "results/IVO_VEN_REL_all_cells_clone_index_oligodT.csv")

# Now, we have everything ready for matrix reformatting! 


plate_indexes <- read_csv("results/IVO_VEN_REL_all_cells_clone_index_oligodT.csv")

plate_indexes <- all

# Import list of OligodT barcodes
oligodt_barcodes <- read_excel("data/raw data per patient/oligodT_384_barcodes.xlsx") %>%
  dplyr::select(OligodT_barcode_name, OligodT_barcode_plate, OligodT_barcode_sequence, Sample_well)

# Join the plate i7 and i7 indexes to the oligodT barcodes
# Note this creates a set of cells for each plate index
transcriptome_indexes <- plate_indexes %>% 
  left_join(oligodt_barcodes, by = c("Sample_well", "OligodT_barcode_plate")) # Now, for each cell we have the information about the oligodT barcode sequence and name, based on the well and oligodT plate info

#### READ IN FILES FROM RNASEQ PIPELINE #####

# Read in the features (genes) file - this will be from the STARSOLO folder! 
features_file <- read_tsv("data/raw data per patient/features.tsv.gz",
                          col_names = c("gene_id", "gene_name", "assay"))

# Read in the tsv file with cell barcodes (one for each pipeline run)
barcodes_file <- read_tsv("data/raw data per patient/barcodes.tsv.gz", col_names = "barcode")

# At this stage, the cell names inside the matrix need to be reformatted to be names after cell_id within transcriptome_indexes. This is done via linkage of HT sequence info. 
# Some additional fiddling is required due to some initial fastq file names being slighlty different from the plate names in the metadata. 

transcriptome_indexes <- transcriptome_indexes %>%
  mutate(Plate_ID = str_extract(Plate, "PL[0-9]+")) 

barcodes_file_adjusted <- barcodes_file %>%
  mutate(new_barcode = str_replace_all(barcode, "-([0-9]+)(?![0-9PL])", "-PL\\1")) %>%
  mutate(new_barcode = str_replace_all(new_barcode, "P([0-9])", "PL\\1")) %>%
  mutate(new_barcode = str_replace_all(new_barcode, "(PL)(?<![0-9])([0-9])(?![0-9])", paste0("\\1", "0", "\\2"))) %>%
  mutate(Sort = str_extract(new_barcode, "S[0-9]+")) %>%
  rowwise() %>%
  mutate(plates = list(str_match_all(new_barcode, "PL[0-9]+"))) %>%
  mutate(Plate_1 = plates[[1]][1], Plate_2 = plates[[1]][2]) %>%
  ungroup() %>%
  mutate(OligodT_barcode_name = str_extract(new_barcode, "HT[0-9]+"))

join1 <- left_join(barcodes_file_adjusted, transcriptome_indexes, by = c("Sort" = "Sort", "Plate_1" = "Plate_ID", "OligodT_barcode_name"))
join2 <- left_join(barcodes_file_adjusted, transcriptome_indexes, by = c("Sort" = "Sort", "Plate_2" = "Plate_ID", "OligodT_barcode_name"))

transcriptome_indexes_merged <- join1 %>%
  mutate(across(everything(), ~ coalesce(.x, join2[[cur_column()]])))

if(sum(!is.na(transcriptome_indexes_merged$cell_id)) == nrow(transcriptome_indexes)){
  "Merge successful"
}else{
  stop("Merge problem")
}

### Reformat the RNA-seq counts to make the reformatted counts matrix ####

# Read in the tsv file with RNA-seq counts: (needs to be loaded from your PC because file too big for github)
raw_counts_file <- as(Matrix::readMM("data/raw data per patient/matrix.mtx.gz"), "CsparseMatrix") #data is in the raw data file in Ivo_Ven onedrive

# Function for adding cell barcodes and gene names
format_counts_matrix <- function(matrix, barcodes, features){
  if (dim(matrix)[2] != length(barcodes$new_barcode)) {
    stop("Mismatch dimension between barcodes file and the matrix file")
  }else{
    print("Cell numbers match :)")
    colnames(matrix) <- barcodes$cell_id
  }
  
  if (dim(matrix)[1] != length(features$gene_name)) {
    stop("Mismatch dimension between gene file and the matrix file")
  }else{
    print("Gene numbers match :)")
    rownames(matrix) <- make.unique(features$gene_name)
  }
  return(matrix)
}

# Add cell barcodes and gene names
if(!identical(transcriptome_indexes_merged$barcode, barcodes_file$barcode)){
  stop("Assembly failed, carefully check")
}
matrix <- format_counts_matrix(raw_counts_file, transcriptome_indexes_merged, features_file)

# This matrix has some NA columns, those correspond to unwanted "cells". Remove them:
dim(matrix)
matrix <- matrix[, !is.na(colnames(matrix))]
dim(matrix)

# Save matrix to assemble 
saveRDS(matrix, "data/0_matrix_reformatted.rds")
saveRDS(transcriptome_indexes_merged, "data/0_transcriptome_indexes.rds")
