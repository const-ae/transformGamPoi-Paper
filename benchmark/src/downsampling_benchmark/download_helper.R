library(tidyverse)
library(SingleCellExperiment)


data_folder <- "/scratch/ahlmanne/benchmark/data"

load_mcSCRBseq_data <- function(){
  file_loc <- file.path(data_folder, "GSE103568_JM8_UMIcounts.txt.gz")
  if(! file.exists(file_loc)){
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103568&format=file&file=GSE103568%5FJM8%5FUMIcounts%2Etxt%2Egz", file_loc)
  }
  raw_data <- as.matrix(read.delim(file_loc))
  SingleCellExperiment(list(counts = raw_data), rowData = DataFrame(gene_id = rownames(raw_data)))
}

load_smartSeq3_fibroblasts <- function(){
  if(! file.exists(file.path(data_folder, "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt"))){
    if(! file.exists(file.path(data_folder, "smart_seq3_data.zip"))){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/E-MTAB-8735.processed.3.zip", file.path(data_folder, "smart_seq3_data.zip"))
    }
    unzip(file.path(data_folder, "smart_seq3_data.zip"), files = "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt", exdir = data_folder)
  }
  if(! file.exists(file.path(data_folder, "Smartseq3.Fibroblasts.sample_annotation.txt"))){
    if(! file.exists(file.path(data_folder, "smart_seq3_data_annotation.zip"))){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/E-MTAB-8735.processed.2.zip", file.path(data_folder, "smart_seq3_data_annotation.zip"))
    }
    unzip(file.path(data_folder, "smart_seq3_data_annotation.zip"), files = "Smartseq3.Fibroblasts.sample_annotation.txt", exdir = data_folder)
  }
  
  col_data <- read.delim(file.path(data_folder, "Smartseq3.Fibroblasts.sample_annotation.txt"))
  raw_data <- read.delim(file.path(data_folder, "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt")) %>%
    as.matrix() %>%
    as("dgCMatrix")
  
  rownames(col_data) <- col_data$BC
  col_data <- col_data[colnames(raw_data), ]
  stopifnot(all(colnames(raw_data) == col_data$BC))
  
  
  SingleCellExperiment(list(counts = raw_data), rowData = DataFrame(gene_id = rownames(raw_data)), 
                       colData = DataFrame(col_data))
}

load_smartSeq3_fibroblasts_alt <- function(){
  if(! file.exists(file.path(data_folder, "Fibroblasts.plate2.umis.ex.txt"))){
    if(! file.exists(file.path(data_folder, "E-MTAB-10148.processed.1.zip"))){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-10148/E-MTAB-10148.processed.1.zip",file.path(data_folder, "E-MTAB-10148.processed.1.zip"))
    }
    if(! file.exists(file.path(data_folder, "E-MTAB-10148.processed.2.zip"))){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-10148/E-MTAB-10148.processed.2.zip", file.path(data_folder, "E-MTAB-10148.processed.2.zip"))
    }
    unzip(file.path(data_folder, "E-MTAB-10148.processed.1.zip"), files = "Fibroblasts.plate1.umis.ex.txt", exdir = data_folder)
    unzip(file.path(data_folder, "E-MTAB-10148.processed.1.zip"), files = "Fibroblasts.plate2.umis.ex.txt", exdir = data_folder)
    unzip(file.path(data_folder, "E-MTAB-10148.processed.2.zip"), files = "Fibroblasts.plate1.annotation.txt", exdir = data_folder)
    unzip(file.path(data_folder, "E-MTAB-10148.processed.2.zip"), files = "Fibroblasts.plate2.annotation.txt", exdir = data_folder)
  }
  
  raw_data_1 <- read.delim(file.path(data_folder, "Fibroblasts.plate1.umis.ex.txt")) %>%
    as.matrix() %>%
    as("dgCMatrix")
  raw_data_2 <- read.delim(file.path(data_folder, "Fibroblasts.plate2.umis.ex.txt")) %>%
    as.matrix() %>%
    as("dgCMatrix")
  annot_1 <- read.delim(file.path(data_folder, "Fibroblasts.plate1.annotation.txt"))
  annot_2 <- read.delim(file.path(data_folder, "Fibroblasts.plate2.annotation.txt"))
  
  common_rows <- intersect(rownames(raw_data_1), rownames(raw_data_2))
  raw_data_1 <- raw_data_1[common_rows, ]
  raw_data_2 <- raw_data_2[common_rows, ]
  
  stopifnot(length(intersect(colnames(raw_data_1), colnames(raw_data_2))) == 0)
  rownames(annot_1) <- annot_1$BC
  annot_1 <- annot_1[colnames(raw_data_1), ]
  stopifnot(all(colnames(raw_data_1) == annot_1$BC))
  rownames(annot_2) <- annot_2$BC
  annot_2 <- annot_2[colnames(raw_data_2), ]
  stopifnot(all(colnames(raw_data_2) == annot_2$BC))
  annot_1 <- transmute(annot_1, plate = 1, quality_check = QC)
  annot_2 <- transmute(annot_2, plate = 2, quality_check = QC)
  
  cbind(SingleCellExperiment(list(counts = raw_data_1), rowData = DataFrame(gene_id = rownames(raw_data_1)), colData = annot_1),
        SingleCellExperiment(list(counts = raw_data_2), rowData = DataFrame(gene_id = rownames(raw_data_2)), colData = annot_2))
}

load_smartSeq3_hek <- function(){
  if(! file.exists(file.path(data_folder, "Smartseq3.HEK.cleanup.UMIcounts.txt"))){
    if(! file.exists(file.path(data_folder, "smart_seq3_data.zip"))){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/E-MTAB-8735.processed.3.zip", file.path(data_folder, "smart_seq3_data.zip"))
    }
    unzip(file.path(data_folder, "smart_seq3_data.zip"), files = "Smartseq3.HEK.cleanup.UMIcounts.txt", exdir = data_folder)
  }
  raw_data <- read.delim(file.path(data_folder, "Smartseq3.HEK.cleanup.UMIcounts.txt")) %>%
    as.matrix() %>%
    as("dgCMatrix")
  
  SingleCellExperiment(list(counts = raw_data), rowData = DataFrame(gene_id = rownames(raw_data)))
}

load_smartSeq3_siRNA_knockdown_experiment <- function(){
  if(! file.exists(file.path(data_folder, "smart_seq3_siRNA_knockdown_full.RDS"))){
    if(! file.exists(file.path(data_folder, "ss3_n4298_fibs_siKD_umiCast.rds"))){
      download.file("https://github.com/sandberg-lab/lncRNAs_bursting/raw/main/data/ss3_n4298_fibs_siKD_umiCast.rds", file.path(data_folder, "ss3_n4298_fibs_siKD_umiCast.rds"))
    }
    if(! file.exists(file.path(data_folder, "ss3_n4298_fibs_siKD_umiC57.rds"))){
      download.file("https://github.com/sandberg-lab/lncRNAs_bursting/raw/main/data/ss3_n4298_fibs_siKD_umiC57.rds", file.path(data_folder, "ss3_n4298_fibs_siKD_umiC57.rds"))
    }
    if(! file.exists(file.path(data_folder, "ss3_n4298_fibs_siKD_meta.rds"))){
      download.file("https://github.com/sandberg-lab/lncRNAs_bursting/raw/main/data/ss3_n4298_fibs_siKD_meta.rds", file.path(data_folder, "ss3_n4298_fibs_siKD_meta.rds"))
    }
    cast_raw <- readRDS(file.path(data_folder, "ss3_n4298_fibs_siKD_umiCast.rds"))
    c57_raw <- readRDS(file.path(data_folder, "ss3_n4298_fibs_siKD_umiC57.rds"))
    meta <- readRDS(file.path(data_folder, "ss3_n4298_fibs_siKD_meta.rds"))
    stopifnot(all(rownames(c57_raw) == rownames(cast_raw)))
    stopifnot(all(colnames(c57_raw) == colnames(cast_raw)))
    total_raw <- cast_raw + c57_raw
    sf <- colSums2(total_raw, na.rm = TRUE)
    sf <- sf / median(sf)
    total_raw <- total_raw[!matrixStats::rowAlls(is.na(total_raw)), ,drop=FALSE]
    is_imputed <- is.na(total_raw)
    # Impute missing values based on row means
    for(row in seq_len(nrow(total_raw))){
      missing_obs <- is.na(total_raw[row,])
      mus <- matrixStats::mean2(total_raw[row,], idxs = ! missing_obs) * sf[missing_obs]
      total_raw[row, missing_obs] <- rpois(length(mus), lambda = mus)
    }
    stopifnot(all(meta$Sample == colnames(total_raw)))
    sce <- SingleCellExperiment(list(counts = total_raw, is_imputed = is_imputed), colData = DataFrame(meta))
    saveRDS(sce, file.path(data_folder, "smart_seq3_siRNA_knockdown_full.RDS"))
  }else{
    sce <- readRDS(file.path(data_folder, "smart_seq3_siRNA_knockdown_full.RDS"))
  }
  sce
}




data_loaders <- list(mcSCRB = load_mcSCRBseq_data, 
                     smartSeq3_fibroblasts = load_smartSeq3_fibroblasts,
                     smartSeq3_fibroblasts_alt = load_smartSeq3_fibroblasts_alt,
                     smartSeq3_hek = load_smartSeq3_hek,
                     smartSeq3_siRNA_knockdown = load_smartSeq3_siRNA_knockdown_experiment)


