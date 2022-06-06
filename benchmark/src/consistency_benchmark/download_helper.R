library(tidyverse)
library(MatrixGenerics)
library(SingleCellExperiment)

data_folder <- "/scratch/ahlmanne/benchmark/data"

download_geo_suppl <- function(id, data_dir = data_folder, file_regex = NULL){
  vapply(id, function(.id){
    file_names <- GEOquery::getGEOSuppFiles(.id, baseDir = data_dir, filter_regex = file_regex, fetch_files = FALSE)$fname
    files_exist <- file.exists(file.path(data_dir, file_names))
    for(idx in seq_along(file_names)){
      if(! files_exist[idx]){
        GEOquery::getGEOSuppFiles(.id, baseDir = data_dir, filter_regex = file_names[idx])
      }
    }
    file.path(data_dir, .id)
  }, FUN.VALUE = "")
}

rename_files <- function(dir, pattern = c("matrix" = "^(.*matrix)\\.mtx(?:\\.gz)?$",
                                          "barcodes" = "^(.*barcodes)\\.tsv(?:\\.gz)?$",
                                          "genes" = "^(.*genes)\\.tsv(?:\\.gz)?$",
                                          "features" = "^(.*features)\\.tsv(?:\\.gz)?$")){
  stopifnot(! is.null(names(pattern)))
  for(.dir in dir){
    for(file in list.files(.dir)){
      matches <- str_match(file, pattern = pattern)
      stopifnot(ncol(matches) == 2)
      stopifnot(nrow(matches) == length(pattern))
      hits <- which(! is.na(matches[,2]))
      if(length(hits) >= 2){
        stop("Too many hits. ", file, " matched pattern at position ", paste0(hits, collapse = ", "))
      }
      if(length(hits) == 1){
        loc <- str_locate(file, matches[hits,2])
        new_file_name <- str_c(str_sub(file, start = 1L, end = loc[1]-1), names(pattern)[hits], str_sub(file, start = loc[2]+1, end = -1))
        file.rename(from = file.path(.dir, file), to = file.path(.dir, new_file_name))
      }
    }
  }
  dir
}


# Ashton et al.
# Homo sapiens
# Illumina HiSeq 2500
get_GSE142647_data <- function(){
  ids <- c("GSM4235299", "GSM4235300")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids) %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 



# Zemmour et al.
# Mus musculus
# Illumina NovaSeq 6000
get_GSE178765_data <- function(){
  ids <- c("GSE178765")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    dir <- download_geo_suppl(ids) %>%
      rename_files() 
    gene_info <- readr::read_tsv(file.path(dir, "genes.tsv.gz"), col_names = FALSE)
    gene_info$X2 <- gene_info$X1
    readr::write_tsv(gene_info, file.path(dir, "genes.tsv.gz"), col_names = FALSE)
    DropletUtils::read10xCounts(dir)
  }
} 


# Qian et al. "ZEB1 promotes pathogenic Th1 and Th17 cell differentiation in multiple sclerosis"
# Homo sapiens
# Illumina HiSeq 2500
get_GSE179831_data <- function(){
  ids <- c("GSM5434863", "GSM5434864")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids) %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 


# Pal et al. "Single cell transcriptome atlas of mouse mammary epithelial cells across development"
# Mus musculus
# Illumina NextSeq 500
get_GSE164017_data <- function(){
  ids <- "GSM4994960"
  if(! file.exists(file.path(data_folder, ids, "features.tsv.gz"))){
    if(! file.exists(file.path(data_folder, ids))){
      dir.create(file.path(data_folder, ids))
    }
    GEOquery::getGEOSuppFiles("GSE164017", baseDir = file.path(data_folder, ids), filter_regex = "features.tsv.gz", makeDirectory = FALSE)
    rename_files(file.path(data_folder, ids))
  }
  
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids) %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 


# Kathiriya et al.
# Homo sapiens
# Illumina NovaSeq 6000
get_GSE150068_data <- function(){
  ids <- c("GSM4522986")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids) %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 

# Bulaeva et al. "MYC-induced human acute myeloid leukemia requires a continuing IL-3/GM-CSF costimulus"
# Homo sapiens
# Illumina NextSeq 500
get_GSE130931_data <- function(){
  ids <- c("GSM4041124", "GSM4041125")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids) %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 





# Ding et al. "Targeting epigenetically maladapted vascular niche alleviates liver fibrosis in nonalcoholic steatohepatitis"
# Sus scrofa
# Illumina NextSeq 6000
get_GSE181483_data <- function(){
  ids <- c("GSM5503300", "GSM5503301", "GSM5503302")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids) %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 



# De Santis et al.
# Homo sapiens
# Illumina NextSeq 6000
get_GSE163505_data <- function(){
  ids <- c("GSM4980292")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids) %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 


# Morrow et al.
# Mus musculus
# Illumina HiSeq 4000
get_GSE158941_data <- function(){
  ids <- c("GSM4816083")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids) %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 


# Iatridis et al.
# Bos taurus
# Illumina NovaSeq 6000
get_GSE179714_data <- function(){
  ids <- c("GSM5429729", "GSM5429730")
  if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    download_geo_suppl(ids, file_regex = "filtered-bc") %>%
      rename_files() %>%
      DropletUtils::read10xCounts()
  }
} 


# Lu et al.
# Homo sapiens
# Illumina NovaSeq 6000
get_GSE184806_data <- function(){
  ids <- c("GSE184806")
  sce <- if(all(file.exists(file.path(data_folder, ids, "matrix.mtx.gz")))){
    DropletUtils::read10xCounts(file.path(data_folder, ids))
  }else{
    dir <- download_geo_suppl(ids, file_regex = "batch12") %>%
      rename_files() 
    rename_files(dir, pattern = c("metadata" = "^(.+metadata).tsv.gz"))
    
    DropletUtils::read10xCounts(dir)
  }
  meta <- readr::read_tsv(file.path(data_folder, ids, "metadata.tsv.gz"))
  colData(sce) <- cbind(colData(sce), as.data.frame(meta))
  sce
} 


data_loaders <- list(GSE142647 = get_GSE142647_data, GSE178765 = get_GSE178765_data, 
                     GSE179831 = get_GSE179831_data, GSE164017 = get_GSE164017_data, 
                     GSE150068 = get_GSE150068_data, 
                     GSE130931 = get_GSE130931_data, GSE181483 = get_GSE181483_data, 
                     GSE163505 = get_GSE163505_data, GSE158941 = get_GSE158941_data, 
                     GSE179714 = get_GSE179714_data, GSE184806 = get_GSE184806_data)





