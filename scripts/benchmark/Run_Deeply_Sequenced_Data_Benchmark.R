


library(tidyverse)
library(SingleCellExperiment)
library(MatrixGenerics)

setwd("/g/huber/users/ahlmanne/projects/transformGamPoi-Paper/src")




script_location <- "./2021-08-11_DSSC_Benchmark_script.R"
script_hash <- read_lines(script_location) %>%
  paste0(collapse = "\n") %>%
  digest::digest() %>%
  str_sub(end = 6)

print(paste0("script_hash: ", script_hash))

output_dir <- file.path("../output/benchmark/deeply_sequence_single_cell_benchmark/")
script_backup_location <- file.path(output_dir, "script_backup", paste0(script_hash, ".R"))

# Copy this script to output folder for backup and reference
if(! file.exists(script_backup_location)){
  file.copy(from = script_location,
            to = script_backup_location)
}



load_msSCRBseq_data <- function(){
  if(! file.exists("../intermediate/deep_single_cell_data/GSE103568_JM8_UMIcounts.txt.gz")){
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103568&format=file&file=GSE103568%5FJM8%5FUMIcounts%2Etxt%2Egz", "../intermediate/deep_single_cell_data/GSE103568_JM8_UMIcounts.txt.gz")
  }
  raw_data <- as.matrix(read.delim("../intermediate/deep_single_cell_data/GSE103568_JM8_UMIcounts.txt.gz"))
  SingleCellExperiment(list(counts = raw_data), rowData = DataFrame(gene_id = rownames(raw_data)))
}

load_smartSeq3_fibroblasts <- function(){
  if(! file.exists("../intermediate/deep_single_cell_data/Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt")){
    if(! file.exists("../intermediate/deep_single_cell_data/smart_seq3_data.zip")){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/E-MTAB-8735.processed.3.zip", "../intermediate/deep_single_cell_data/smart_seq3_data.zip")
    }
    unzip("../intermediate/deep_single_cell_data/smart_seq3_data.zip", files = "Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt", exdir = "../intermediate/deep_single_cell_data")
  }
  if(! file.exists("../intermediate/deep_single_cell_data/Smartseq3.Fibroblasts.sample_annotation.txt")){
    if(! file.exists("../intermediate/deep_single_cell_data/smart_seq3_data_annotation.zip")){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/E-MTAB-8735.processed.2.zip", "../intermediate/deep_single_cell_data/smart_seq3_data_annotation.zip")
    }
    unzip("../intermediate/deep_single_cell_data/smart_seq3_data_annotation.zip", files = "Smartseq3.Fibroblasts.sample_annotation.txt", exdir = "../intermediate/deep_single_cell_data")
  }
  
  col_data <- read.delim("../intermediate/deep_single_cell_data/Smartseq3.Fibroblasts.sample_annotation.txt")
  raw_data <- read.delim("../intermediate/deep_single_cell_data/Smartseq3.Fibroblasts.NovaSeq.UMIcounts.txt") %>%
    as.matrix() %>%
    as("dgCMatrix")
  
  rownames(col_data) <- col_data$BC
  col_data <- col_data[colnames(raw_data), ]
  stopifnot(all(colnames(raw_data) == col_data$BC))
  
  
  SingleCellExperiment(list(counts = raw_data), rowData = DataFrame(gene_id = rownames(raw_data)), 
                       colData = DataFrame(col_data))
}

load_smartSeq3_fibroblasts_alt <- function(){
  if(! file.exists("../intermediate/deep_single_cell_data/Fibroblasts.plate2.umis.ex.txt")){
    if(! file.exists("../intermediate/deep_single_cell_data/E-MTAB-10148.processed.1.zip")){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-10148/E-MTAB-10148.processed.1.zip", "../intermediate/deep_single_cell_data/E-MTAB-10148.processed.1.zip")
    }
    if(! file.exists("../intermediate/deep_single_cell_data/E-MTAB-10148.processed.2.zip")){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-10148/E-MTAB-10148.processed.2.zip", "../intermediate/deep_single_cell_data/E-MTAB-10148.processed.2.zip")
    }
    unzip("../intermediate/deep_single_cell_data/E-MTAB-10148.processed.1.zip", files = "Fibroblasts.plate1.umis.ex.txt", exdir = "../intermediate/deep_single_cell_data")
    unzip("../intermediate/deep_single_cell_data/E-MTAB-10148.processed.1.zip", files = "Fibroblasts.plate2.umis.ex.txt", exdir = "../intermediate/deep_single_cell_data")
    unzip("../intermediate/deep_single_cell_data/E-MTAB-10148.processed.2.zip", files = "Fibroblasts.plate1.annotation.txt", exdir = "../intermediate/deep_single_cell_data")
    unzip("../intermediate/deep_single_cell_data/E-MTAB-10148.processed.2.zip", files = "Fibroblasts.plate2.annotation.txt", exdir = "../intermediate/deep_single_cell_data")
  }
  
  raw_data_1 <- read.delim("../intermediate/deep_single_cell_data/Fibroblasts.plate1.umis.ex.txt") %>%
    as.matrix() %>%
    as("dgCMatrix")
  raw_data_2 <- read.delim("../intermediate/deep_single_cell_data/Fibroblasts.plate2.umis.ex.txt") %>%
    as.matrix() %>%
    as("dgCMatrix")
  annot_1 <- read.delim("../intermediate/deep_single_cell_data/Fibroblasts.plate1.annotation.txt")
  annot_2 <- read.delim("../intermediate/deep_single_cell_data/Fibroblasts.plate2.annotation.txt")
  
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
  if(! file.exists("../intermediate/deep_single_cell_data/Smartseq3.HEK.cleanup.UMIcounts.txt")){
    if(! file.exists("../intermediate/deep_single_cell_data/smart_seq3_data.zip")){
      download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8735/E-MTAB-8735.processed.3.zip", "../intermediate/deep_single_cell_data/smart_seq3_data.zip")
    }
    unzip("../intermediate/deep_single_cell_data/smart_seq3_data.zip", files = "Smartseq3.HEK.cleanup.UMIcounts.txt", exdir = "../intermediate/deep_single_cell_data")
  }
  raw_data <- read.delim("../intermediate/deep_single_cell_data/Smartseq3.HEK.cleanup.UMIcounts.txt") %>%
    as.matrix() %>%
    as("dgCMatrix")
  
  SingleCellExperiment(list(counts = raw_data), rowData = DataFrame(gene_id = rownames(raw_data)))
}




sce_list <- list(smartSeq3_fibro = load_smartSeq3_fibroblasts(),
                 smartSeq3_fibro_alt = load_smartSeq3_fibroblasts_alt(),
                 msSCRB = load_msSCRBseq_data(),
                 smartSeq3_hek = load_smartSeq3_hek())


# Run Sanity
run_sanity <- function(x, variance_min = 0.01, variance_max = 20, n_bins = 116, n_threads = 1){
  stopifnot(is.matrix(x))
  stopifnot(all(matrixStats::colSums2(x) > 0))
  stopifnot(all(matrixStats::rowSums2(x) > 0))
  # write x to a file
  input <- tempfile(pattern = "matrix", fileext = ".tsv")
  output_dir <- file.path("/scratch/ahlmanne/Sanity_output/", 
                          paste0("sanity_output_", paste0(sample(letters, 7, replace = TRUE), collapse = "")))
  dir.create(output_dir)
  print(paste0("Output dir: ", output_dir))
  
  if(is.null(rownames(x))){
    rownames(x) <- paste0("Gene_", seq_len(nrow(x)))
  }
  if(is.null(colnames(x))){
    colnames(x) <- paste0("Cell_", seq_len(ncol(x)))
  }
  readr::write_tsv(tibble::as_tibble(x, rownames = "GeneID"), input)
  
  code <- system2("/home/ahlmanne/prog/Sanity/bin/Sanity", args = c(paste0("--file ", input),
                                                                    paste0("--destination ", output_dir),
                                                                    paste0("--variance_min ", variance_min),
                                                                    paste0("--variance_max ", variance_max),
                                                                    paste0("--number_of_bins ", n_bins),
                                                                    paste0("--n_threads ", n_threads),
                                                                    "--extended_output true"))
  if(code != 0){
    stop("Sanity failed with status ", code)
  }
  
  
  res <- readr::read_tsv(file.path(output_dir, "log_transcription_quotients.txt"), col_types = readr::cols()) 
  mat <- as.matrix(res[,-1,drop=FALSE])
  rownames(mat) <- res[[1]]
  colnames(mat) <- colnames(res)[-1]
  
  errror_bars_df <- readr::read_tsv(file.path(output_dir, "ltq_error_bars.txt"), col_types = readr::cols()) 
  error_bars <- as.matrix(errror_bars_df[,-1,drop=FALSE])
  rownames(error_bars) <- errror_bars_df[[1]]
  colnames(error_bars) <- colnames(errror_bars_df)[-1]
  
  
  list(mean = mat, sd = error_bars, outdir = output_dir)
}

sanity_distance <- function(directory, take_error_into_account = TRUE, signal_to_noise_cutoff = 1, n_threads = 1){
  if(signal_to_noise_cutoff < 0.000001){
    stop("signal_to_noise_cutoff is too small")
  }
  code <- system2("/home/ahlmanne/prog/Sanity/bin/Sanity_distance", args = c(paste0("--folder ", directory),
                                                                             paste0("--with_error_bars ", if(take_error_into_account) "true" else "false"),
                                                                             paste0("--signal_to_noise_cutoff ", signal_to_noise_cutoff),
                                                                             paste0("--n_threads ", n_threads)))
  if(code != 0){
    stop("Sanity failed with status ", code)
  }
  
  s2n_string <- stringr::str_remove(sprintf("%f", signal_to_noise_cutoff), "[0.]*$")
  
  res_file <- file.path(directory, paste0("cell_cell_distance_", if(take_error_into_account) "with_errorbar" else "euclidean", "_s2n_gt_", s2n_string, ".txt"))
  
  dist_vec <- readr::read_tsv(res_file, col_names = "dist", col_types = readr::cols())[[1]] 
  n_cells <- 1/2 + sqrt(1/4 + length(dist_vec) * 2)
  mat <- matrix(NA, nrow = n_cells, ncol = n_cells)
  mat[lower.tri(mat)] <- dist_vec
  as.dist(mat)
}

knn_for_dist <- function(x, k){
  tmp_mat <- as.matrix(x)
  stopifnot(k < ncol(tmp_mat)-1)
  diag(tmp_mat) <- Inf
  proDA:::mply_dbl(seq_len(nrow(tmp_mat)), function(idx){
    order(tmp_mat[idx, ])[seq_len(k)]
  }, ncol = k)
}


pca_dims <- 15
alphas <- list(msSCRB = 0.1, smartSeq3_hek = 0.1,
               smartSeq3_fibro = TRUE, smartSeq3_fibro_alt = TRUE)
k_nearest_neighbors <- 30
n_iter <- 10
thresholds <- list(msSCRB = c(1e4, 1e5),
                   smartSeq3_fibro = c(-Inf, Inf),
                   smartSeq3_fibro_alt = c(8e4, Inf),
                   smartSeq3_hek = c(1e4, Inf))


result <- tibble(dataset = character(0),
                 iteration = integer(0L),
                 overlap = list(),
                 common_knns = list())

knns_per_dataset <- list()

for(dataset_name in names(sce_list)){
  full_sce <- sce_list[[dataset_name]]
  alpha <- alphas[[dataset_name]]
  colsums <- colSums2(counts(full_sce))
  thres <- thresholds[[dataset_name]]
  sce <- full_sce[, colsums > thres[1] & colsums < thres[2]]
  sce <- sce[rowSums2(counts(sce)) > 0, ]
  print(paste0("Dataset ", dataset_name, ". Cells ", ncol(sce), " and genes ", nrow(sce)))
  
  colsums <- colSums2(counts(sce))
  sf <- colsums / mean(colsums)
  downsample_proportion <- 5000 / mean(colsums)
  
  cat(paste0("Overlap for ", dataset_name, " with ", k_nearest_neighbors, " knns. PCA=", pca_dims, " and alpha=", alpha, "\n"))
  
  assay(sce, "logp1")  <- transformGamPoi::shifted_log_transform(counts(sce), pseudo_count = 1, size_factors = sf)
  assay(sce, "logp_cpm")  <- log1p(t(t(counts(sce)) / colsums * 1e6))
  assay(sce, "logp_alpha")  <- transformGamPoi::shifted_log_transform(counts(sce), overdispersion = alpha, size_factors = sf, on_disk = FALSE)
  assay(sce, "acosh")  <- transformGamPoi::acosh_transform(counts(sce), overdispersion = alpha, size_factors = sf, on_disk = FALSE)
  assay(sce, "pearson")  <- transformGamPoi::residual_transform(counts(sce), residual_type = "pearson",
                                                                size_factor = sf, overdispersion = alpha, on_disk = FALSE)
  assay(sce, "pearson_clipping")  <- transformGamPoi::residual_transform(counts(sce), residual_type = "pearson", 
                                                                         size_factor = sf, overdispersion = alpha, clipping = TRUE, on_disk = FALSE)
  assay(sce, "rand_quantile")  <- transformGamPoi::residual_transform(counts(sce), residual_type = "randomized_quantile",
                                                              size_factor = sf, overdispersion = alpha, on_disk = FALSE)
  sanity_main <- run_sanity(as.matrix(counts(sce)), n_threads = 10)
  assay(sce, "sanity_vst") <- sanity_main$mean

  # KNNs for each transformation (except 'counts')
  knns <- lapply(assays(sce)[-1], function(as){
    pca <- irlba::prcomp_irlba(t(as), n = pca_dims)
    BiocNeighbors::findAnnoy(pca$x, k = k_nearest_neighbors, warn.ties = FALSE, 
                             BPPARAM = BiocParallel::MulticoreParam(workers = 10))$index
  })
  
  san_dist <- sanity_distance(sanity_main$outdir, n_threads = 10)
  san_dist_knn <- knn_for_dist(san_dist, k_nearest_neighbors)
  knns <- c(knns, list(sanity_dist = san_dist_knn))
  
  knns_per_dataset <- c(knns_per_dataset, list(knns))
  
  # Nearest neighbors found with all transformations
  colData(sce)$common_knns <- lapply(seq_len(ncol(sce)), function(idx){
    merged_nn <- lapply(knns, function(knn) knn[idx, ])
    colnames(sce)[purrr::reduce(merged_nn, intersect)]
  })
  print(summary(lengths(sce$common_knns) / k_nearest_neighbors))
  
  for(iter in seq_len(n_iter)){
    # Downsample Data       
    ds_sce <- SingleCellExperiment(list(
      counts = scuttle::downsampleMatrix(counts(sce), prop = downsample_proportion, bycol = FALSE)
    ), colData = colData(sce))
    ds_sce <- ds_sce[rowSums2(counts(ds_sce)) > 0, colSums2(counts(ds_sce)) > 0]
    ds_colsums <- colSums2(counts(ds_sce))
    ds_sf <- ds_colsums / mean(ds_colsums)
    
    # Calculate Transformations on Downsampled Data
    assay(ds_sce, "logp1")  <- transformGamPoi::shifted_log_transform(counts(ds_sce), pseudo_count = 1, size_factors = ds_sf)
    assay(ds_sce, "logp_cpm")  <- log1p(t(t(counts(ds_sce)) / ds_colsums * 1e6))
    assay(ds_sce, "logp_alpha")  <- transformGamPoi::shifted_log_transform(counts(ds_sce), overdispersion = alpha, size_factors = ds_sf, on_disk = FALSE)
    assay(ds_sce, "acosh")  <- transformGamPoi::acosh_transform(counts(ds_sce), overdispersion = alpha, size_factors = ds_sf, on_disk = FALSE)
    assay(ds_sce, "pearson")  <- transformGamPoi::residual_transform(counts(ds_sce), residual_type = "pearson", 
                                                                     size_factor = ds_sf, overdispersion = alpha, on_disk = FALSE)
    assay(ds_sce, "pearson_clipping")  <- transformGamPoi::residual_transform(counts(ds_sce), residual_type = "pearson", 
                                                                              size_factor = ds_sf, overdispersion = alpha, clipping = TRUE, on_disk = FALSE)
    assay(ds_sce, "rand_quantile")  <- transformGamPoi::residual_transform(counts(ds_sce), residual_type = "randomized_quantile",
                                                                           size_factor = ds_sf, overdispersion = alpha, on_disk = FALSE)
    sanity_main_downsampled <- run_sanity(as.matrix(counts(ds_sce)), n_threads = 10)
    assay(ds_sce, "sanity_vst") <- sanity_main_downsampled$mean
    
    # Calculate KNNs after downsampling and applying PCA 
    knns_after_downsampling <- lapply(assays(ds_sce)[-1], function(as){
      pca <- irlba::prcomp_irlba(t(as), n = pca_dims)
      knn_indices <- BiocNeighbors::findAnnoy(pca$x, k = k_nearest_neighbors, warn.ties = FALSE, 
                                              BPPARAM = BiocParallel::MulticoreParam(workers = 10))$index
      matrix(colnames(ds_sce)[knn_indices], nrow = ncol(ds_sce), ncol = k_nearest_neighbors)
    })
    san_dist <- sanity_distance(sanity_main_downsampled$outdir, n_threads = 10)
    san_dist_knn <- knn_for_dist(san_dist, k_nearest_neighbors)
    knns_after_downsampling <- c(knns_after_downsampling,
                                 list(sanity_dist = matrix(colnames(ds_sce)[san_dist_knn], nrow = ncol(ds_sce), ncol = k_nearest_neighbors)))
    
    
    overlapping_knns_per_trans <- lapply(knns_after_downsampling, function(knn){
      stopifnot(nrow(knn) == length(ds_sce$common_knns))
      overlap <- sapply(seq_len(ncol(ds_sce)), function(idx){
        com_knn <- ds_sce$common_knns[[idx]]
        sum(knn[idx, ] %in% com_knn) 
      })
    })
    df <- enframe(round(map_dbl(overlapping_knns_per_trans, function(x) mean(x/lengths(ds_sce$common_knns), na.rm = TRUE)), 2))
    cat(paste0(format(df), collapse = "\n"))
    cat("\n")
    
    result <- bind_rows(result,
                        tibble(dataset = dataset_name,
                               iteration = iter,
                               overlap = list(overlapping_knns_per_trans),
                               common_knns = list(ds_sce$common_knns)))
  }
}
names(knns_per_dataset) <- names(sce_list)

saveRDS(result, file.path(output_dir, paste0("deep_seq-result-", script_hash, "-overlap.Rds")))
saveRDS(knns_per_dataset, file.path(output_dir, paste0("deep_seq-result-", script_hash, "-raw_knns.Rds")))

# Print software versions
sessionInfo()
