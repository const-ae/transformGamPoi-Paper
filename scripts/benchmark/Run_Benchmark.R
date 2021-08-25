#!/g/easybuild/x86_64/CentOS/7/haswell/software/R/4.0.0-foss-2020a/bin/Rscript


library(tidyverse)

# argv <- c("--id", "uddxn", "--knn", 20, "--test_run", "--linear_random_tree")

parsed_args <- argparser::arg_parser("Specify the parameters for the overlap calculation") %>%
  argparser::add_argument("--id", type = "character", help = "A list of ids to calculate the overlap for") %>%
  argparser::add_argument("--linear_random_tree", flag = TRUE, help = "Flag to signal that instead of running on the 'random_tree/' folder, the 'linear_random_tree/' folder is used") %>%
  argparser::add_argument("--test_run", flag = TRUE, help = "Run on a significantly reduced set of genes") %>%
  argparser::add_argument("--knn", type = "character", help = "The number of nearest neighbors to overlap", default = "100") %>%
  argparser::add_argument("--pca_dimensions", type = "character", help = "The number of dimensions to use for the PCA", default = "50") %>%
  # argparser::parse_args(argv = argv)
  argparser::parse_args()

if(is.na(parsed_args$id)){
  stop("--id must not be NA")
}

setwd("/home/ahlmanne/projects/transformGamPoi-Paper/src")

script_location <- "./2021-07-05_Run_Benchmark.R"
script_hash <- read_lines(script_location) %>%
  paste0(collapse = "\n") %>%
  digest::digest() %>%
  str_sub(end = 6)

print(paste0("script_hash: ", script_hash))


random_tree_folder <- if(! parsed_args$linear_random_tree){
  "random_tree"
}else{
  "linear_random_tree"
}

input_folder <- file.path("../intermediate/benchmark/", random_tree_folder, parsed_args$id)
output_dir <- file.path("../output/benchmark/", random_tree_folder)

script_backup_location <- file.path(output_dir, "script_backup", paste0(script_hash, ".R"))

# Copy this script to output folder for backup and reference
if(! file.exists(script_backup_location)){
  file.copy(from = script_location,
            to = script_backup_location)
}



print(paste0("input_folder: ", input_folder))



# Read input files
delta_true <- readRDS(file.path(input_folder, "delta_true.RDS"))
UMI <- readRDS(file.path(input_folder, "UMI.RDS"))
stopifnot(nrow(delta_true) == nrow(UMI) && ncol(delta_true) == ncol(UMI))

if(parsed_args$test_run){
  n_cells <- 500
  n_genes <- 200
  sel_cells <- sample(which(matrixStats::colSums2(UMI) > 10), n_cells)
  sel_genes <- sample(which(matrixStats::rowSums2(UMI) > 0), n_genes)
  UMI <- UMI[sel_genes, sel_cells]
  delta_true <- delta_true[sel_genes, sel_cells]
}
expressed_cells <- matrixStats::colSums2(UMI) > 0
expressed_genes <- matrixStats::rowSums2(UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]
delta_true <- delta_true[expressed_genes, expressed_cells]


###########################################################################
#                APPLY TRANSFORMATIONS
###########################################################################

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





sf <- matrixStats::colSums2(UMI)
sf <- sf / mean(sf)
summary(sf)

logp1_timing <- system.time(
  logp1 <- transformGamPoi::shifted_log_transform(UMI, pseudo_count = 1, size_factors = sf)
)
print("Finished log transformation")

colsums <- matrixStats::colSums2(UMI)
logp_cpm_timing <- system.time(
  logp_cpm <- log1p(t(t(UMI) / colsums * 1e6))
)
print("Finished log cpm transformation")

alpha <- 0.05
logp_alpha_timing <- system.time(
  logp_alpha <- transformGamPoi::shifted_log_transform(UMI, overdispersion = alpha, size_factors = sf)
)
print("Finished log plus alpha transformation")

acosh_timing <- system.time(
  acosh <- transformGamPoi::acosh_transform(UMI, overdispersion = alpha, size_factors = sf)
)
print("Finished acosh")

pearson_timing <- system.time(
  pearson <- transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "pearson", size_factor = sf)
)
print("Finished pearson")

pearson_clip_timing <- system.time(
  pearson_clip <- transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "pearson", size_factor = sf, clipping = TRUE)
)
print("Finished pearson clipping")

rand_quantile_timing <- system.time(
  rand_quantile <- transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "randomized_quantile", size_factor = sf)
)
# This is necessary because to a bug in pnbinom which diverges for some unfortunate combination of input values
rand_quantile[matrixStats::rowAnyMissings(rand_quantile), ] <- 0
print("Finished randomized quantile")

colnames(UMI) <- paste0("Cell_", seq_len(ncol(UMI)))
rownames(UMI) <- paste0("Gene_", seq_len(nrow(UMI)))
sctransform_timing <- system.time(
  # Installed via `devtools::install_github("saketkc/sctransform@ff4422")`
  sctransform <- sctransform::vst(UMI, vst.flavor = "v2")$y
)

sanity_main_timing <- system.time(
  sanity_main <- run_sanity(UMI, n_threads = 10)
)
print("Finished sanity main")

sanity_dists_timing <- system.time(
  sanity_dists <- sanity_distance(sanity_main$outdir, n_threads = 10)
)
print("Finished sanity dists")


timing_df <- bind_rows(
  mutate(method = 'logp1', enframe(logp1_timing)),
  mutate(method = 'logp_cpm', enframe(logp_cpm_timing)),
  mutate(method = 'logp_alpha', enframe(logp_alpha_timing)),
  mutate(method = 'acosh', enframe(acosh_timing)),
  mutate(method = 'pearson', enframe(pearson_timing)),
  mutate(method = 'pearson_clip', enframe(pearson_clip_timing)),
  mutate(method = 'rand_quantile', enframe(rand_quantile_timing)),
  mutate(method = 'sctransform', enframe(sctransform_timing)),
  mutate(method = 'sanity_main', enframe(sanity_main_timing)),
  mutate(method = 'sanity_dists', enframe(sanity_dists_timing))
)


###########################################################################
#                Benchmark with ground truth
###########################################################################


cell_cell_distances <- list(sanity_dists = sanity_dists)
transformed_data <- list(logp1 = logp1, logp_cpm = logp_cpm, logp_alpha = logp_alpha, acosh = acosh,
                         pearson = pearson, pearson_clip = pearson_clip, rand_quantile = rand_quantile, 
                         sctransform = sctransform, sanity_main = sanity_main$mean)



knn_for_dist <- function(x, k){
  tmp_mat <- as.matrix(x)
  stopifnot(k < ncol(tmp_mat)-1)
  diag(tmp_mat) <- Inf
  proDA:::mply_dbl(seq_len(nrow(tmp_mat)), function(idx){
    order(tmp_mat[idx, ])[seq_len(k)]
  }, ncol = k)
}


for(k_nearest_neighbors in as.integer(str_split(parsed_args$knn, ",")[[1]])){
  print(paste0("k_nearest_neighbors: ", k_nearest_neighbors))
  system.time(
    ref_knn <- BiocNeighbors::findAnnoy(t(delta_true), k = k_nearest_neighbors, warn.ties = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))
  )
  
  annoy_nearest_neighbors <- bind_rows(lapply(transformed_data, function(dat){
    dur <- system.time({
      knn <- BiocNeighbors::findAnnoy(t(dat), k = k_nearest_neighbors, warn.ties = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))$index
      overlap <- sapply(seq_len(nrow(knn)), function(idx){
        sum(ref_knn$index[idx,] %in% knn[idx,])
      })
    })
    tibble(mean = mean(overlap), median = median(overlap), min = min(overlap), max = max(overlap),
           duration_elapsed = unname(dur["elapsed"]), 
           duration_system = unname(dur["sys.self"]), 
           duration_user = unname(dur["user.self"]),
           duration_sys_child = unname(dur["sys.child"]),
           duration_user_child = unname(dur["user.child"]))
  })) %>%
    mutate(method = names(transformed_data))
  
  distance_based_nearest_neighbors <- bind_rows(lapply(cell_cell_distances, function(dat){
    dur <- system.time({
      knn <- knn_for_dist(dat, k = k_nearest_neighbors)
      overlap <- sapply(seq_len(nrow(knn)), function(idx){
        sum(ref_knn$index[idx,] %in% knn[idx,])
      })
    })
    tibble(mean = mean(overlap), median = median(overlap), min = min(overlap), max = max(overlap),
           duration_elapsed = unname(dur["elapsed"]), 
           duration_system = unname(dur["sys.self"]), 
           duration_user = unname(dur["user.self"]),
           duration_sys_child = unname(dur["sys.child"]),
           duration_user_child = unname(dur["user.child"]))
  })) %>%
    mutate(method = names(cell_cell_distances))
  
  for(pca_dims in as.integer(str_split(parsed_args$pca_dimensions, ",")[[1]])){
    print(paste0("pca_dims: ", pca_dims))
    
    output_file <- file.path(output_dir, 
                             paste0(random_tree_folder, "-", if(parsed_args$test_run) "test_run-" else "full_run-",
                                    "id_", parsed_args$id, "-knn_", k_nearest_neighbors, "-pcadim_", pca_dims,
                                    "-script_", script_hash, "-results.tsv"))
    print(paste0("output_file: ", output_file))
    
    
    
    pca_nearest_neighbors <- bind_rows(lapply(transformed_data, function(dat){
      dur <- system.time({
        pca <- irlba::prcomp_irlba(t(dat), n = pca_dims)
        knn <- BiocNeighbors::findAnnoy(pca$x, k = k_nearest_neighbors, warn.ties = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 10))$index
        overlap <- sapply(seq_len(nrow(knn)), function(idx){
          sum(ref_knn$index[idx,] %in% knn[idx,])
        })
      })
      tibble(mean = mean(overlap), median = median(overlap), min = min(overlap), max = max(overlap),
             duration_elapsed = unname(dur["elapsed"]), 
             duration_system = unname(dur["sys.self"]), 
             duration_user = unname(dur["user.self"]),
             duration_sys_child = unname(dur["sys.child"]),
             duration_user_child = unname(dur["user.child"]))
    })) %>%
      mutate(method = paste0("pca_", names(transformed_data)))
    
    
    
    res <- bind_rows(annoy_nearest_neighbors, distance_based_nearest_neighbors, pca_nearest_neighbors) %>%
      mutate(id = parsed_args$id,
             script_hash = script_hash,
             knn = k_nearest_neighbors,
             pcadims = pca_dims) %>%
      mutate(preproc_method = str_remove(method, "pca_")) %>%
      left_join({
        timing_df %>%
          pivot_wider(method, names_from = name, values_from = value) %>%
          dplyr::select(method, preproc_duration_elapsed = elapsed,
                        preproc_duration_system = sys.self, 
                        preproc_duration_user = user.self,
                        preproc_duration_sys_child = sys.child,
                        preproc_duration_user_child = user.child)
      }, by = c("preproc_method" = "method"))
    
    
    write_tsv(res, output_file)
    print(paste0("Finished knn=", k_nearest_neighbors, " and pcadim=", pca_dims))
  }
}

print("Finished looping")

sessionInfo()

print("Finished")
