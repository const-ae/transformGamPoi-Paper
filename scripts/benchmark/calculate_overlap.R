#!/g/easybuild/x86_64/CentOS/7/haswell/software/R/4.0.0-foss-2020a/bin/Rscript


library(tidyverse)
set.seed(1)

# argv <- c("--id", "gbagq", "--knn", 20)

parsed_args <- argparser::arg_parser("Specify the parameters for the overlap calculation") %>%
  argparser::add_argument("--id", type = "character", help = "A list of ids to calculate the overlap for") %>%
  argparser::add_argument("--simpler_random_tree", flag = TRUE, help = "Flag to signal that instead of running on the 'random_tree/' folder, the 'simpler_random_tree/' folder is used") %>%
  argparser::add_argument("--linear_random_tree", flag = TRUE, help = "Flag to signal that instead of running on the 'random_tree/' folder, the 'linear_random_tree/' folder is used") %>%
  argparser::add_argument("--test_run", flag = TRUE, help = "Run on a significantly reduced set of genes") %>%
  argparser::add_argument("--knn", type = "integer", help = "The number of nearest neighbors to overlap", default = 100) %>%
  argparser::add_argument("--pca_dimensions", type = "integer", help = "The number of dimensions to use for the PCA", default = 50) %>%
  # argparser::parse_args(argv = argv)
  argparser::parse_args()

if(is.na(parsed_args$id)){
  stop("--id must not be NA")
}

setwd("/home/ahlmanne/projects/transformGamPoi-Paper/src")

script_location <- "./2021-05-18_Benchmark_Random_Walk_Tree_overlap.R"
script_hash <- read_lines(script_location) %>%
  paste0(collapse = "\n") %>%
  digest::digest() %>%
  str_sub(end = 6)

random_tree_folder <- if(! parsed_args$simpler_random_tree && ! parsed_args$linear_random_tree){
  "random_tree"
}else if(parsed_args$simpler_random_tree){
  "simpler_random_tree"
}else{
  "linear_random_tree"
}

input_folder <- file.path("../intermediate/benchmark/", random_tree_folder, parsed_args$id)
output_dir <- file.path("../output/benchmark/", random_tree_folder)
output_file <- file.path(output_dir, 
                         paste0(random_tree_folder, "-", if(parsed_args$test_run) "test_run-" else "full_run-",
                                "id_", parsed_args$id, "-knn_", parsed_args$knn, "-pcadim_", parsed_args$pca_dimensions,
                                "-script_", script_hash, "-results.tsv"))
script_backup_location <- file.path(output_dir, "script_backup", paste0(script_hash, ".R"))

# Copy this script to output folder for backup and reference
if(! file.exists(script_backup_location)){
  file.copy(from = script_location,
            to = script_backup_location)
}



print(paste0("input_folder: ", input_folder))
print(paste0("output_file: ", output_file))

if(parsed_args$test_run){
  n_cells <- 500
  n_genes <- 200
  print(paste0("Subset to ", n_genes, " genes and ", n_cells, " cells"))
  delta_true <- readRDS(file.path(input_folder, "delta_true.RDS"))[seq_len(n_genes), seq_len(n_cells)]
  sanity_dists <- as.dist(as.matrix(readRDS(file.path(input_folder, "sanity_dists.RDS")))[seq_len(n_cells), seq_len(n_cells)])
  parents <- readRDS(file.path(input_folder, "parents.RDS"))
  logp1 <- readRDS(file.path(input_folder, "logp1.RDS"))[seq_len(n_genes), seq_len(n_cells)]
  logp_alpha <- readRDS(file.path(input_folder, "logp_alpha.RDS"))[seq_len(n_genes), seq_len(n_cells)]
  pearson <- readRDS(file.path(input_folder, "pearson.RDS"))[seq_len(n_genes), seq_len(n_cells)]
  rand_quantile <- readRDS(file.path(input_folder, "rand_quantile.RDS"))[seq_len(n_genes), seq_len(n_cells)]
  acosh <- readRDS(file.path(input_folder, "acosh.RDS"))[seq_len(n_genes), seq_len(n_cells)]
}else{
  delta_true <- readRDS(file.path(input_folder, "delta_true.RDS"))
  sanity_dists <- readRDS(file.path(input_folder, "sanity_dists.RDS"))
  parents <- readRDS(file.path(input_folder, "parents.RDS"))
  logp1 <- readRDS(file.path(input_folder, "logp1.RDS"))
  logp_alpha <- readRDS(file.path(input_folder, "logp_alpha.RDS"))
  pearson <- readRDS(file.path(input_folder, "pearson.RDS"))
  rand_quantile <- readRDS(file.path(input_folder, "rand_quantile.RDS"))
  acosh <- readRDS(file.path(input_folder, "acosh.RDS"))
}


cell_cell_distances <- list(sanity_dists = sanity_dists)
transformed_data <- list(logp1 = logp1, logp_alpha = logp_alpha, pearson = pearson, rand_quantile = rand_quantile, acosh = acosh)



knn_for_dist <- function(x, k){
  tmp_mat <- as.matrix(x)
  stopifnot(k < ncol(tmp_mat)-1)
  diag(tmp_mat) <- Inf
  proDA:::mply_dbl(seq_len(nrow(tmp_mat)), function(idx){
    order(tmp_mat[idx, ])[seq_len(k)]
  }, ncol = k)
}


k_nearest_neighbors <- parsed_args$knn
print(paste0("k_nearest_neighbors: ", k_nearest_neighbors))
pca_dims <- parsed_args$pca_dimensions
print(paste0("pca_dims: ", pca_dims))

system.time(
  ref_knn <- BiocNeighbors::findAnnoy(t(delta_true), k = k_nearest_neighbors, warn.ties = FALSE)
)

annoy_nearest_neighbors <- bind_rows(lapply(transformed_data, function(dat){
  dur <- system.time({
    knn <- BiocNeighbors::findAnnoy(t(dat), k = k_nearest_neighbors, warn.ties = FALSE)$index
    overlap <- sapply(seq_len(nrow(knn)), function(idx){
      sum(ref_knn$index[idx,] %in% knn[idx,])
    })
  })
  tibble(mean = mean(overlap), median = median(overlap), min = min(overlap), max = max(overlap),
         duration = unname(dur["elapsed"]))
})) %>%
  mutate(method = names(transformed_data))

distance_based_nearest_neighbors <- bind_rows(lapply(cell_cell_distances, function(dat){
  dur <- system.time({
    knn <- knn_for_dist(sanity_dists, k = k_nearest_neighbors)
    overlap <- sapply(seq_len(nrow(knn)), function(idx){
      sum(ref_knn$index[idx,] %in% knn[idx,])
    })
  })
  tibble(mean = mean(overlap), median = median(overlap), min = min(overlap), max = max(overlap),
         duration = unname(dur["elapsed"]))
})) %>%
  mutate(method = names(cell_cell_distances))


pca_nearest_neighbors <- bind_rows(lapply(transformed_data, function(dat){
  dur <- system.time({
    pca <- irlba::prcomp_irlba(t(dat), n = pca_dims)
    knn <- BiocNeighbors::findAnnoy(pca$x, k = k_nearest_neighbors, warn.ties = FALSE)$index
    overlap <- sapply(seq_len(nrow(knn)), function(idx){
      sum(ref_knn$index[idx,] %in% knn[idx,])
    })
  })
  tibble(mean = mean(overlap), median = median(overlap), min = min(overlap), max = max(overlap),
         duration = unname(dur["elapsed"]))
})) %>%
  mutate(method = paste0("pca_", names(transformed_data)))



res <- bind_rows(annoy_nearest_neighbors, distance_based_nearest_neighbors, pca_nearest_neighbors) %>%
  mutate(id = parsed_args$id,
         script_hash = script_hash)


write_tsv(res, output_file)
