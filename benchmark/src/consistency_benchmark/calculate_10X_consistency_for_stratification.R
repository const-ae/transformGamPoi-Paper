library(tidyverse)

pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--dataset_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--input_ids", type = "character", nargs = Inf, help = "The ids of a file in output/results") 
pa <- argparser::add_argument(pa, "--transformation_names", type = "character", nargs = Inf, help = "A readable identifier of the transformation")
pa <- argparser::add_argument(pa, "--dataset_name", type = "character", help = "[Just for documentation purposes] A readable identifier of the data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "[Just for documentation purposes] The seed used to tame randomness")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "[Just for documentation purposes] The number of PCA dimensions used before KNN graph construction")
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "[Just for documentation purposes] The number of nearest neighbors considered")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)
stopifnot(length(pa$input_ids) == length(pa$transformation_names))
n_transformations <- length(pa$transformation_names)

overlaps <- lapply(seq_len(n_transformations), function(idx){
  KNNs <- readRDS(file.path(pa$working_dir, "results", pa$input_ids[idx]))
  stopifnot(all(dim(KNNs[[1]]) == dim(KNNs[[2]])))
  n_cells <- nrow(KNNs[[1]])
  
  cons <- sapply(seq_len(n_cells), function(cell_idx){
    length(intersect(KNNs[[1]][cell_idx,], KNNs[[2]][cell_idx,]))
  })
})
names(overlaps) <- pa$transformation_names

res <- bind_cols(tibble(dataset = pa$dataset_name, seed = pa$seed, pca_dim = pa$pca_dim, knn = pa$knn),
                 overlaps)

write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))


