library(tidyverse)

pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--input_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset", type = "character", help = "[Just for documentation purposes] A readable identifier of the data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "[Just for documentation purposes] The seed used to tame randomness")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "[Just for documentation purposes] The number of PCA dimensions used before KNN graph construction")
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "[Just for documentation purposes] The number of nearest neighbors considered")
pa <- argparser::add_argument(pa, "--transformation", type = "character", help = "[Just for documentation purposes] A readable identifier of the transformation")
pa <- argparser::add_argument(pa, "--alpha", type = "character", help = "[Just for documentation purposes] Specification of the overdispersion.")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)


KNNs <- readRDS(file.path(pa$working_dir, "results", pa$input_id))
stopifnot(all(dim(KNNs[[1]]) == dim(KNNs[[2]])))
n_cells <- nrow(KNNs[[1]])


cons <- mean(sapply(seq_len(n_cells), function(cell_idx){
  length(intersect(KNNs[[1]][cell_idx,], KNNs[[2]][cell_idx,]))
}))

res <- tibble(mean_overlap = cons, transformation_id = pa$input_id, dataset = pa$dataset, seed = pa$seed,
              pca_dim = pa$pca_dim, knn = pa$knn, transformation = pa$transformation, alpha = pa$alpha)

write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))