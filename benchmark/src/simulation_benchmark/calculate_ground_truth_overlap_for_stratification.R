library(tidyverse)

pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--knn_ids", type = "character", nargs = Inf, help = "The id of a file in output/results")
pa <- argparser::add_argument(pa, "--transformation_names", type = "character", nargs = Inf, help = "A readable identifier of the transformation")
pa <- argparser::add_argument(pa, "--simulator", type = "character", help = "[Just for documentation purposes] A readable identifier of the data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "[Just for documentation purposes] The seed used to tame randomness")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "[Just for documentation purposes] The number of PCA dimensions used before KNN graph construction")
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "[Just for documentation purposes] The number of nearest neighbors considered")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)
stopifnot(length(pa$knn_ids) == length(pa$transformation_names))

ground_truth <- readRDS(file.path(pa$working_dir, "results", pa$data_id))$ground_truth
groung_truth_knn <- BiocNeighbors::findAnnoy(t(ground_truth), k = pa$knn, warn.ties = FALSE, get.distance = FALSE)$index

overlaps <- lapply(pa$knn_ids, function(id){
  knn <- readRDS(file.path(pa$working_dir, "results", id))
  stopifnot(ncol(ground_truth) == nrow(knn))
  knn_overlap <- sapply(seq_len(nrow(knn)), function(cell_idx){
    length(intersect(groung_truth_knn[cell_idx,], knn[cell_idx,]))
  })
  knn_overlap
})
names(overlaps) <- pa$transformation_names

res <- bind_cols(tibble(ground_truth_id = pa$data_id, simulator = pa$simulator, 
                        seed = pa$seed, pca_dim = pa$pca_dim, knn = pa$knn),
                 overlaps)
write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))
