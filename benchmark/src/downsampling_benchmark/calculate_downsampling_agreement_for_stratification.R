library(tidyverse)

pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--full_knn_result_ids", type = "character", nargs = Inf, help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--reduced_knn_result_ids", type = "character", nargs = Inf, help = "The id of a file in output/results")
pa <- argparser::add_argument(pa, "--dataset", type = "character", help = "[Just for documentation purposes] A readable identifier of the data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "[Just for documentation purposes] The seed used to tame randomness")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "[Just for documentation purposes] The number of PCA dimensions used before KNN graph construction")
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "[Just for documentation purposes] The number of nearest neighbors considered")
pa <- argparser::add_argument(pa, "--transformations", type = "character", nargs = Inf, help = "[Just for documentation purposes] A readable identifier of the transformation")
pa <- argparser::add_argument(pa, "--alphas", type = "character", nargs = Inf, help = "[Just for documentation purposes] Specification of the overdispersion.")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)
stopifnot(length(pa$transformations) == length(pa$full_knn_result_ids))
stopifnot(length(pa$transformations) == length(pa$reduced_knn_result_ids))
stopifnot(length(pa$alphas) == length(pa$full_knn_result_ids))
stopifnot(length(pa$alphas) == length(pa$reduced_knn_result_ids))


full_KNNs <- lapply(pa$full_knn_result_ids, function(id) readRDS(file.path(pa$working_dir, "results", id)))
# Filter out negative controls
full_KNNs <- full_KNNs[! pa$transformations %in% c("raw_counts", "scaled_raw_counts")]

reduced_KNNs <- lapply(pa$reduced_knn_result_ids, function(id) readRDS(file.path(pa$working_dir, "results", id)))

stopifnot(nrow(full_KNNs[[1]]) == sapply(full_KNNs, nrow))
stopifnot(ncol(full_KNNs[[1]]) == sapply(full_KNNs, ncol))
stopifnot(nrow(full_KNNs[[1]]) == sapply(reduced_KNNs, nrow))
stopifnot(ncol(full_KNNs[[1]]) == sapply(reduced_KNNs, ncol))

n_cells <- nrow(full_KNNs[[1]])

# Nearest neighbors found with all transformations
common_knns <- lapply(seq_len(n_cells), function(idx){
  merged_nn <- lapply(full_KNNs, function(knn) knn[idx, ])
  purrr::reduce(merged_nn, intersect)
})
n_common_knns <- lengths(common_knns)

overlaps <- lapply(reduced_KNNs, function(knn){
  overlap <- sapply(seq_len(n_cells), function(idx){
    com_knn <- common_knns[[idx]]
    sum(knn[idx, ] %in% com_knn) 
  })
  overlap
})
names(overlaps) <- paste0(pa$transformations, "-", pa$alphas)


res <- bind_cols(tibble(dataset = pa$dataset, seed = pa$seed, pca_dim = pa$pca_dim, knn = pa$knn, n_common_knns = n_common_knns),
                 overlaps)
write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))


