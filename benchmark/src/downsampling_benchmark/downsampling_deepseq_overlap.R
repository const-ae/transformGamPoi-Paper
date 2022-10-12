library(tidyverse)

pa <- argparser::arg_parser("Take the set of KNN graphs and convert them to a long format")
pa <- argparser::add_argument(pa, "--full_knn_result_ids", type = "character", nargs = Inf, help = "The id of a file in output/results") 
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
stopifnot(length(pa$alphas) == length(pa$full_knn_result_ids))



full_KNNs <- lapply(pa$full_knn_result_ids, function(id) readRDS(file.path(pa$working_dir, "results", id)))
# Filter out negative controls
# full_KNNs <- full_KNNs[! pa$transformations %in% c("raw_counts", "scaled_raw_counts")]

stopifnot(nrow(full_KNNs[[1]]) == sapply(full_KNNs, nrow))
stopifnot(ncol(full_KNNs[[1]]) == sapply(full_KNNs, ncol))

n_cells <- nrow(full_KNNs[[1]])
names(full_KNNs) <- pa$transformations

long_neighbor_df <- map_df(full_KNNs, function(index_mat){
  as_tibble(index_mat %>% magrittr::set_colnames(paste0("Nei_", seq_len(ncol(index_mat))))) %>%
    mutate(origin = seq_len(n())) %>%
    pivot_longer(-origin, names_to = "tmp", values_to = "neighbor") %>%
    dplyr::transmute(origin = as.integer(origin), neighbor = as.integer(neighbor)) 
}, .id = "transformation")

write_tsv(long_neighbor_df, file.path(pa$working_dir, "results", pa$result_id))


