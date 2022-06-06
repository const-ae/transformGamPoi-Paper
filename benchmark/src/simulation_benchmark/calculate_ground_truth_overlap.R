library(tidyverse)

pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--knn_id", type = "character", help = "The id of a file in output/results")
pa <- argparser::add_argument(pa, "--simulator", type = "character", help = "[Just for documentation purposes] A readable identifier of the data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "[Just for documentation purposes] The seed used to tame randomness")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "[Just for documentation purposes] The number of PCA dimensions used before KNN graph construction")
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "[Just for documentation purposes] The number of nearest neighbors considered")
pa <- argparser::add_argument(pa, "--transformation", type = "character", help = "[Just for documentation purposes] A readable identifier of the transformation")
pa <- argparser::add_argument(pa, "--alpha", type = "character", help = "[Just for documentation purposes] Specification of the overdispersion.")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)

ground_truth <- readRDS(file.path(pa$working_dir, "results", pa$data_id))$ground_truth
knn <- readRDS(file.path(pa$working_dir, "results", pa$knn_id))
stopifnot(ncol(ground_truth) == nrow(knn))

groung_truth_knn <- BiocNeighbors::findAnnoy(t(ground_truth), k = pa$knn, warn.ties = FALSE, get.distance = FALSE)$index

knn_overlap <- mean(sapply(seq_len(nrow(knn)), function(cell_idx){
  length(intersect(groung_truth_knn[cell_idx,], knn[cell_idx,]))
}))


ground_truth_knn_graph <- bluster::neighborsToKNNGraph(groung_truth_knn)
ground_truth_clustering_obj <- igraph::cluster_walktrap(ground_truth_knn_graph)
ground_truth_clustering <- factor(igraph::membership(ground_truth_clustering_obj))

n_clusters <- length(levels(ground_truth_clustering))
print(paste0("n_clusters = ", n_clusters))

knn_graph <- bluster::neighborsToKNNGraph(knn)
tmp_clustering_obj <- igraph::cluster_walktrap(knn_graph)
cluster_assignment <- factor(igraph::membership(tmp_clustering_obj))
n_clusters_counts <- length(levels(cluster_assignment))
print(paste0("n_clusters_counts = ", n_clusters_counts))

ARI <- aricode::ARI(ground_truth_clustering, cluster_assignment)
AMI <- aricode::AMI(ground_truth_clustering, cluster_assignment)
NMI <- aricode::NMI(ground_truth_clustering, cluster_assignment)

res <- tibble(ARI = ARI, AMI = AMI, NMI = NMI, mean_knn_overlap = knn_overlap,
              n_clusters = n_clusters, n_clusters_counts = n_clusters_counts,
              ground_truth_id = pa$data_id, transformation_id = pa$knn_id, simulator = pa$simulator, 
              seed = pa$seed, pca_dim = pa$pca_dim, knn = pa$knn, transformation = pa$transformation, alpha = pa$alpha)
write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))
