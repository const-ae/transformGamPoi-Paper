
pa <- argparser::arg_parser("Take simulated data and generate the tSNE and UMAP for plotting")
pa <- argparser::add_argument(pa, "--name", type = "character", help = "The name of the dataset")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "The seed that was used for the downsampling")
pa <- argparser::add_argument(pa, "--input_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--simulator", "muscat", "--input_id", "1d86f65b3f270-2e5de9d89e26e", "--working_dir", "/scratch/ahlmanne/transformation_benchmark/output/", "--result_id", "000000")
pa <- argparser::parse_args(pa)

print(pa)


# Get ground truth and simulated data
counts_full <- readRDS(file.path(pa$working_dir, "/results", pa$input_id))$full
counts_reduced <- readRDS(file.path(pa$working_dir, "/results", pa$input_id))$reduced

log_counts_full <- scuttle::normalizeCounts(counts_full)
log_counts_reduced <- scuttle::normalizeCounts(counts_reduced)

pca_log_counts_full <- scater::calculatePCA(log_counts_full, ncomponents = 2)
pca_log_counts_reduced <- scater::calculatePCA(log_counts_reduced, ncomponents = 2)

tsne_log_counts_full <- scater::calculateTSNE(log_counts_full)
tsne_log_counts_reduced <- scater::calculateTSNE(log_counts_reduced)

clustering <- bluster::clusterRows(t(log_counts_full), bluster::KNNGraphParam(), full = TRUE)
if(length(unique(clustering$clusters)) > 15){
  clusters <- bluster::mergeCommunities(clustering$objects$graph, clustering$clusters, number = 15)
}else{
  clusters <- clustering$clusters
}


res <- data.frame(name = pa$name, cluster = clusters, 
                  col_sums_full = MatrixGenerics::colSums2(counts_full),
                  col_sums_reduced = MatrixGenerics::colSums2(counts_reduced),
                  tsne_log_counts_full_axis1 = tsne_log_counts_full[,1],       tsne_log_counts_full_axis2 = tsne_log_counts_full[,2],
                  tsne_log_counts_reduced_axis1 = tsne_log_counts_reduced[,1], tsne_log_counts_reduced_axis2 = tsne_log_counts_reduced[,2],
                  pca_log_counts_full_axis1 = pca_log_counts_full[,1],         pca_log_counts_full_axis2 = pca_log_counts_full[,2],
                  pca_log_counts_reduced_axis1 = pca_log_counts_reduced[,1],   pca_log_counts_reduced_axis2 = pca_log_counts_reduced[,2])

write.table(res, file.path(pa$working_dir, "results", pa$result_id), sep = "\t", row.names = FALSE, quote = FALSE)