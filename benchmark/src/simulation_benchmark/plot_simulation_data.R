
pa <- argparser::arg_parser("Take simulated data and generate the tSNE and UMAP for plotting")
pa <- argparser::add_argument(pa, "--simulator", type = "character", help = "The name of the simulator")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = "The seed that was used for the simulator")
pa <- argparser::add_argument(pa, "--input_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--simulator", "muscat", "--input_id", "1d86f65b3f270-2e5de9d89e26e", "--working_dir", "/scratch/ahlmanne/transformation_benchmark/output/", "--result_id", "000000")
pa <- argparser::parse_args(pa)

print(pa)


# Get ground truth and simulated data
ground_truth <- readRDS(file.path(pa$working_dir, "/results", pa$input_id))$ground_truth
counts <- readRDS(file.path(pa$working_dir, "/results", pa$input_id))$UMI

source("src/transformations/transformation_helper.R")
sf <- MatrixGenerics::colSums2(counts)
sf <- sf / mean(sf)
log_counts <- logp1_fnc(counts, sf, alpha = FALSE)
pca_log_counts <- BiocSingular::runPCA(t(log_counts), rank = 2, get.rotation = FALSE, BSPARAM = BiocSingular::FastAutoParam())$x



pca_gt <- scater::calculatePCA(ground_truth, ncomponents = 2)
tsne_gt <- scater::calculateTSNE(ground_truth)
tsne_log_counts <- scater::calculateTSNE(log_counts)

clustering <- bluster::clusterRows(t(ground_truth), bluster::KNNGraphParam(), full = TRUE)
if(length(unique(clustering$clusters)) > 15){
  clusters <- bluster::mergeCommunities(clustering$objects$graph, clustering$clusters, number = 15)
}else{
  clusters <- clustering$clusters
}


res <- data.frame(simulator = pa$simulator, cluster = clusters, col_sums = MatrixGenerics::colSums2(counts),
                  tsne_ground_truth_axis1 = tsne_gt[,1],       tsne_ground_truth_axis2 = tsne_gt[,2],
                  tsne_log_counts_axis1 = tsne_log_counts[,1], tsne_log_counts_axis2 = tsne_log_counts[,2],
                  pca_ground_truth_axis1 = pca_gt[,1],         pca_ground_truth_axis2 = pca_gt[,2],
                  pca_log_counts_axis1 = pca_log_counts[,1],   pca_log_counts_axis2 = pca_log_counts[,2])

write.table(res, file.path(pa$working_dir, "results", pa$result_id), sep = "\t", row.names = FALSE, quote = FALSE)