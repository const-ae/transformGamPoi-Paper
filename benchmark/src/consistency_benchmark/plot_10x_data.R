
pa <- argparser::arg_parser("Take simulated data and generate the tSNE and UMAP for plotting")
pa <- argparser::add_argument(pa, "--name", type = "character", help = "The name of the dataset")
pa <- argparser::add_argument(pa, "--input_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)

counts <- readRDS(file.path(pa$working_dir, "/results", pa$input_id))

source("src/transformations/transformation_helper.R")
sf <- MatrixGenerics::colSums2(counts)
sf <- sf / mean(sf)
log_counts <- logp1_fnc(counts, sf, alpha = FALSE)
pca_log_counts <- BiocSingular::runPCA(t(log_counts), rank = 20, get.rotation = FALSE, BSPARAM = BiocSingular::FastAutoParam())$x

tsne_log_counts <- scater::calculateTSNE(log_counts)

clustering <- bluster::clusterRows(pca_log_counts, bluster::KNNGraphParam(), full = TRUE)
if(length(unique(clustering$clusters)) > 15){
  clusters <- bluster::mergeCommunities(clustering$objects$graph, clustering$clusters, number = 15)
}else{
  clusters <- clustering$clusters
}



res <- data.frame(name = pa$name, cluster = clusters, col_sums = MatrixGenerics::colSums2(counts),
                  tsne_log_counts_axis1 = tsne_log_counts[,1], tsne_log_counts_axis2 = tsne_log_counts[,2],
                  pca_log_counts_axis1 = pca_log_counts[,1],   pca_log_counts_axis2 = pca_log_counts[,2])

write.table(res, file.path(pa$working_dir, "results", pa$result_id), sep = "\t", row.names = FALSE, quote = FALSE)