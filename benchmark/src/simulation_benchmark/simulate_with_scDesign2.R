

library(SingleCellExperiment)

pa <- argparser::arg_parser("Simulate some data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = 'the seed used for the simulation')
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--seed", "1", "--result_id", "000000")
pa <- argparser::parse_args(pa)

set.seed(pa$seed)

source("src/consistency_benchmark/download_helper.R")


sce <- get_GSE130931_data()

colData(sce)$cluster_id <- scran::quickCluster(sce, min.size = 20)


n_genes <- nrow(sce)
n_cells <- ncol(sce)
sce <- sce[seq_len(n_genes), seq_len(n_cells)]
sce <- sce[,colSums2(assay(sce)) > 0]
sce <- sce[rowSums2(assay(sce)) > 0,]

# For plotting
sce <- scuttle::logNormCounts(sce)


# More Clusters, make sure that there are no singlets
mat <- counts(sce)
colnames(mat) <- as.character(sce$cluster_id)
cluster_names <- unique(as.character(sce$cluster_id))
fit <- scDesign2::fit_model_scDesign2(mat, cell_type_sel = cluster_names, sim_method = "copula", marginal = "nb")
UMI <- scDesign2::simulate_count_scDesign2(fit, n_cell_new = ncol(mat), cell_type_prop = table(colnames(mat))[cluster_names] / ncol(mat), sim_method = "copula")

# Reconstruct ground truth
ground_truth <- do.call(cbind, lapply(fit[cluster_names], function(fi) {
  mat <- matrix(NA, nrow = nrow(mat), ncol = fi$n_cell)
  mat[fi$gene_sel1, ] <- fi$marginal_param1[,3]
  mat[fi$gene_sel2, ] <- fi$marginal_param2[,3]
  mat[fi$gene_sel3, ] <- 1e-8
  mat
}))

sim <- SingleCellExperiment(assays = list(counts = UMI, LogExpressionMean = log10(ground_truth)), colData = data.frame(cluster = colnames(UMI)))
sim$seq_depth <- colSums2(counts(sim))

# Save simulated data
saveRDS(list(ground_truth = assay(sim, "LogExpressionMean"), UMI = assay(sim, "counts")), 
        file.path(pa$working_dir, "results/", pa$result_id))

