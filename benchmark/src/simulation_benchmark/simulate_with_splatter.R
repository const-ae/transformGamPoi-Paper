
library(SingleCellExperiment)

pa <- argparser::arg_parser("Simulate some data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = 'the seed used for the simulation')
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--seed", "1", "--result_id", "000000")
pa <- argparser::parse_args(pa)

set.seed(pa$seed)

source("src/consistency_benchmark/download_helper.R")

sce <- get_GSE150068_data()

sce <- scuttle::logNormCounts(sce)
colData(sce)$cluster_id <- scran::quickCluster(sce)

n_genes <- nrow(sce)
n_cells <- ncol(sce)
# n_genes <- 4000
# n_cells <- 5000
reference_data_counts <- assay(sce)[seq_len(n_genes), seq_len(n_cells)]
reference_data_counts <- reference_data_counts[,colSums2(reference_data_counts) > 0]
reference_data_counts <- reference_data_counts[rowSums2(reference_data_counts) > 0,]


params <- splatter::splatEstimate(as.matrix(reference_data_counts))
# Simulate Paths (with strong batch effect)
sim <- splatter::splatSimulatePaths(params, de.prob = 0.2, de.facLoc = 1.8, 
                                    batchCells = c(3000, 2200, 4800), batch.facLoc = 0.15,
                                    lib.loc = 8, lib.scale = 0.1)
# Ground truth (get rid of the exact zeros)
assay(sim, "LogCellMeans") <- log10(t(t(assay(sim, "CellMeans")) / sim$ExpLibSize) + 1e-10)
UMI <- as.matrix(assay(sim, "counts"))
# Filter out problematic columns and rows
expressed_cells <- colSums2(UMI) > 10
expressed_genes <- rowSums2(UMI) > 0
sim <- sim[, expressed_cells]
sim <- sim[expressed_genes, ]

saveRDS(list(ground_truth = assay(sim, "LogCellMeans"), UMI = assay(sim, "counts")), 
        file.path(pa$working_dir, "results/", pa$result_id))