

library(SingleCellExperiment)

pa <- argparser::arg_parser("Simulate some data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = 'the seed used for the simulation')
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--seed", "1", "--result_id", "000000")
pa <- argparser::parse_args(pa)

set.seed(pa$seed)





options(Ncpus = 1L) # change this to the number of cores in your system
options(dyngen_download_cache_dir = "/scratch/ahlmanne/benchmark/data")

num_cells <- 5000
num_features <- 1000
backbone <- dyngen::backbone_consecutive_bifurcating()
num_tfs <- nrow(backbone$module_info)
num_targets <- round((num_features - num_tfs) / 2)
num_hks <- num_features - num_targets - num_tfs
config <- dyngen::initialise_model(backbon = backbone,
                                   num_tfs = num_tfs,
                                   num_targets = num_targets,
                                   num_hks = num_hks,
                                   num_cells = num_cells)
# calls generate_tf_network, generate_feature_network, generate_kinetics, 
# generate_gold_standard, generate_cells, generate_experiment
library(Matrix)
res <- dyngen::generate_dataset(config, format = "sce")
sim <- res$dataset

# Filter out problematic columns and rows
UMI <- counts(sim)
expressed_cells <- colSums2(UMI) > 10
expressed_genes <- rowSums2(UMI) > 0
sim <- sim[, expressed_cells]
sim <- sim[expressed_genes, ]

# Save simulated data
saveRDS(list(ground_truth = t(reducedDim(sim, "MDS")), UMI = assay(sim, "counts")), 
        file.path(pa$working_dir, "results/", pa$result_id))