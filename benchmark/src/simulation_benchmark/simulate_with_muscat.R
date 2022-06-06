

library(SingleCellExperiment)
library(tidyverse)

pa <- argparser::arg_parser("Simulate some data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = 'the seed used for the simulation')
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--seed", "1", "--result_id", "000000")
pa <- argparser::parse_args(pa)

set.seed(pa$seed)

library(muscat)
data(example_sce)
sce_preped <- prepSim(example_sce, verbose = FALSE)


# Make dataset with 
# * 5000 cells
# * 100 genes
# * 2 samples
# * 70% equal expression, 30% DE 
# * and 4 cluster
sim <- muscat::simData(sce_preped, rel_lfc =c(1, 0.5, 0.1, 0.05), nc = 5e3, nk = 4, p_dd = c(0.7, 0, 0.3, 0, 0, 0), lfc = 2, ng = 1e3, force = TRUE)

# Reconstruct ground truth matrix
tmp <- as_tibble(metadata(sim)$gene_info) %>%
  pivot_longer(starts_with("sim_mean"), names_sep = "\\.", names_to = c(".value", "group_id")) %>%
  dplyr::select(gene, cluster_id, group_id, sim_mean) %>%
  full_join(colData(sim) %>% as_tibble(rownames = "cell_id") %>% mutate(cell_id = factor(cell_id, levels = cell_id)),
            by = c("cluster_id", "group_id"))

ground_truth_mat <- tmp %>%
  arrange(cell_id) %>%
  dplyr::select(gene, cell_id, sim_mean) %>%
  pivot_wider(id_cols = gene, names_from = cell_id, values_from = sim_mean) %>%
  column_to_rownames("gene") %>%
  as.matrix()

assay(sim, "LogExpressionMean") <- log10(ground_truth_mat + 1e-4)

saveRDS(list(ground_truth = assay(sim, "LogExpressionMean"), UMI = assay(sim, "counts")), 
        file.path(pa$working_dir, "results/", pa$result_id))