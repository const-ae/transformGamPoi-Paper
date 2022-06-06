library(SingleCellExperiment)

pa <- argparser::arg_parser("Download 10X datasets from GEO dataset")
pa <- argparser::add_argument(pa, "--dataset", type = "character", help = "The GSE id of the dataset") 
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = 'the seed used for the simulation')
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--dataset", "mcSCRB", "--result_id", "000000")
pa <- argparser::parse_args(pa)

set.seed(pa$seed)

source("src/downsampling_benchmark/download_helper.R")
message("Get ", pa$dataset)
sce <- data_loaders[[pa$dataset]]()
UMI <- as.matrix(assay(sce))


colsums <- colSums2(UMI)
downsample_proportion <- 5000 / median(colsums)
downsampled_UMI = as.matrix(scuttle::downsampleMatrix(UMI, prop = downsample_proportion, bycol = FALSE))

expressed_cells <- matrixStats::colSums2(downsampled_UMI) > 0
expressed_genes <- matrixStats::rowSums2(downsampled_UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]
downsampled_UMI <- downsampled_UMI[expressed_genes, expressed_cells]

saveRDS(list(full = UMI, reduced = downsampled_UMI), file.path(pa$working_dir, "results", pa$result_id))
