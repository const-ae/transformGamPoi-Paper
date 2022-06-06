library(SingleCellExperiment)

pa <- argparser::arg_parser("Download 10X datasets from GEO dataset")
pa <- argparser::add_argument(pa, "--dataset", type = "character", help = "The GSE id of the dataset") 
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

source("src/consistency_benchmark/download_helper.R")
message("Get ", pa$dataset)
sce <- data_loaders[[pa$dataset]]()
UMI <- as.matrix(assay(sce))

saveRDS(UMI, file.path(pa$working_dir, "results", pa$result_id))
