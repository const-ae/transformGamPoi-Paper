


pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--transformation", type = "character", help = "The name of the transformation")
pa <- argparser::add_argument(pa, "--input_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--data_mode", type = "character", help = "Either 'full' or 'reduced'") 
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "The number of k nearest neighbors that are considered")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "The dimensions for the pca transformation")
pa <- argparser::add_argument(pa, "--alpha", type = "character", default = "FALSE", help = "The alpha parameter. Ignored by some transformations.")
pa <- argparser::add_argument(pa, "--output_folder", type = "character", help = "The folder where the results and everything else is stored")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--transformation", "acosh", "--input_id", "5c007d-f6f0c7", "--knn", 5, "--pca_dim", 10, "--alpha", "0.05", "--result_id", "000000")
pa <- argparser::parse_args(pa)

print(pa)
stopifnot(pa$data_mode %in% c("full", "reduced"))

# get a list called `all_transformations` that contains 
# all transformations as ready to call functions
# Furthermore, it initializes the `make_knn_graph` function
source("src/transformations/transformation_helper.R")


######### Start Transformation #######

if(pa$data_mode == "full"){
  UMI <- readRDS(file.path(pa$working_dir, "results", pa$input_id))$full
}else{
  UMI <- readRDS(file.path(pa$working_dir, "results", pa$input_id))$reduced
}
expressed_cells <- matrixStats::colSums2(UMI) > 0
expressed_genes <- matrixStats::rowSums2(UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]

alpha <- pa$alpha
if(pa$alpha == "global"){
  alpha <- "global"
}else if(! is.na(suppressWarnings(readr::parse_double(pa$alpha, na = character(0L))))){
  alpha <- readr::parse_double(pa$alpha)
}else if(! is.na(suppressWarnings(readr::parse_logical(pa$alpha, na = character(0L))))){
  alpha <- readr::parse_logical(pa$alpha)
}else{
  stop("Cannot parse alpha=", alpha)
}

sf <- MatrixGenerics::colSums2(UMI)
sf <- sf / mean(sf)

duration <- system.time({
  trans_dat <- all_transformations[[pa$transformation]](UMI, sf, alpha)
  KNN <- make_knn_graph(pa$transformation, trans_dat, pa$pca_dim, pa$knn)
})

write.table(data.frame(name = names(duration), seconds = as.vector(duration)),
            file.path(pa$working_dir, "duration", pa$result_id), 
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(KNN, file.path(pa$working_dir, "results", pa$result_id))
