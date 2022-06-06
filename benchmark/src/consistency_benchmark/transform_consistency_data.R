


pa <- argparser::arg_parser("Take a matrix and generate a KNN graph")
pa <- argparser::add_argument(pa, "--transformation", type = "character", help = "The name of the transformation")
pa <- argparser::add_argument(pa, "--input_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--knn", type = "numeric", help = "The number of k nearest neighbors that are considered")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "The dimensions for the pca transformation")
pa <- argparser::add_argument(pa, "--alpha", type = "character", default = "FALSE", help = "The alpha parameter. Ignored by some transformations.")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = 'the seed used for the simulation')
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--transformation", "acosh", "--input_id", "5c007d-f6f0c7", "--knn", 5, "--pca_dim", 10, "--alpha", "0.05", "--result_id", "000000")
pa <- argparser::parse_args(pa)

print(pa)

set.seed(pa$seed)

# get a list called `all_transformations` that contains 
# all transformations as ready to call functions
# Furthermore, it initializes the `make_knn_graph` function
source("src/transformations/transformation_helper.R")


######### Start Transformation #######

UMI <- readRDS(file.path(pa$working_dir, "results", pa$input_id))
expressed_cells <- matrixStats::colSums2(UMI) > 0
expressed_genes <- matrixStats::rowSums2(UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells]

first_gene_half <- sample(seq_len(nrow(UMI)), round(nrow(UMI)/2))
second_gene_half <- setdiff(seq_len(nrow(UMI)), first_gene_half)
UMI_1 <- UMI[first_gene_half,,drop=FALSE]
UMI_2 <- UMI[second_gene_half,,drop=FALSE]

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

sf_1 <- MatrixGenerics::colSums2(UMI_1)
sf_1 <- sf_1 / mean(sf_1)

sf_2 <- MatrixGenerics::colSums2(UMI_2)
sf_2 <- sf_2 / mean(sf_2)

duration <- system.time({
  trans_dat1 <- all_transformations[[pa$transformation]](UMI_1, sf_1, alpha)
  trans_dat2 <- all_transformations[[pa$transformation]](UMI_2, sf_2, alpha)
  
  KNN_1 <- make_knn_graph(pa$transformation, trans_dat1, pa$pca_dim, pa$knn)
  KNN_2 <- make_knn_graph(pa$transformation, trans_dat2, pa$pca_dim, pa$knn)
})

write.table(data.frame(name = names(duration), seconds = as.vector(duration)),
            file.path(pa$working_dir, "duration/", pa$result_id), 
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(list(KNN_1 = KNN_1, KNN_2 = KNN_2), file.path(pa$working_dir, "results/", pa$result_id))
