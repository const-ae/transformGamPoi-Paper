# This was executed on the server
print("Make Linear Simulated Tree")

random_id <- paste0(sample(letters, size = 5, replace = TRUE), collapse = "")
print(paste0("id: ", random_id))


library(SingleCellExperiment)
library(tidyverse)

output_dir <- paste0("/home/ahlmanne/projects/transformGamPoi-Paper/intermediate/benchmark/linear_random_tree/", random_id)

if(! dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
  
}

# Copy this script to output folder for backup and reference
file.copy(from = "/home/ahlmanne/projects/transformGamPoi-Paper/src/2021-05-14_Benchmark_Linear_Tree.R",
          to = file.path(output_dir, "generating_script_backup.R"), overwrite = TRUE)

# Convert ID to 
tmp <- as.numeric(charToRaw(random_id))/10
seed <- round(prod(tmp) + sum(cumsum(tmp)))
print(paste0("seed: ", seed))
set.seed(seed)

# Load Baron Dataset as a reference
sce <- scRNAseq::BaronPancreasData("human")
n_genes <- nrow(sce)
n_cells <- ncol(sce)
# n_genes <- 100
# n_cells <- 75
baron_counts <- assay(sce)[seq_len(n_genes), seq_len(n_cells)]


# Make random walk tree
delta_true <- matrix(NA, n_genes, n_cells)
parents <- rep(NA, n_cells)

branch_length <- 500
branch_idx <- 0
start_point <- NULL
end_point <- NULL

for(idx in seq_len(n_cells)-1){
  if(idx == 0){
    parents[idx+1] <- -1
    start_point <- rnorm(n_genes, mean = 0, sd = 1)
    end_point <- rnorm(n_genes, mean = 0, sd = 1)
  }else if(idx %% branch_length == 0){
    start_id <- sample.int(idx, size = 1)
    parents[idx+1] <- start_id - 1
    # Make new start and end point
    branch_idx <- 0
    start_point <- rnorm(n_genes, mean = delta_true[,  parents[idx]], sd = 1)
    end_point <- rnorm(n_genes, mean = 0, sd = 1)
  }else{
    # Continue with idx-1 as parent
    branch_idx <- branch_idx + 1
    parents[idx+1] <- idx-1
  }
  delta <- rnorm(n_genes, mean = start_point + (end_point - start_point) * branch_idx / branch_length, sd = 0.1)
  delta_true[,idx+1] <- delta
}



# Copied from https://github.com/jmbreda/Sanity/blob/94e7063027cb1cd0368134395bfb501e1f8b8377/reproducibility/run_Simulations.m
# N_c = sum(T{:,:},1);
N_c <- colSums2(baron_counts)
# ng <- ng / mean(ng)

# mu_tilde_g = log(sum(SC{:,:},2)./sum(sum(SC{:,:})));
mu_tilde_g <- log(rowSums2(baron_counts) / sum(baron_counts))

# Note the confusing rate vs scale parametrization...
# sig2_g = exprnd(2,N_gene,1);
sig2_g <- rexp(n_genes, rate = 1/2)


# lambda = sqrt( sig2_g./var(delta_true,0,2) );
lambda <- sqrt( sig2_g / rowVars(delta_true) )


# delta_true <- matrix(rnorm(n_genes * n_cells, mean = 0, sd = 1), nrow = n_genes, ncol = n_cells)
# delta_true = lambda.*(delta_true-mean(delta_true,2));
delta_true <- (delta_true - rowMeans2(delta_true)) * lambda

# mu_g = mu_tilde_g - sig2_g/2;
mu_g <- mu_tilde_g - sig2_g / 2

# UMI = poissrnd( N_c.*exp(E) );
UMI <- matrix(rnbinom(n_genes * n_cells, mu = t(t(exp(mu_g + delta_true)) * N_c), size = 1/0.01), n_genes, n_cells)

# Compare properties of generated dataset and Baron
print("UMI")
unname(quantile(c(UMI), c(0.5, 0.99)))
median(colSums2(UMI))
median(rowMeans(UMI))

print("Baron")
unname(quantile(c(as.vector(baron_counts)), c(0.5, 0.99)))
median(colSums2(baron_counts))
median(rowMeans(baron_counts))


# Filter out problematic columns and rows
UMI <- UMI[, colSums2(UMI) > 10]
UMI <- UMI[rowSums2(UMI) > 0, ]

# Save simulated data
saveRDS(parents, file.path(output_dir, "parents.RDS"))
saveRDS(delta_true, file.path(output_dir, "delta_true.RDS"))
saveRDS(UMI, file.path(output_dir, "UMI.RDS"))

print("Saved simulation data")

# Run Sanity
run_sanity <- function(x, variance_min = 0.01, variance_max = 20, n_bins = 116, n_threads = 1){
  stopifnot(is.matrix(x))
  stopifnot(all(colSums2(x) > 0))
  stopifnot(all(rowSums2(x) > 0))
  # write x to a file
  input <- tempfile(pattern = "matrix", fileext = ".tsv")
  output_dir <- file.path("/scratch/ahlmanne/Sanity_output/", 
                          paste0("sanity_output_", paste0(sample(letters, 7, replace = TRUE), collapse = "")))
  dir.create(output_dir)
  print(paste0("Output dir: ", output_dir))
  
  if(is.null(rownames(x))){
    rownames(x) <- paste0("Gene_", seq_len(nrow(x)))
  }
  if(is.null(colnames(x))){
    colnames(x) <- paste0("Cell_", seq_len(ncol(x)))
  }
  readr::write_tsv(tibble::as_tibble(x, rownames = "GeneID"), input)
  
  code <- system2("/home/ahlmanne/prog/Sanity/bin/Sanity", args = c(paste0("--file ", input),
                                                                    paste0("--destination ", output_dir),
                                                                    paste0("--variance_min ", variance_min),
                                                                    paste0("--variance_max ", variance_max),
                                                                    paste0("--number_of_bins ", n_bins),
                                                                    paste0("--n_threads ", n_threads),
                                                                    "--extended_output true"))
  if(code != 0){
    stop("Sanity failed with status ", code)
  }
  
  
  res <- readr::read_tsv(file.path(output_dir, "log_transcription_quotients.txt"), col_types = readr::cols()) 
  mat <- as.matrix(res[,-1,drop=FALSE])
  rownames(mat) <- res[[1]]
  colnames(mat) <- colnames(res)[-1]
  
  errror_bars_df <- readr::read_tsv(file.path(output_dir, "ltq_error_bars.txt"), col_types = readr::cols()) 
  error_bars <- as.matrix(errror_bars_df[,-1,drop=FALSE])
  rownames(error_bars) <- errror_bars_df[[1]]
  colnames(error_bars) <- colnames(errror_bars_df)[-1]
  
  
  list(mean = mat, sd = error_bars, outdir = output_dir)
}

sanity_distance <- function(directory, take_error_into_account = TRUE, signal_to_noise_cutoff = 1, n_threads = 1){
  if(signal_to_noise_cutoff < 0.000001){
    stop("signal_to_noise_cutoff is too small")
  }
  code <- system2("/home/ahlmanne/prog/Sanity/bin/Sanity_distance", args = c(paste0("--folder ", directory),
                                                                             paste0("--with_error_bars ", if(take_error_into_account) "true" else "false"),
                                                                             paste0("--signal_to_noise_cutoff ", signal_to_noise_cutoff),
                                                                             paste0("--n_threads ", n_threads)))
  if(code != 0){
    stop("Sanity failed with status ", code)
  }
  
  s2n_string <- stringr::str_remove(sprintf("%f", signal_to_noise_cutoff), "[0.]*$")
  
  res_file <- file.path(directory, paste0("cell_cell_distance_", if(take_error_into_account) "with_errorbar" else "euclidean", "_s2n_gt_", s2n_string, ".txt"))
  
  dist_vec <- readr::read_tsv(res_file, col_names = "dist", col_types = readr::cols())[[1]] 
  n_cells <- 1/2 + sqrt(1/4 + length(dist_vec) * 2)
  mat <- matrix(NA, nrow = n_cells, ncol = n_cells)
  mat[lower.tri(mat)] <- dist_vec
  as.dist(mat)
}

system.time(
  sane <- run_sanity(UMI, n_threads = 10)
)
system.time(
  sanity_dists <- sanity_distance(sane$outdir, n_threads = 10)
)

saveRDS(sane, file = file.path(output_dir, "sane.RDS"))
saveRDS(sanity_dists, file = file.path(output_dir, "sanity_dists.RDS"))
print("Finished sanity")



sf <- colSums2(UMI)
sf <- sf / mean(sf)
summary(sf)

system.time(
  logp1 <- log(t(t(UMI) / sf) + 1)
)
print("Finished log transformation")

alpha <- 0.05
system.time(
  logp_alpha <- log(t(t(UMI) / sf) + 1/(4 * alpha))
)
print("Finished log plus alpha transformation")

system.time(
  acosh <- 1 / sqrt(alpha) * acosh(2 * alpha * t(t(UMI) / sf) + 1)
)
print("Finished acosh")

system.time(
  pearson <- transformGamPoi::transformGamPoi(UMI, residual_type = "pearson", size_factor = sf)
)
print("Finished pearson")

system.time(
  rand_quantile <- transformGamPoi::transformGamPoi(UMI, residual_type = "randomized_quantile", size_factor = sf)
)
print("Finished randomized quantile")


saveRDS(logp1, file.path(output_dir, "logp1.RDS"))
saveRDS(logp_alpha, file.path(output_dir, "logp_alpha.RDS"))
saveRDS(acosh, file.path(output_dir, "acosh.RDS"))
saveRDS(pearson, file.path(output_dir, "pearson.RDS"))
saveRDS(rand_quantile, file.path(output_dir, "rand_quantile.RDS"))

print("Saved everything")

print("Finished")

