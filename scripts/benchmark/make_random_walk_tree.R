# This was executed on the server
print("Make Random Walk Tree")

random_id <- paste0(sample(letters, size = 5, replace = TRUE), collapse = "")
print(paste0("id: ", random_id))


library(SingleCellExperiment)
library(tidyverse)

output_dir <- paste0("/home/ahlmanne/projects/transformGamPoi-Paper/intermediate/benchmark/random_tree/", random_id)

if(! dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
  
}

# Copy this script to output folder for backup and reference
file.copy(from = "/home/ahlmanne/projects/transformGamPoi-Paper/src/2021-07-05_Make_Random_Walk_Tree.R",
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
for(idx in seq_len(n_cells)){
  # See section S1.2 in Breda paper
  if(idx == 1){
    delta <- rnorm(n_genes, mean = 0, sd = 1)
    parents[idx] <- 0
  }else if(idx %% 13 == 0){
    # Choose new parent
    parents[idx] <- sample.int(idx - 1, size = 1)
    parent_delta <- delta_true[,  parents[idx]]
    delta <- rnorm(n_genes, mean = parent_delta, sd = 1)
  }else{
    # Continue with idx-1 as parent
    parents[idx] <- idx-1
    parent_delta <- delta_true[, parents[idx]]
    delta <- rnorm(n_genes, mean = parent_delta, sd = 1)
  }
  delta_true[,idx] <- delta
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
UMI <- matrix(rpois(n_genes * n_cells, lambda = t(t(exp(mu_g + delta_true)) * N_c)), n_genes, n_cells)

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
expressed_cells <- colSums2(UMI) > 10
expressed_genes <- rowSums2(UMI) > 0
UMI <- UMI[, expressed_cells]
UMI <- UMI[expressed_genes, ]
delta_true <- delta_true[, expressed_cells]
delta_true <- delta_true[expressed_genes, ]

# Save simulated data
saveRDS(delta_true, file.path(output_dir, "delta_true.RDS"))
saveRDS(UMI, file.path(output_dir, "UMI.RDS"))

print("Saved simulation data")

sessionInfo()

print("Finished")

