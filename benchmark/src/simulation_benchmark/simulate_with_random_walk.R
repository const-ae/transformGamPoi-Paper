
library(SingleCellExperiment)

pa <- argparser::arg_parser("Simulate some data")
pa <- argparser::add_argument(pa, "--seed", type = "numeric", help = 'the seed used for the simulation')
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--seed", "1", "--result_id", "000000")
pa <- argparser::parse_args(pa)

set.seed(pa$seed)




# Load Baron Dataset as a reference
sce <- scRNAseq::BaronPancreasData("human")

sce <- scuttle::logNormCounts(sce)
colData(sce)$cluster_id <- scran::quickCluster(sce)



n_genes <- nrow(sce)
n_cells <- ncol(sce)
# n_genes <- 4000
# n_cells <- 500
reference_data_counts <- assay(sce)[seq_len(n_genes), seq_len(n_cells)]
reference_data_counts <- reference_data_counts[,colSums2(reference_data_counts) > 0]
reference_data_counts <- reference_data_counts[rowSums2(reference_data_counts) > 0,]
n_genes <- nrow(reference_data_counts)
n_cells <- ncol(reference_data_counts)

## Simulation


# Make random walk tree
delta_true <- matrix(NA, n_genes, n_cells)
parents <- rep(NA, n_cells)

branch_length <- 13

# Make random walk tree
for(idx in seq_len(n_cells)){
  # See section S1.2 in Breda paper
  if(idx == 1){
    delta <- rnorm(n_genes, mean = 0, sd = 1)
    parents[idx] <- 0
  }else if(idx %% branch_length == 0){
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


# Copied and adapted from https://github.com/jmbreda/Sanity/blob/94e7063027cb1cd0368134395bfb501e1f8b8377/reproducibility/run_Simulations.m
# N_c = sum(T{:,:},1);
N_c <- colSums2(reference_data_counts)
# ng <- ng / mean(ng)

# mu_tilde_g = log(sum(SC{:,:},2)./sum(sum(SC{:,:})));
mu_tilde_g <- log(rowSums2(reference_data_counts) / sum(reference_data_counts))

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

# Save simulated data
saveRDS(list(ground_truth = delta_true, UMI = UMI), 
        file.path(pa$working_dir, "results/", pa$result_id))



