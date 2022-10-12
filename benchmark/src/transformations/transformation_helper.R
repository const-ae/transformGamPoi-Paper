library(Matrix)


.SANITY_FOLDER <- "/home/ahlmanne/prog/Sanity/"


mply_dbl <- function(x, FUN, ncol=1, ...){
  if(is.vector(x)){
    res <- vapply(x, FUN, FUN.VALUE=rep(0.0, times=ncol), ...)
  }else{
    res <- apply(x, 1, FUN, ...) * 1.0
    if(nrow(x) > 0 && length(res) == 0){
      # Empty result, make matrix
      res <- matrix(numeric(0),nrow=0, ncol=nrow(x))
    }else if(nrow(x) == 0){
      res <- matrix(numeric(0), nrow=ncol, ncol=0)
    }
    if((ncol == 1 && ! is.vector(res)) || (ncol > 1 && nrow(res) != ncol)){
      stop(paste0("values must be length ", ncol,
                  ", but result is length ", nrow(res)))
    }
  }
  
  if(ncol == 1){
    as.matrix(res, nrow=length(res), ncol=1)
  }else{
    t(res)
  }
}



knn_for_dist <- function(x, k){
  tmp_mat <- as.matrix(x)
  stopifnot(k < ncol(tmp_mat)-1)
  diag(tmp_mat) <- Inf
  mply_dbl(seq_len(nrow(tmp_mat)), function(idx){
    order(tmp_mat[idx, ])[seq_len(k)]
  }, ncol = k)
}

make_knn_graph <- function(transformation, dat, pca_dim, k_nearest_neighbors){
  if(transformation ==  "sanity_dists"){
    knn_for_dist(dat, k = k_nearest_neighbors)
  }else if(transformation == "glmpca"){
    pca_dim <- min(pca_dim, nrow(dat), ncol(dat))
    pca_res <- glmpca::glmpca(Y = dat, L = pca_dim, fam = if(isFALSE(attr(dat, "alpha"))) "poi" else "nb",
                              sz = attr(dat, "size_factor"), nb_theta = 1/attr(dat, "alpha"),
                              # ctl = list(verbose = TRUE, maxIter = 2, minIter = 1)
    )
    red_dat <- as.matrix(pca_res$factors)
    BiocNeighbors::findAnnoy(red_dat, k = k_nearest_neighbors, warn.ties = FALSE)$index
  }else if(transformation == "newwave"){
    pca_dim <- min(pca_dim, nrow(dat)-1, ncol(dat)-1)
    se <- SummarizedExperiment::SummarizedExperiment(assay = list(counts = dat))
    pca_res <- NewWave::newWave(Y = se, K = pca_dim, n_gene_disp = 100, children = 1)
    red_dat <- SingleCellExperiment::reducedDim(pca_res, "newWave")
    BiocNeighbors::findAnnoy(red_dat, k = k_nearest_neighbors, warn.ties = FALSE)$index
  }else{
    if(pca_dim >= nrow(dat) || pca_dim >= ncol(dat)){ 
      red_dat <- t(dat)
    }else if(pca_dim >= nrow(dat) / 2 || pca_dim >= ncol(dat) / 2){
      red_dat <- BiocSingular::runPCA(t(dat), rank = pca_dim, get.rotation = FALSE, BSPARAM = BiocSingular::ExactParam())$x
    }else{
      red_dat <- BiocSingular::runPCA(t(dat), rank = pca_dim, get.rotation = FALSE, BSPARAM = BiocSingular::FastAutoParam())$x
    }
    BiocNeighbors::findAnnoy(red_dat, k = k_nearest_neighbors, warn.ties = FALSE)$index
  }
}



# Run Sanity
run_sanity <- function(x, variance_min = 0.01, variance_max = 20, n_bins = 116, n_threads = 1){
  # This is helpful to avoid problems with big files
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 20)
  stopifnot(is.matrix(x))
  stopifnot(all(matrixStats::colSums2(x) > 0))
  stopifnot(all(matrixStats::rowSums2(x) > 0))
  # write x to a file and make sure it is deleted after the function finishes
  on.exit({
    if(exists("input")){
      input <- get("input")
      if(file.exists(input)){
        file.remove(input)
      }
    }
  })
  input <- tempfile(pattern = "matrix", fileext = ".tsv")
  output_dir <- tempfile(pattern = "sanity_dir")
  dir.create(output_dir)
  print(paste0("Output dir: ", output_dir))
  
  if(is.null(rownames(x))){
    rownames(x) <- paste0("Gene_", seq_len(nrow(x)))
  }
  if(is.null(colnames(x))){
    colnames(x) <- paste0("Cell_", seq_len(ncol(x)))
  }
  readr::write_tsv(tibble::as_tibble(x, rownames = "GeneID"), input)
  
  code <- system2(file.path(.SANITY_FOLDER, "bin/Sanity"), args = c(paste0("--file ", input),
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
  code <- system2(file.path(.SANITY_FOLDER, "bin/Sanity_distance"), args = c(paste0("--folder ", directory),
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


z_score_rows <- function(x, center = TRUE, scale = TRUE){
  t(scale(t(x), center = center, scale = scale))
}

select_hvg <- function(x, n = 1000){
  row_var <- MatrixGenerics::rowVars(x)
  x[rank(-row_var, ties.method = "first") <= 1000,,drop=FALSE]
}

  
logp1_fnc <- function(UMI, sf, alpha){
  transformGamPoi::shifted_log_transform(UMI, pseudo_count = 1, size_factors = sf)
}

logp1_zscore_fnc <- function(UMI, sf, alpha){
  res <- logp1_fnc(UMI, sf, alpha)
  z_score_rows(res)
}

logp1_hvg_zscore_fnc <- function(UMI, sf, alpha){
  res <- logp1_fnc(UMI, sf, alpha)
  z_score_rows(select_hvg(res))
}

logp1_hvg_fnc <- function(UMI, sf, alpha){
  res <- logp1_fnc(UMI, sf, alpha)
  select_hvg(res)
}

logp1_size_normed_fnc <- function(UMI, sf, alpha){
  x <- transformGamPoi::shifted_log_transform(UMI, pseudo_count = 1, size_factors = sf)
  # This differs from the code in Suppl. Fig. 12 of the 
  # _Depth normalization for single-cell genomics count data_ paper
  # because (I think) they confused the row and colSums function
  # But as we want to normalize the depth per cell, we need to works with the 
  # *colSums*.
  pf <- MatrixGenerics::colSums2(x)
  pf <- pf / mean(pf)
  sweep(x, MARGIN = 2, STATS = pf, FUN = "/")
}

logp_cpm_fnc  <- function(UMI, sf, alpha){
  colsums <- MatrixGenerics::colSums2(UMI)
  log1p(t(t(UMI) / colsums * 1e6))
}

logp_alpha_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::shifted_log_transform(UMI, overdispersion = alpha, size_factors = sf, on_disk = FALSE)
}

acosh_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::acosh_transform(UMI, overdispersion = alpha, size_factors = sf, on_disk = FALSE)
}

pearson_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "pearson", size_factor = sf, on_disk = FALSE)
}

pearson_clip_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "pearson", size_factor = sf, clipping = TRUE, on_disk = FALSE)
}

pearson_clip_zscore_fnc <- function(UMI, sf, alpha){
  res <- pearson_clip_fnc(UMI, sf, alpha)
  z_score_rows(res)
}

pearson_clip_hvg_zscore_fnc <- function(UMI, sf, alpha){
  res <- pearson_clip_fnc(UMI, sf, alpha)
  z_score_rows(select_hvg(res))
}

pearson_clip_hvg_fnc <- function(UMI, sf, alpha){
  res <- pearson_clip_fnc(UMI, sf, alpha)
  select_hvg(res)
}

pearson_analytic_fnc  <- function(UMI, sf, alpha){
  transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "analytic_pearson", size_factor = sf, clipping = TRUE, on_disk = FALSE)
}

sctransform_fnc <- function(UMI, sf, alpha){
  UMI <- as(UMI, "dgCMatrix")
  colnames(UMI) <- paste0("Cell_", seq_len(ncol(UMI)))
  rownames(UMI) <- paste0("Gene_", seq_len(nrow(UMI)))
  res <- sctransform::vst(UMI, vst.flavor = "v2")$y
  stopifnot(ncol(res) == ncol(UMI))
  res
}

rand_quantile_fnc  <- function(UMI, sf, alpha){
  rand_quantile <- transformGamPoi::residual_transform(UMI, overdispersion = alpha, residual_type = "randomized_quantile", size_factor = sf, on_disk = FALSE)
  rand_quantile[matrixStats::rowAnyMissings(rand_quantile), ] <- 0
  rand_quantile
}

dino_fnc <- function(UMI, sf, alpha){
  UMI <- as.matrix(UMI)
  colnames(UMI) <- paste0("Cell_", seq_len(ncol(UMI)))
  rownames(UMI) <- paste0("Gene_", seq_len(nrow(UMI)))
  as.matrix(log1p(Dino::Dino(UMI, nCores = 1)))
}



normalisr_norm_fnc <- function(UMI, sf, alpha){
  UMI <- as.matrix(UMI)
  reticulate::use_virtualenv("normalisr_python_env", required = TRUE)
  norm <- reticulate::import("normalisr.normalisr")
  norm_res <- norm$lcpm(UMI, nth = 1)
  
  covars <- norm$normcov(norm_res[[4]])
  scaling_factors <- norm$scaling_factor(as.matrix(UMI))
  weight <- norm$compute_var(norm_res[[1]],covars)
  # Changed `normmean = FALSE` (the default) to `normmean = TRUE` after
  # private correspondence with Lingfei Wang, the author of `normalisr`.
  normvar_res <- norm$normvar(norm_res[[1]], covars, weight, scaling_factors, nth=1, normmean = TRUE)
  as.matrix(normvar_res[[1]])
}

sanity_map_fnc  <- function(UMI, sf, alpha){
  UMI <- as.matrix(UMI)
  colnames(UMI) <- paste0("Cell_", seq_len(ncol(UMI)))
  rownames(UMI) <- paste0("Gene_", seq_len(nrow(UMI)))
  on.exit({
    if(exists("sanity_map")){
      sanity_map <- get("sanity_map")
      if(file.exists(sanity_map$outdir)){
        file.remove(list.files(sanity_map$outdir, full.names = TRUE))
        file.remove(sanity_map$outdir)
      }
    }
  })
  sanity_map <- run_sanity(UMI, n_threads = 40)
  
  
  sanity_map$mean
}

sanity_dists_fnc <- function(UMI, sf, alpha){
  UMI <- as.matrix(UMI)
  colnames(UMI) <- paste0("Cell_", seq_len(ncol(UMI)))
  rownames(UMI) <- paste0("Gene_", seq_len(nrow(UMI)))
  on.exit({
    if(exists("sanity_map")){
      sanity_map <- get("sanity_map")
      if(file.exists(sanity_map$outdir)){
        file.remove(list.files(sanity_map$outdir, full.names = TRUE))
        file.remove(sanity_map$outdir)
      }
    }
  })
  sanity_map <- run_sanity(UMI, n_threads = 40)
  res <- sanity_distance(sanity_map$outdir, n_threads = 40)
  
  res
}


glmpca_fnc <- function(UMI, sf, alpha){
  attr(UMI, "alpha") <- alpha
  attr(UMI, "size_factor") <- sf
  UMI
}

newwave_fnc <- function(UMI, sf, alpha){
  UMI
}

raw_counts_fnc <- function(UMI, sf, alpha){
  UMI
}

scaled_raw_counts_fnc <- function(UMI, sf, alpha){
  t(t(UMI) / sf)
}




all_transformations <- list(logp1 = logp1_fnc, logp_cpm = logp_cpm_fnc, acosh = acosh_fnc, logp_alpha = logp_alpha_fnc, 
                            logp1_hvg = logp1_hvg_fnc, logp1_zscore = logp1_zscore_fnc, logp1_hvg_zscore = logp1_hvg_zscore_fnc, logp1_size_normed = logp1_size_normed_fnc,
                            pearson = pearson_fnc, pearson_clip = pearson_clip_fnc, pearson_analytic = pearson_analytic_fnc, sctransform = sctransform_fnc, rand_quantile = rand_quantile_fnc,
                            pearson_clip_hvg = pearson_clip_hvg_fnc, pearson_clip_zscore = pearson_clip_zscore_fnc, pearson_clip_hvg_zscore = pearson_clip_hvg_zscore_fnc,
                            dino = dino_fnc, normalisr_norm = normalisr_norm_fnc, sanity_map = sanity_map_fnc, sanity_dists = sanity_dists_fnc,
                            glmpca = glmpca_fnc, newwave = newwave_fnc, 
                            raw_counts = raw_counts_fnc, scaled_raw_counts = scaled_raw_counts_fnc)

