library(tidyverse)

pa <- argparser::arg_parser("Make a table of the consistency results")
pa <- argparser::add_argument(pa, "--dataset_names", type = "character", nargs = Inf, help = "A list of the names of the datasets")
pa <- argparser::add_argument(pa, "--dataset_ids", type = "character", nargs = Inf, help = "A list of ids that link to the datasets")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)

stopifnot(length(pa$dataset_names) == length(pa$dataset_ids))


make_summary <- function(UMI, id, name, benchmark){
  UMI <- UMI[MatrixGenerics::rowSums2(UMI) > 0, MatrixGenerics::colSums2(UMI) > 0]
  glm_fit_global <- glmGamPoi::glm_gp(UMI, overdispersion = "global", overdispersion_shrinkage = FALSE, on_disk = ncol(UMI) >= 20000, verbose = TRUE)
  
  quantiles <- quantile(UMI, c(0.25, 0.5, 0.75, 0.9, 0.99), names = FALSE)
  col_sums <- MatrixGenerics::colSums2(UMI)
  
  sf_perc_rank <- dplyr::percent_rank(glm_fit_global$size_factors)
  UMI_sf_filtered <- UMI[,0.25 < sf_perc_rank & sf_perc_rank <= 0.75]
  row_means <- MatrixGenerics::rowMeans2(UMI_sf_filtered)
  row_vars <- MatrixGenerics::rowVars(UMI_sf_filtered)
  
  tibble(id = id, name = name, benchmark = benchmark,
         n_cells = ncol(UMI), n_genes = nrow(UMI),
         proportion_zeros = sum(UMI == 0) / (n_cells * n_genes),
         quantile025_count = quantiles[1], quantile050_count = quantiles[2], quantile075_count = quantiles[3],
         quantile090_count = quantiles[4], quantile099_count = quantiles[5],
         max_count = max(UMI), 
         `mean_sequencing-depth` = mean(col_sums),
         `median_sequencing-depth` = median(col_sums),
         `min_sequencing-depth` = min(col_sums),
         `max_sequencing-depth` = max(col_sums),
         `overdispersion-global` = glm_fit_global$overdispersions[1],
         row_means = list(row_means),
         row_vars = list(row_vars)
  )
}

res <- map_dfr(seq_along(pa$dataset_names), function(idx){
  name <- pa$dataset_names[idx]
  id <- pa$dataset_ids[idx]
  input <- readRDS(file.path(pa$working_dir, "/results", id))
  if(is.null(names(input))){
    make_summary(input, id = id, name = name, benchmark = "consistency")
  }else if(names(input)[1] == "full"){
    bind_rows(make_summary(input$full, id = id, name = name, benchmark = "downsampling_full"),
              make_summary(input$reduced, id = id, name = name, benchmark = "downsampling_reduced"))
  }else{
    make_summary(input$UMI, id = id, name = name, benchmark = "simulation")
  }
})

# write_tsv does not like list columns
res <- res %>%
  mutate(row_means = map_chr(row_means, ~ paste0(.x, collapse = ",")),
         row_vars = map_chr(row_vars, ~ paste0(.x, collapse = ",")))


write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))

