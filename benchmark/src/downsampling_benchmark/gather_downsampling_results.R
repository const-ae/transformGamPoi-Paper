library(tidyverse)

pa <- argparser::arg_parser("Make a table of the consistency results")
pa <- argparser::add_argument(pa, "--downsampling_results", type = "character", nargs = Inf, help = "A list of consistency results")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--downsampling_results", "b18d75fb35378-3f68221c07e97", "--working_dir", "/scratch/ahlmanne/transformation_benchmark/output/")
pa <- argparser::parse_args(pa)

print(pa)

res <- lapply(pa$downsampling_results, function(res_id){
  res <- read_tsv(file.path(pa$working_dir, "results", res_id), show_col_types = FALSE)
  duration_full <- vapply(res$transformation_full_data_ids, function(id){
    sum(read_tsv(file.path(pa$working_dir, "duration", id), show_col_types = FALSE)$seconds)
  }, numeric(1L))
  duration_reduced <- vapply(res$transformation_reduced_data_ids, function(id){
    sum(read_tsv(file.path(pa$working_dir, "duration", id), show_col_types = FALSE)$seconds)
  }, numeric(1L))
  bind_cols(res, list(duration_full_sec = duration_full, duration_reduced_sec = duration_reduced))
}) %>%
  bind_rows()



write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))
