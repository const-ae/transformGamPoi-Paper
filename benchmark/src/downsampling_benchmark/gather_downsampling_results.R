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
  duration_full <- map_df(res$transformation_full_data_ids, function(id){
    durf <- read_tsv(file.path(pa$working_dir, "duration", id), show_col_types = FALSE)
    tibble(full_cputime_sec = sum(deframe(durf)[c('user.self', "sys.self", "user.child", "sys.child")]), full_elapsed_sec = deframe(durf)['elapsed'])
  })
  duration_reduced <- map_df(res$transformation_reduced_data_ids, function(id){
    durr <- read_tsv(file.path(pa$working_dir, "duration", id), show_col_types = FALSE)
    tibble(reduced_cputime_sec = sum(deframe(durr)[c('user.self', "sys.self", "user.child", "sys.child")]), reduced_elapsed_sec = deframe(durr)['elapsed'])
  })
  bind_cols(res, duration_full, duration_reduced)
}) %>%
  bind_rows()



write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))
