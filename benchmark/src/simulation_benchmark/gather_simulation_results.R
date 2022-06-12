library(tidyverse)

pa <- argparser::arg_parser("Make a table of the consistency results")
pa <- argparser::add_argument(pa, "--simulation_results", type = "character", help = "A list of simulation results")
pa <- argparser::add_argument(pa, "--file_id", type = "character", help = "The hash of the file")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--simulation_results", "5fd46c2395015-671e2cf01815b", "--working_dir", "/scratch/ahlmanne/transformation_benchmark/output/")
pa <- argparser::parse_args(pa)

print(pa)

simulation_results <- read_lines(pa$simulation_results)
stopifnot(digest::digest(simulation_results) == pa$file_id)

res <- lapply(simulation_results, function(res_id){
  res <- read_tsv(file.path(pa$working_dir, "results", res_id), col_types = list(alpha = "character"), show_col_types = FALSE)
  duration <- read_tsv(file.path(pa$working_dir, "duration", res$transformation_id), show_col_types = FALSE)
  bind_cols(res, list(cputime_sec = sum(deframe(duration)[c('user.self', "sys.self", "user.child", "sys.child")]), elapsed_sec = deframe(duration)['elapsed']))
}) %>%
  bind_rows()



write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))
