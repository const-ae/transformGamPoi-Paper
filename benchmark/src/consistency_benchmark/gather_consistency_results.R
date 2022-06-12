library(tidyverse)

pa <- argparser::arg_parser("Make a table of the consistency results")
pa <- argparser::add_argument(pa, "--consistency_results", type = "character", help = "A list of consistency results")
pa <- argparser::add_argument(pa, "--file_id", type = "character", help = "The hash of the file")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
# argv <- c("--consistency_results", "040c23eaf6a17-7b8f45a06acb5", "040c23eaf6a17-caf0d73819b69", "--working_dir", "/scratch/ahlmanne/transformation_benchmark/output/")
pa <- argparser::parse_args(pa)

print(pa)
consistency_results <- read_lines(pa$consistency_results)
stopifnot(digest::digest(consistency_results) == pa$file_id)

res <- lapply(consistency_results, function(res_id){
  res <- read_tsv(file.path(pa$working_dir, "results", res_id), col_types = list(alpha = "character"), show_col_types = FALSE)
  duration <- read_tsv(file.path(pa$working_dir, "duration", res$transformation_id), show_col_types = FALSE)
  bind_cols(res, list(cputime_sec = sum(deframe(duration)[c('user.self', "sys.self", "user.child", "sys.child")]), elapsed_sec = deframe(duration)['elapsed']))
}) %>%
  bind_rows()



write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))
