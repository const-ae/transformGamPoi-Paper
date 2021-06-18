library(tidyverse)
library(tidylog)
`%|%` <- rlang::`%|%`

# Get Random_tree

dat <- list.files("../output/benchmark/random_tree", pattern = "random_tree-full_run", full.names = TRUE) %>%
  map_df(~ dplyr::mutate(filename = .x, read_tsv(.x, col_types = cols()))) %>%
  dplyr::select(filename, everything()) %>%
  mutate(filename = as_tibble(str_match(filename, 
                                        pattern = ".+random_tree-full_run-id_([a-z]{5})-knn_(\\d+)-pcadim_(\\d+)-script_([0-9a-z]{6})-results.tsv")[,-1,drop=FALSE], 
                              .name_repair = ~ c("id.x", "knn", "pcadim", "script_hash.x"))) %>%
  unpack(filename) %>%
  filter(id.x == id, script_hash.x == script_hash) %>%
  dplyr::select(- ends_with(".x")) %>%
  mutate(knn = as.integer(knn), pcadim = as.integer(pcadim))


method_annotation <- tibble(method = c("logp1", "logp_alpha", "acosh", "pearson", "rand_quantile" ,
                                       "pca_logp1", "pca_logp_alpha", "pca_acosh", "pca_pearson", "pca_rand_quantile", 
                                       "sanity_dists")) %>%
  mutate(group = case_when(
    str_starts(method, "pca_") ~ "PCA",
    str_starts(method, "sanity_") ~ "Sanity",
    TRUE ~ "Transformation"
  ))



# Get Linear Random_tree

dat3 <- list.files("../output/benchmark/linear_random_tree", pattern = "linear_random_tree-full_run", full.names = TRUE) %>%
  map_df(~ dplyr::mutate(filename = .x, read_tsv(.x, col_types = cols()))) %>%
  dplyr::select(filename, everything()) %>%
  mutate(filename = as_tibble(str_match(filename, 
                                        pattern = ".+linear_random_tree-full_run-id_([a-z]{5})-knn_(\\d+)-pcadim_(\\d+)-script_([0-9a-z]{6})-results.tsv")[,-1,drop=FALSE], 
                              .name_repair = ~ c("id.x", "knn", "pcadim", "script_hash.x"))) %>%
  unpack(filename) %>%
  filter(id.x == id, script_hash.x == script_hash) %>%
  dplyr::select(- ends_with(".x")) %>%
  mutate(knn = as.integer(knn), pcadim = as.integer(pcadim))


# Duration


extract_user_time_for_method <- function(content, method){
  as.numeric(str_match(content, paste0('\n\\s*user\\s+system\\s+elapsed\\s*\n\\s*(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s*\n\\[1\\] \"', method, '\"'))[,2])
}


duration_df <- tibble(files = list.files("../intermediate/benchmark/", pattern = "slurm-.+.out", full.names = TRUE)) %>%
  mutate(modified_date = as.Date(file.info(files)$mtime)) %>%
  # filter(lubridate::ymd(modified_date) == "2021-05-27" | lubridate::ymd(modified_date) == "2021-05-28") %>%
  mutate(content = map_chr(files, ~ paste0(read_lines(.x), collapse ="\n"))) %>%
  mutate(id = str_match(content, '\n\\[1\\] \"id: ([a-z]{5})\"')[,2]) %>%
  mutate(durations = map(content, function(co){
    tibble(method = c("sanity_dists", "logp1", "logp_alpha", "acosh", "pearson", "rand_quantile",  
                      "pca_logp1", "pca_logp_alpha", "pca_acosh", "pca_pearson", "pca_rand_quantile"),
           preproc_duration = extract_user_time_for_method(co, c(paste0("Finished ", c("sanity", "log transformation", "log plus alpha transformation", "acosh", "pearson", "randomized quantile")), 
                                                                 paste0("Finished ", c(          "log transformation", "log plus alpha transformation", "acosh", "pearson", "randomized quantile")))))
  })) %>%
  unnest(durations) %>%
  mutate(sanity_transformation_dur = as.numeric(str_match(content, '\nPrint extended output\\s*\n\\s*user\\s+system\\s+elapsed\\s*\n\\s*(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s+(\\d+\\.\\d+)\\s*\nSanity folder')[,2])) %>%
  mutate(preproc_duration = preproc_duration +  ifelse(str_starts(method, "sanity_"), sanity_transformation_dur, 0)) %>%
  dplyr::select(-c(content, sanity_transformation_dur))




bind_rows(mutate(dat, data = "random_tree"),
          mutate(dat3, data = "linear_random_tree")) %>%
  left_join(method_annotation) %>%
  mutate(method = factor(method, levels = method_annotation$method)) %>%
  filter(script_hash == "6e2dd3") %>%
  left_join(duration_df, by = c("id", "method")) %>%
  mutate(total_duration = duration + preproc_duration) %>%
  dplyr::select(-c(preproc_duration, duration)) %>%
  write_tsv("../output/benchmark_results.tsv")






