# This is a simple script to aggregate all the results from the benchmark into one convenient table

library(tidyverse)
library(tidylog, warn.conflicts = FALSE)

method_annotation <- tibble(
  method = c("logp1", "logp_cpm", "logp_alpha", "acosh", "pearson", "pearson_clip", "rand_quantile", "sctransform", "sanity_main",
             "pca_logp1", "pca_logp_cpm", "pca_logp_alpha", "pca_acosh", "pca_pearson", "pca_pearson_clip", "pca_rand_quantile", "pca_sctransform", "pca_sanity_main",
             "sanity_dists")) %>%
  mutate(group = case_when(
    str_starts(method, "pca_") ~ "PCA",
    method == "sanity_dists" ~ "Sanity",
    TRUE ~ "Transformation"
  ))

bind_rows(
  map_df(list.files("data/benchmark/random_tree/", full.names = TRUE), ~ read_tsv(.x, col_types = cols())) %>%
    mutate(dataset = "random_tree"),
  map_df(list.files("data/benchmark/linear_random_tree/", full.names = TRUE), ~ read_tsv(.x, col_types = cols())) %>%
    mutate(dataset = "linear_random_tree")
) %>%
  left_join(method_annotation) %>%
  # I do a lot of double saving of all non-PCA results, so I have to get rid of them again
  mutate(pcadims = ifelse(group == "PCA", pcadims, NA)) %>%
  distinct() %>%
  # Fix order of method
  mutate(method = factor(method, levels = method_annotation$method)) %>%
  # Pick columns
  transmute(script_hash, id, method, transformation_method = preproc_method, dataset, group, preproc_method, pcadims, knn, mean, median, min, max,
            transformation_cputime = preproc_duration_user + preproc_duration_user_child + preproc_duration_system + preproc_duration_sys_child,
            knnsearch_cputime = duration_user + duration_user_child + duration_system + duration_sys_child,
            transformation_elapsed = preproc_duration_elapsed,
            knnsearch_elapsed = duration_elapsed) %>% 
  # Save
  write_tsv("data/benchmark_results.tsv")



# Session Info
sessionInfo()

