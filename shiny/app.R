library(shiny)
library(tidyverse)
library(gggroupedscale)
library(zeallot)
source("../notebooks/utils.R")
source("../notebooks/annotation_helper.R")


res <- bind_rows(
  read_tsv("../benchmark/output/benchmark_results/simulation_results.tsv") %>% 
    transmute(benchmark = "simulation", overlap = mean_knn_overlap, knn, pca_dim, alpha = as.character(alpha), transformation, dataset = simulator, replicate = seed, cputime_sec, elapsed_sec),
  read_tsv("../benchmark/output/benchmark_results/consistency_results.tsv") %>% 
    transmute(benchmark = "consistency", overlap = mean_overlap, knn, pca_dim, alpha = as.character(alpha), transformation, dataset, replicate = seed, cputime_sec, elapsed_sec),
  read_tsv("../benchmark/output/benchmark_results/downsampling_results.tsv") %>% 
    transmute(benchmark = "downsampling", overlap = overlap, knn, pca_dim, alpha = as.character(alpha), transformation, dataset, replicate = seed, cputime_sec = full_cputime_sec, elapsed_sec = full_elapsed_sec)
) %>%
  mutate(transformation = factor(transformation, levels = trans_families$transformation)) %>%
  mutate(alpha = ifelse(alpha == "FALSE", "0", alpha))

tmp <- read_rds("../benchmark/output/benchmark_results/dataset_plot_data.RDS")
reduced_dim_data <- bind_rows(
  bind_rows(tmp$downsampling) %>%
    pivot_longer(-c(name, cluster, col_sums_full, col_sums_reduced), names_pattern = "(tsne|pca)_(ground_truth|log_counts)_(full|reduced)_(axis\\d)", names_to = c("dim_red_method", "origin", "downsampling", ".value")) %>%
    pivot_longer(c(col_sums_full, col_sums_reduced), names_pattern = "(.+)_(full|reduced)", names_to = c(".value", "downsampling2")) %>%
    filter(downsampling == downsampling2) %>% 
    dplyr::select(-downsampling2)%>%
    dplyr::rename(dataset = name),
  bind_rows(tmp$consistency) %>%
    pivot_longer(-c(name, cluster, col_sums), names_pattern = "(tsne|pca)_(ground_truth|log_counts)_(axis\\d)", names_to = c("dim_red_method", "origin", ".value")) %>%
    dplyr::rename(dataset = name),
  bind_rows(tmp$simulation) %>%
    pivot_longer(-c(simulator, cluster, col_sums), names_pattern = "(tsne|pca)_(ground_truth|log_counts)_(axis\\d)", names_to = c("dim_red_method", "origin", ".value")) %>%
    dplyr::rename(dataset = simulator)
) %>%
  left_join(enframe(dataset_benchmark, name = "dataset", value = "benchmark"))



source("utils.R")
source("option_pane.R")
source("benchmark_plot.R")
source("duration_plot.R")
source("contrast_plot.R")
source("reduced_dim_plot.R")

ui <- fluidPage(
  includeCSS("www/main.css"),
  h1("Transformation and Preprocessing of Single-Cell RNA-Seq Data"),
  h2("Online Supplementary Information"),
  h3("Benchmark Results"),
  benchmarkPlotUI("benchmark"),
  optionPaneUI("benchmark_options"),
  h3("Contrasts"),
  contrastPlotUI("contrasts"),
  optionPaneUI("contrast_options"),
  h3("Computational Expenses"),
  durationPlotUI("duration"),
  optionPaneUI("duration_options", show_detailed_pcadim_selector = FALSE),
  h3("Dataset Exploration"),
  reducedDimPlotUI("reducedDimPlots")
)

server <- function(input, output, session) {

  op_bench <- optionPaneServer("benchmark_options", data = res)
  benchmarkPlotServer("benchmark", data = res, pcadim_sel = op_bench$pca_sel,
                      knn_sel = op_bench$knn_sel, alpha_sel = op_bench$alpha_sel, dataset_sel = op_bench$dataset_sel)
  
  op_contr <- optionPaneServer("contrast_options", data = res)
  contrastPlotServer("contrasts", data = res, pcadim_sel = op_contr$pca_sel, 
                      knn_sel = op_contr$knn_sel, dataset_sel = op_contr$dataset_sel)

  op_dur <- optionPaneServer("duration_options", data = res)
  durationPlotServer("duration", data = res, pcadim_sel = op_dur$pca_sel,
                      knn_sel = op_dur$knn_sel, alpha_sel = op_dur$alpha_sel, dataset_sel = op_dur$dataset_sel)

  reducedDimPlotServer("reducedDimPlots", data = reduced_dim_data)
  
}
shinyApp(ui, server)
