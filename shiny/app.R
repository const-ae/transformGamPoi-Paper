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

div_with_floating_gear <- function(panel, menu_content){
  div(class = "floating_gear_holder",
      div(class = "floating_gear", shinyWidgets::dropdownButton(menu_content, icon = icon("gear"), size = "sm", right = TRUE)),
      panel)
}

source("option_pane.R")

cons_res <- filter(res, benchmark == "consistency") %>%
  group_by(dataset, transformation, knn, pca_dim, alpha) %>%
  summarize(overlap = mean(overlap))

ui <- fluidPage(
  includeCSS("www/main.css"),
  pcadimPlotUI("pca_sel"),
)
server <- function(input, output, session) {
  z <- pcadimPlotServer("pca_sel", data = cons_res, pca_sel = reactive(20))[[1]]
  observeEvent(z(), {
    print(paste0("The value of z has changed: ", z()))
  })
}
shinyApp(ui, server)
