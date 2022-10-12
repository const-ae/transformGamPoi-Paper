
library(tidyverse)
source("src/job_management_utils.R")

benchmark_complete <- TRUE

##### Benchmark Results #####
consistency <- readRDS(file.path("output/job_storage/consistency_job.RDS"))
simulation <- readRDS(file.path("output/job_storage/simulation_job.RDS"))
downsampling <- readRDS(file.path("output/job_storage/downsampling_job.RDS"))
downsampling_best_of <- readRDS(file.path("output/job_storage/downsampling_best_of_job.RDS"))


# Which consistency jobs failed and why?
if(job_status(consistency) == "done"){
  message("save consistency")
  file.copy(file.path(.OUTPUT_FOLDER, "results", consistency$result_id), "output/benchmark_results/consistency_results.tsv", overwrite = TRUE)
}else{
  benchmark_complete <- FALSE
  job_stats <- sapply(consistency$dependencies, job_status)
  message("Job stats for Consistency")
  print(table(job_stats))
}
if(job_status(simulation) == "done"){
  message("save simulation")
  file.copy(file.path(.OUTPUT_FOLDER, "results", simulation$result_id), "output/benchmark_results/simulation_results.tsv", overwrite = TRUE)
}else{
  benchmark_complete <- FALSE
  job_stats <- sapply(simulation$dependencies, job_status)
  message("Job stats for Simulation")
  print(table(job_stats))
}
if(job_status(downsampling) == "done"){
  message("save downsampling")
  file.copy(file.path(.OUTPUT_FOLDER, "results", downsampling$result_id), "output/benchmark_results/downsampling_results.tsv", overwrite = TRUE)
}else{
  benchmark_complete <- FALSE
  job_stats <- sapply(downsampling$dependencies, job_status)
  message("Job stats for downsampling")
  print(table(job_stats))
}
if(job_status(downsampling_best_of) == "done"){
  message("save downsampling_best_of")
  file.copy(file.path(.OUTPUT_FOLDER, "results", downsampling_best_of$result_id), "output/benchmark_results/downsampling_best_of_results.tsv", overwrite = TRUE)
}else{
  benchmark_complete <- FALSE
  job_stats <- sapply(downsampling_best_of$dependencies, job_status)
  message("Job stats for downsampling_best_of")
  print(table(job_stats))
}


##### Plots #####
con_plots <- readRDS(file.path("output/job_storage/consistency_plotting_jobs.RDS"))
sim_plots <- readRDS(file.path("output/job_storage/simulation_plotting_job.RDS"))
dow_plots <- readRDS(file.path("output/job_storage/downsampling_plotting_jobs.RDS"))

if(all(sapply(sim_plots, job_status) == "done") && all(sapply(con_plots, job_status) == "done") && all(sapply(dow_plots, job_status) == "done")){
  message("save plot_data")
  plots <- list(
    simulation = lapply(sim_plots, function(job){
      read_tsv(file.path(.OUTPUT_FOLDER, "results", job$result_id), show_col_types = FALSE)
    }),
    consistency = lapply(con_plots, function(job){
      read_tsv(file.path(.OUTPUT_FOLDER, "results", job$result_id), show_col_types = FALSE)
    }),
    downsampling = lapply(dow_plots, function(job){
      read_tsv(file.path(.OUTPUT_FOLDER, "results", job$result_id), show_col_types = FALSE)
    })
  )
  saveRDS(plots, file.path("output/benchmark_results/dataset_plot_data.RDS"))
}else{
  benchmark_complete <- FALSE
  message("Job stats for Consistency Plots")
  print(table(sapply(con_plots, job_status)))
  message("Job stats for Simulation Plots")
  print(table(sapply(sim_plots, job_status)))
  message("Job stats for Downsampling Plots")
  print(table(sapply(dow_plots, job_status)))
}

##### Deeply sequenced data overlap #####
deepseq_overlap_jobs <- readRDS(file.path("output/job_storage/deepseq_overlap_jobs.RDS"))
if(all(sapply(deepseq_overlap_jobs, job_status) == "done")){
  lapply(deepseq_overlap_jobs, function(job){
    read_tsv(file.path(.OUTPUT_FOLDER, "results", job$result_id), show_col_types = FALSE, lazy = FALSE)
  }) %>%
    magrittr::set_names(sapply(deepseq_overlap_jobs, \(x) get_params(x)$dataset)) %>%
    write_rds("output/benchmark_results/downsampling_deepseq_overlap.RDS", compress = "gz")
}else{
  benchmark_complete <- FALSE
  message("Job stats for Deeply sequenced data overlap")
  print(table(sapply(deepseq_overlap_jobs, job_status)))
}

##### Dataset summaries #####
dataset_summary <- readRDS(file.path("output/job_storage/dataset_summary.RDS"))
if(job_status(dataset_summary) == "done"){
  message("save dataset_summary")
  file.copy(file.path(.OUTPUT_FOLDER, "results", dataset_summary$result_id), "output/benchmark_results/dataset_summaries.tsv", overwrite = TRUE)
}else{
  message("Job stats for dataset_summary:")
  print(job_status(dataset_summary))
}

#### Overlap per cell #####
simulation_strat_jobs <- readRDS(file.path("output/job_storage/simulation_for_stratification_jobs.RDS"))
if(all(sapply(simulation_strat_jobs, job_status) == "done")){
  message("save simulation_strat_jobs")
  bind_rows(map_df(simulation_strat_jobs, function(x){
    read_tsv(file.path(.OUTPUT_FOLDER, "results", x$result_id), show_col_types = FALSE)
  })) %>%
    write_tsv("output/benchmark_results/simulation_stratefication_results.tsv")
}else{
  benchmark_complete <- FALSE
  message("Job stats for simulation_strat_jobs:")
  print(table(sapply(simulation_strat_jobs, job_status)))
}

consistency_strat_jobs <- readRDS(file.path("output/job_storage/consistency_for_stratification_jobs.RDS"))
if(all(sapply(consistency_strat_jobs, job_status) == "done")){
  message("save consistency_strat_jobs")
  bind_rows(map_df(consistency_strat_jobs, function(x){
    read_tsv(file.path(.OUTPUT_FOLDER, "results", x$result_id), show_col_types = FALSE)
  })) %>%
    write_tsv("output/benchmark_results/consistency_stratefication_results.tsv")
}else{
  benchmark_complete <- FALSE
  message("Job stats for consistency_strat_jobs:")
  print(table(sapply(consistency_strat_jobs, job_status)))
}

downsampling_strat_jobs <- readRDS(file.path("output/job_storage/downsampling_for_stratification_jobs.RDS"))
if(all(sapply(downsampling_strat_jobs, job_status) == "done")){
  message("save downsampling_strat_jobs")
  bind_rows(map_df(downsampling_strat_jobs, function(x){
    read_tsv(file.path(.OUTPUT_FOLDER, "results", x$result_id), show_col_types = FALSE)
  })) %>%
    write_tsv("output/benchmark_results/downsampling_stratefication_results.tsv")
}else{
  benchmark_complete <- FALSE
  message("Job stats for downsampling_strat_jobs:")
  print(table(sapply(downsampling_strat_jobs, job_status)))
}


if(benchmark_complete){
  print("The benchmark is complete. All results were stored in 'output/benchmark_results'")
}else{
  print("!!!!!!!!The benchmark is not complete!!!!!\nOnly some results were stored in 'output/benchmark_results'")
}

