library(tidyverse)

.OUTPUT_FOLDER <- "/scratch/ahlmanne/transformation_benchmark2/output/"
source("src/job_management_utils.R")

instructions <- yaml::read_yaml("job_overview.yaml")
si <- instructions[["simulation"]]
si2 <- instructions[["simulation_for_stratification"]]
ci <- instructions[["consistency"]]
ci2 <- instructions[["consistency_for_stratification"]]
di <- instructions[["downsampling"]]
di2 <- instructions[["downsampling_best_of"]]
di3 <- instructions[["downsampling_for_stratification"]]

##############################
#        Helper Function     #
##############################
get_data_for_simulation_benchmark <- function(simulator, seed){
  wrap_script(paste0("src/simulation_benchmark/simulate_with_", simulator, ".R"), 
              params = list(seed = as.integer(seed)),
              duration = "0-10:00:00", memory = "40GB")
}

get_10X_data_for_consistency_benchmark <- function(dataset_name){
  wrap_script("src/consistency_benchmark/download_10X_datasets.R", 
              params = list(dataset = dataset_name),
              duration = "0-01:00:00", memory = "40GB")
}

get_data_for_downsampling_benchmark <- function(dataset_name, seed){
  wrap_script("src/downsampling_benchmark/download_deeply_sequenced_datasets.R", 
              params = list(dataset = dataset_name, seed = seed),
              duration = "0-01:00:00", memory = "40GB")
}


plot_simulation_data <- function(simulator, seed, dataset){
  params <- list(simulator = simulator, seed = seed, input_id = dataset$result_id)
  wrap_script("src/simulation_benchmark/plot_simulation_data.R", params = params,
              dependencies = list(dataset), duration = "0-05:00:00", memory = "40GB")
}

plot_10x_data <- function(name, dataset){
  params <- list(name = name, input_id = dataset$result_id)
  wrap_script("src/consistency_benchmark/plot_10x_data.R", params = params,
              dependencies = list(dataset), duration = "0-05:00:00", memory = "120GB")
}

plot_smartseq3_data <- function(name, seed, dataset){
  params <- list(name = name, seed = seed, input_id = dataset$result_id)
  wrap_script("src/downsampling_benchmark/plot_smartseq3_data.R", params = params,
              dependencies = list(dataset), duration = "0-05:00:00", memory = "40GB")
}

collect_dataset_summary_statistics <- function(names, datasets){
  params <- list(dataset_names = names, dataset_ids = map_chr(datasets, "result_id"))
  wrap_script("src/transformations/collect_dataset_summary_statistics.R", params = params,
              dependencies = datasets, duration = "0-05:00:00", memory = "120GB")
}


transform_simulated_data <- function(transformation, data, knn, pca_dim, alpha){
  params <- list(transformation = transformation, input_id = data$result_id, pca_dim = pca_dim, knn = knn, alpha = alpha)
  if(transformation %in% c("sanity_map", "sanity_dists")){
    duration <- "1-00:00:00"
    n_cpus <- 40
    memory <- "120GB"
  }else if(transformation %in% c("dino", "glmpca", "newwave")){
    duration <- "2-00:00:00"
    n_cpus <- 1
    memory <- "80GB"
  }else{
    duration <- "0-05:00:00"
    n_cpus <- 1
    memory <- "80GB"
  }
  wrap_script("src/simulation_benchmark/transform_simulated_data.R", 
              params = params,
              dependencies = list(data),
              duration = duration, memory = memory, n_cpus = n_cpus)
}

transform_deeply_sequenced_data <- function(transformation, data, data_mode = c("full", "reduced"), knn, pca_dim, alpha){
  data_mode <- match.arg(data_mode)
  params <- list(transformation = transformation, input_id = data$result_id, data_mode = data_mode, pca_dim = pca_dim, knn = knn, alpha = alpha)
  if(transformation %in% c("sanity_map", "sanity_dists")){
    duration <- "1-00:00:00"
    n_cpus <- 40
  }else if(transformation %in% c("dino", "glmpca", "newwave")){
    duration <- "2-00:00:00"
    n_cpus <- 1
  }else{
    duration <- "0-05:00:00"
    n_cpus <- 1
  }
  wrap_script("src/downsampling_benchmark/transform_downsampling_data.R", 
              params = params,
              dependencies = list(data),
              duration = duration, memory = "80GB", n_cpus = n_cpus)
}

transform_consistency_data <- function(transformation, data, knn, pca_dim, alpha, seed){
  params <- list(transformation = transformation, input_id = data$result_id, pca_dim = pca_dim, knn = knn, alpha = alpha, seed = seed)
  if(transformation %in% c("sanity_map", "sanity_dists")){
    duration <- "8-00:00:00"
    n_cpus <- 40
    memory <- "150GB"
  }else if(transformation %in% c("dino", "glmpca", "newwave")){
    duration <- "8-00:00:00"
    n_cpus <- 1
    memory <- "150GB"
  }else{
    duration <- "0-05:00:00"
    n_cpus <- 1
    memory <- "80GB"
  }
  wrap_script("src/consistency_benchmark/transform_consistency_data.R", 
              params = params,
              dependencies = list(data),
              duration = duration, memory = memory, n_cpus = n_cpus)
}

calculate_10X_consistency <- function(knn_construction_job,
                                       dataset, seed, pca_dim, knn, transformation, alpha){
  wrap_script("src/consistency_benchmark/calculate_10X_consistency.R", 
              params = list(input_id = knn_construction_job$result_id,
                            dataset = dataset, seed = seed, pca_dim = pca_dim, knn = knn, 
                            transformation = transformation, alpha = alpha),
              dependencies = list(knn_construction_job),
              duration = "0-00:10:00", memory = "60GB")
}

calculate_10X_consistency_for_stratification <- function(data_job, transformed_data, seed, pca_dim, knn, dataset_name){
  wrap_script("src/consistency_benchmark/calculate_10X_consistency_for_stratification.R", 
              params = list(dataset_id = data_job$result_id,
                            input_ids = sapply(transformed_data, function(x) x$id),
                            transformation_names = sapply(transformed_data, function(x) x$name),
                            seed = seed, pca_dim = pca_dim, knn = knn, dataset_name = dataset_name),
              dependencies = c(list(data_job), lapply(transformed_data, function(x) x$dep)),
              duration = "0-01:00:00", memory = "60GB")
}

calculate_downsampling_agreement <- function(full_knns, reduced_knns, dataset, seed, pca_dim, knn, transformations, alphas){
  wrap_script("src/downsampling_benchmark/calculate_downsampling_agreement.R", 
              params = list(full_knn_result_ids = vapply(full_knns, function(e) e$result_id, character(1L)),
                            reduced_knn_result_ids = vapply(reduced_knns, function(e) e$result_id, character(1L)),
                            dataset = dataset, seed = seed, 
                            pca_dim = pca_dim, knn = knn, 
                            transformations = transformations, alphas = alphas),
              dependencies = c(full_knns, reduced_knns),
              duration = "0-00:10:00", memory = "60GB")
}

calculate_downsampling_agreement_for_stratification <- function(full_knns, reduced_knns, dataset, seed, pca_dim, knn, transformations, alphas){
  wrap_script("src/downsampling_benchmark/calculate_downsampling_agreement_for_stratification.R", 
              params = list(full_knn_result_ids = vapply(full_knns, function(e) e$result_id, character(1L)),
                            reduced_knn_result_ids = vapply(reduced_knns, function(e) e$result_id, character(1L)),
                            dataset = dataset, seed = seed, 
                            pca_dim = pca_dim, knn = knn, 
                            transformations = transformations, alphas = alphas),
              dependencies = c(full_knns, reduced_knns),
              duration = "0-01:00:00", memory = "60GB")
}

calculate_ground_truth_overlap <- function(knn_job, data_job,
                                           simulator, seed, pca_dim, knn, transformation, alpha){
  wrap_script("src/simulation_benchmark/calculate_ground_truth_overlap.R", 
              params = list(data_id = data_job$result_id, knn_id = knn_job$result_id,
                            simulator = simulator, seed = seed, pca_dim = pca_dim, knn = knn, 
                            transformation = transformation, alpha = alpha),
              dependencies = list(knn_job, data_job),
              duration =  "1-00:00:00", memory = "60GB")
}

calculate_ground_truth_overlap_for_stratification <- function(data_job, transformed_data,
                                                             simulator, seed, pca_dim, knn){
  wrap_script("src/simulation_benchmark/calculate_ground_truth_overlap_for_stratification.R", 
              params = list(data_id = data_job$result_id, 
                            knn_ids = sapply(transformed_data, function(x) x$id),
                            transformation_names = sapply(transformed_data, function(x) x$name),
                            simulator = simulator, seed = seed, pca_dim = pca_dim, knn = knn),
              dependencies = c(list(data_job), lapply(transformed_data, function(x) x$dep)),
              duration =  "0-01:00:00", memory = "60GB")
}

gather_simulation_results <- function(simulation_jobs){
  # As there is an input length limit, I cannot call the script with 24,000 arguments.
  # Instead I save them to a file and pass the location
  
  result_ids <- vapply(simulation_jobs, function(e) e$result_id, character(1L))
  id <- digest::digest(result_ids)
  dep_id_file <- paste0("output/tmp/simulation_result_file_names_", id)
  write_lines(result_ids, dep_id_file)
  
  
  wrap_script("src/simulation_benchmark/gather_simulation_results.R", 
              params = list(simulation_results = dep_id_file,
                            file_id = id),
              dependencies = simulation_jobs,
              duration = "0-03:00:00", memory = "40GB")
}

gather_consistency_results <- function(consistency_jobs){
  # As there is an input length limit, I cannot call the script with 24,000 arguments.
  # Instead I save them to a file and pass the location
  result_ids <- vapply(consistency_jobs, function(e) e$result_id, character(1L))
  id <- digest::digest(result_ids)
  dep_id_file <- paste0("output/tmp/consistency_result_file_names_", id)
  write_lines(result_ids, dep_id_file)
  
  
  wrap_script("src/consistency_benchmark/gather_consistency_results.R", 
              params = list(consistency_results = dep_id_file, file_id = id),
              dependencies = consistency_jobs,
              duration = "0-03:00:00", memory = "40GB")
}

gather_downsampling_results <- function(downsampling_jobs){
  wrap_script("src/downsampling_benchmark/gather_downsampling_results.R", 
              params = list(downsampling_results = vapply(downsampling_jobs, function(e) e$result_id, character(1L))),
              dependencies = downsampling_jobs,
              duration = "0-00:20:00", memory = "10GB")
}


datasets_to_summarize <- list()
##############################
#  Run Simulation benchmarks #
##############################

simulation_jobs <- list()
sim_plotting_jobs <- list()
# Loop over all simulations, seeds, transformations, alphas, knn and pca settings
for(simulator in si$input_data$simulator){
  for(seed in si$input_data$seed){
    dataset <- get_data_for_simulation_benchmark(simulator, seed)
    if(seed == 1) {
      sim_plotting_jobs <- append(sim_plotting_jobs, list(plot_simulation_data(simulator, seed, dataset)))
      datasets_to_summarize <- append(datasets_to_summarize, list(list(name = simulator, dataset = dataset)))
    }
    for(pca_dim in si$knn_construction$pca){
      for(knn in si$knn_construction$knn){
        for(trans in si$knn_construction$transformations){
          for(a in trans$alpha){
            if(trans$name == "newwave" && pca_dim > 300){
              # skip
            }else{
              knn_constr <- transform_simulated_data(trans$name, dataset, pca_dim = pca_dim, knn = knn, alpha = a)
              res <- calculate_ground_truth_overlap(knn_constr, dataset,
                                                    simulator = simulator, seed = seed, pca_dim = pca_dim,
                                                    knn = knn, transformation = trans$name, alpha = a)
              simulation_jobs <- append(simulation_jobs, list(res))
            }
          }
        }
      }
    }
  }
  message("Added all jobs for ", simulator)
}

simulation <- gather_simulation_results(simulation_jobs)
table(sapply(simulation$dependencies, job_status))
# Store jobs and run them
saveRDS(simulation, file.path("output/job_storage/simulation_job.RDS"))
simulation <- run_job(simulation, "normal")
saveRDS(sim_plotting_jobs, file.path("output/job_storage/simulation_plotting_job.RDS"))
for(job in sim_plotting_jobs) run_job(job)





##################################
# Run 10X Consistency benchmarks #
##################################

consistency_jobs <- list()
consistency_plotting_jobs <- list()
# Loop over all datasets, transformations, knn and pca settings
for(dataset_name in ci$input_data$dataset){
  for(seed in ci$input_data$seed){
    dataset <- get_10X_data_for_consistency_benchmark(dataset_name)
    if(seed == 1) {
      consistency_plotting_jobs <- append(consistency_plotting_jobs, list(plot_10x_data(dataset_name, dataset)))
      datasets_to_summarize <- append(datasets_to_summarize, list(list(name = dataset_name, dataset = dataset)))
    }
    for(pca_dim in ci$knn_construction$pca){
      for(knn in ci$knn_construction$knn){
        for(trans in ci$knn_construction$transformations){
          for(a in trans$alpha){
            if(trans$name == "newwave" && dataset_name == "GSE184806" && pca_dim >  50){
              # Skip (newwave doesn't handle large pca_dim well)              
            }else{
              knn_constr <- transform_consistency_data(trans$name, dataset, pca_dim = pca_dim, knn = knn, alpha = a, seed = seed)
              consistency_metric <- calculate_10X_consistency(knn_constr,
                                                              dataset = dataset_name, seed = seed, pca_dim = pca_dim,
                                                              knn = knn, transformation = trans$name, alpha = a)
              consistency_jobs <- append(consistency_jobs, list(consistency_metric))
            }
          }
        }
      }
    }
  }
  message("Added all jobs for ", dataset_name)
}

consistency <- gather_consistency_results(consistency_jobs)
table(sapply(consistency$dependencies, job_status))
# Store jobs and run them
saveRDS(consistency, file.path("output/job_storage", "consistency_job.RDS"))
consistency <- run_job(consistency, "normal")
saveRDS(consistency_plotting_jobs, file.path("output/job_storage/consistency_plotting_jobs.RDS"))
for(job in consistency_plotting_jobs) run_job(job)

##############################
# Run Downsampling benchmarks#
##############################

downsampling_plotting_jobs <- list()
downsampling_jobs <- list()
# Loop over all datasets, seeds, transformations, knn and pca settings
for(dataset_name in di$input_data$dataset){
  for(seed in di$input_data$seed){
    dataset <- get_data_for_downsampling_benchmark(dataset_name, seed)
    if(seed == 1) {
      downsampling_plotting_jobs <- append(downsampling_plotting_jobs, list(plot_smartseq3_data(dataset_name, seed, dataset)))
      datasets_to_summarize <- append(datasets_to_summarize, list(list(name = dataset_name, dataset = dataset)))
    }
    for(pca_dim in di$knn_construction$pca){
      for(knn in di$knn_construction$knn){
          full_knn_transformations <- list()
          reduced_knn_transformations <- list()
          for(trans in di$knn_construction$transformations){
            stopifnot(length(trans$alpha) == 1)
            knn_full <- transform_deeply_sequenced_data(trans$name, dataset, data_mode = "full", pca_dim = pca_dim, knn = knn, alpha = trans$alpha)
            knn_reduced <- transform_deeply_sequenced_data(trans$name, dataset, data_mode = "reduced", pca_dim = pca_dim, knn = knn, alpha = trans$alpha)
            full_knn_transformations <- append(full_knn_transformations, list(knn_full))
            reduced_knn_transformations <- append(reduced_knn_transformations, list(knn_reduced))
          }
          res <- calculate_downsampling_agreement(full_knn_transformations, reduced_knn_transformations,
                                                  dataset = dataset_name, seed = seed, pca_dim = pca_dim, knn = knn,
                                                  alphas = vapply(di$knn_construction$transformations, function(e) e$alpha, character(1L)),
                                                  transformations = vapply(di$knn_construction$transformations, function(e) e$name, character(1L)))
          downsampling_jobs <- append(downsampling_jobs, list(res))
      }
    }
  }
  message("Added all jobs for ", dataset_name)
}

downsampling <- gather_downsampling_results(downsampling_jobs)
table(sapply(downsampling$dependencies, job_status))
# Store jobs and run them
saveRDS(downsampling, file.path("output/job_storage", "downsampling_job.RDS"))
downsampling <- run_job(downsampling)
saveRDS(downsampling_plotting_jobs, file.path("output/job_storage/downsampling_plotting_jobs.RDS"))
for(job in downsampling_plotting_jobs) run_job(job)


########################################
# Run Downsampling Betst-of benchmarks #
########################################

downsampling_best_of_jobs <- list()
# Loop over all datasets, seeds, transformations, knn and pca settings
for(dataset_name in di2$input_data$dataset){
  for(seed in di2$input_data$seed){
    dataset <- get_data_for_downsampling_benchmark(dataset_name, seed)
    for(pca_dim in di2$knn_construction$pca){
      for(knn in di2$knn_construction$knn){
        full_knn_transformations <- list()
        reduced_knn_transformations <- list()
        for(trans in di2$knn_construction$transformations){
          stopifnot(length(trans$alpha) == 1)
          knn_full <- transform_deeply_sequenced_data(trans$name, dataset, data_mode = "full", pca_dim = pca_dim, knn = knn, alpha = trans$alpha)
          knn_reduced <- transform_deeply_sequenced_data(trans$name, dataset, data_mode = "reduced", pca_dim = pca_dim, knn = knn, alpha = trans$alpha)
          full_knn_transformations <- append(full_knn_transformations, list(knn_full))
          reduced_knn_transformations <- append(reduced_knn_transformations, list(knn_reduced))
        }
        res <- calculate_downsampling_agreement(full_knn_transformations, reduced_knn_transformations,
                                                dataset = dataset_name, seed = seed, pca_dim = pca_dim, knn = knn,
                                                alphas = vapply(di2$knn_construction$transformations, function(e) e$alpha, character(1L)),
                                                transformations = vapply(di2$knn_construction$transformations, function(e) e$name, character(1L)))
        downsampling_best_of_jobs <- append(downsampling_best_of_jobs, list(res))
      }
    }
  }
  message("Added all jobs for ", dataset_name)
}

downsampling_best_of <- gather_downsampling_results(downsampling_best_of_jobs)
table(sapply(downsampling_best_of$dependencies, job_status))
# Store jobs and run them
saveRDS(downsampling_best_of, file.path("output/job_storage", "downsampling_best_of_job.RDS"))
downsampling_best_of <- run_job(downsampling_best_of)



### Start Dataset summarization job ###
dataset_summary <- collect_dataset_summary_statistics(names = map_chr(datasets_to_summarize, "name"),
                                                      dataset = map(datasets_to_summarize, "dataset"))
saveRDS(dataset_summary, file.path("output/job_storage", "dataset_summary.RDS"))
dataset_summary <- run_job(dataset_summary)

#################################################
#  Run Simulation for stratification benchmarks #
#################################################

simulation_strat_jobs <- list()
# Loop over all simulations, seeds, transformations, alphas, knn and pca settings
for(simulator in si2$input_data$simulator){
  for(seed in si2$input_data$seed){
    dataset <- get_data_for_simulation_benchmark(simulator, seed)
    for(pca_dim in si2$knn_construction$pca){
      for(knn in si2$knn_construction$knn){
        transformations <- list()
        for(trans in si2$knn_construction$transformations){
          for(a in trans$alpha){
            knn_constr <- transform_simulated_data(trans$name, dataset, pca_dim = pca_dim, knn = knn, alpha = a)
            transformed_data <- list(name = paste0(trans$name, "-", a), id = knn_constr$result_id, dep = knn_constr)
            transformations <- append(transformations, list(transformed_data))
          }
        }
        simulation_strat_job <- calculate_ground_truth_overlap_for_stratification(dataset, transformations,
                                                                                   simulator = simulator, seed = seed, pca_dim = pca_dim, knn = knn)
        simulation_strat_jobs <- append(simulation_strat_jobs, list(simulation_strat_job))
      }
    }
  }
  message("Added all jobs for ", simulator)
}

table(sapply(simulation_strat_jobs, job_status))
# Store jobs and run them
saveRDS(simulation_strat_jobs, file.path("output/job_storage/simulation_for_stratification_jobs.RDS"))
for(job in simulation_strat_jobs) run_job(job, "normal")

#####################################################
# Run 10X Consistency for stratification benchmarks #
#####################################################

consistency_strat_jobs <- list()
# Loop over all datasets, transformations, knn and pca settings
for(dataset_name in ci2$input_data$dataset){
  for(seed in ci2$input_data$seed){
    dataset <- get_10X_data_for_consistency_benchmark(dataset_name)
    for(pca_dim in ci2$knn_construction$pca){
      for(knn in ci2$knn_construction$knn){
        transformations <- list()
        for(trans in ci2$knn_construction$transformations){
          for(a in trans$alpha){
            knn_constr <- transform_consistency_data(trans$name, dataset, pca_dim = pca_dim, knn = knn, alpha = a, seed = seed)
            transformed_data <- list(name = paste0(trans$name, "-", a), id = knn_constr$result_id, dep = knn_constr)
            transformations <- append(transformations, list(transformed_data))
          }
        }
        consistency_strat_job <- calculate_10X_consistency_for_stratification(dataset, transformations, seed = seed, pca_dim = pca_dim, knn = knn, dataset_name = dataset_name)
        consistency_strat_jobs <- append(consistency_strat_jobs, list(consistency_strat_job))
      }
    }
  }
  message("Added all jobs for ", dataset_name)
}

table(sapply(consistency_strat_jobs, job_status))
# Store jobs and run them
saveRDS(consistency_strat_jobs, file.path("output/job_storage/consistency_for_stratification_jobs.RDS"))
for(job in consistency_strat_jobs) run_job(job, "high")

#################################################
# Run Downsampling for stratification benchmarks#
#################################################

downsampling_strat_jobs <- list()
# Loop over all datasets, seeds, transformations, knn and pca settings
for(dataset_name in di3$input_data$dataset){
  for(seed in di3$input_data$seed){
    dataset <- get_data_for_downsampling_benchmark(dataset_name, seed)
    for(pca_dim in di3$knn_construction$pca){
      for(knn in di3$knn_construction$knn){
        full_knn_transformations <- list()
        reduced_knn_transformations <- list()
        for(trans in di3$knn_construction$transformations){
          stopifnot(length(trans$alpha) == 1)
          knn_full <- transform_deeply_sequenced_data(trans$name, dataset, data_mode = "full", pca_dim = pca_dim, knn = knn, alpha = trans$alpha)
          knn_reduced <- transform_deeply_sequenced_data(trans$name, dataset, data_mode = "reduced", pca_dim = pca_dim, knn = knn, alpha = trans$alpha)
          full_knn_transformations <- append(full_knn_transformations, list(knn_full))
          reduced_knn_transformations <- append(reduced_knn_transformations, list(knn_reduced))
        }
        res <- calculate_downsampling_agreement_for_stratification(full_knn_transformations, reduced_knn_transformations,
                                                dataset = dataset_name, seed = seed, pca_dim = pca_dim, knn = knn,
                                                alphas = vapply(di3$knn_construction$transformations, function(e) e$alpha, character(1L)),
                                                transformations = vapply(di3$knn_construction$transformations, function(e) e$name, character(1L)))
        downsampling_strat_jobs <- append(downsampling_strat_jobs, list(res))
      }
    }
  }
  message("Added all jobs for ", dataset_name)
}

table(sapply(downsampling_strat_jobs, job_status))
# Store jobs and run them
saveRDS(downsampling_strat_jobs, file.path("output/job_storage/downsampling_for_stratification_jobs.RDS"))
for(job in downsampling_strat_jobs) run_job(job, "normal")

