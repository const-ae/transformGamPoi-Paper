# Online Supplement of _Comparison of Transformations for Single-Cell RNA-Seq Data_ 

[![DOI](https://zenodo.org/badge/378103949.svg)](https://zenodo.org/badge/latestdoi/378103949)


This repository contains the code to generate the figures of our paper 'Comparison of Transformations for Single-Cell RNA-Seq Data' in Nature Methods https://www.nature.com/articles/s41592-023-01814-1. 

The accompanying R package `transformGamPoi` is available via [Bioconductor](https://bioconductor.org/packages/transformGamPoi/) or at https://github.com/const-ae/transformGamPoi. 

At https://shiny-portal.embl.de/shinyapps/app/08_single-cell_transformation_benchmark, we provide an interactive website to explore the benchmark results. The underlying source code is available at https://github.com/const-ae/transformGamPoi-ShinyApp.


# Repository content

The repository is organized into four folders. The root directory contains some additional files for documentation and data plumbing purposes.

* **benchmark** contains the infrastructure for our benchmark of 22 different transformation methods. The [job_overview.yaml](https://github.com/const-ae/transformGamPoi-Paper/blob/master/benchmark/job_overview.yaml) is a declarative specification of all the parameter combinations (more than 58,000) that we benchmarked. The [renv.lock](https://github.com/const-ae/transformGamPoi-Paper/blob/master/benchmark/renv.lock) file contains detailed information about each R package that was used in the benchmark using the [`renv`](https://rstudio.github.io/renv/articles/renv.html) system.
* **notebooks** contains the Rmd scripts to generate the plots shown in the manuscript. Links to the rendered HTML notebooks are below
* **illustrations** and **output** contain the figures

# Code to generate figures

* [Figure 1](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_concept_overview_figure.html): Conceptual differences between variance-stabilizing transformations
* [Figure 2](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_benchmark_results.html): Benchmark results
* [Figure 3](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_benchmark_results.html): Computational expense
* [Figure 4](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_benchmark_contrasts.html): Comparison of selected transformations
--------
* [plot_size_factor_effect_on_dim_reduction.Rmd](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_size_factor_effect_on_dim_reduction.html) generates Extended Data Fig. S1.
* [plot_variance_stabilization_results.Rmd](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_variance_stabilization_results.html) generates Extended Data Fig. S2.
* [plot_bimodal_gene_histograms.Rmd](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_bimodal_gene_histograms.html) generates Extended Data Fig. S3.
* [plot_dataset_structure.Rmd](https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_dataset_structure.html) generates Suppl. Figs. S1, Extended Data Fig. S5, and Suppl. Tabl. S1.
* [plot_benchmark_results.Rmd](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_benchmark_results.html) generates Extended Data Figs. S6, S7, in addition to Fig. 2 and 3.
* [plot_deep_sequenced_overlap.Rmd](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_deep_sequenced_overlap.html) generates Extended Data Fig. S8.
* [plot_clustering_results.Rmd](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_clustering_results.html) generates Extended Data Fig. S9.
* [plot_stratified_performance.Rmd](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/notebooks/plot_stratified_performance.html) generates Extended Data Figs. S10.
