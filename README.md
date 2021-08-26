# transformGamPoi-Paper

This repository contains the code to generate the figures of our paper 'Transformation and Preprocessing of Single-Cell RNA-Seq Data' on bioRxiv https://www.biorxiv.org/content/10.1101/2021.06.24.449781. 

The accompanying R package `transformGamPoi` is at https://github.com/const-ae/transformGamPoi.

Each figure was created using an Rmarkdown notebook. To see the respective code and output go to
* [Figure 1](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/delta_method_vst_graph.nb.html): Example of heteroskedasticty
* [Figure from _explainer box_](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/delta_method_schematic.nb.html)
* [Figure 2](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/marker_gene_histograms.nb.html): Distribution after transformation for three marker genes
* [Figure 3, S5, S6](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/benchmark_results.nb.html): Benchmark on simulated data (results, duration, and PCA dimension dependence)
* [Figure 4](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/benchmark_deep_seq_data.nb.html): Benchmark on real-world data
* [Figure S1](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/mean_variance_relation.nb.html): Log-log mean variance plots for seven different single-cell datasets
* [Figure S2](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/size_factor_confounding_effect.nb.html): Effect of varying size factor on dimension reduction
* [Figure S3](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/variance_after_transformation.nb.html): Variance stabilization plots
* [Figure S4](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/randomized_quantile_residuals.nb.html): Schematic how the randomized quantile residuals are constructed
* [Figure S7](https://htmlpreview.github.io/?https://github.com/const-ae/transformGamPoi-Paper/blob/master/scripts/describe_deep_seq_data.nb.html): Upset plot of the overlap of neighbors across transformations on deeply sequenced data


