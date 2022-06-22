
trans_families <- list(delta_method = c("logp1", "acosh", "logp_alpha", "logp_cpm", "logp1_size_normed", "logp1_hvg", "logp1_zscore",  "logp1_hvg_zscore"),
                       glm_residual = c("pearson_clip", "sctransform", "pearson_analytic", "rand_quantile", "pearson",  "pearson_clip_hvg", "pearson_clip_zscore", "pearson_clip_hvg_zscore"),
                       latent_expr = c("sanity_map", "sanity_dists", "dino", "normalisr_normvar"),
                       count_model = c("glmpca", "newwave")) %>%
  enframe() %>%
  unnest(value) %>%
  transmute(transformation = value, family = name)


# trans_labels <- expression("acosh" = acosh(2*alpha*y/s+1), "logp_alpha" = log(y/s+1/(4*alpha)), "logp_cpm" = log(CPM+1),
#                            "logp1" = log(y/s+1), "logp1_hvg" = log(y/s+1)%->%" HVG", "logp1_hvg_zscore" = log(y/s+1)%->%" HVG"%->%"Z", 
#                            "logp1_zscore" = log(y/s+1)%->%"Z", "logp1_size_normed" = log(y/s+1)/u,
#                            "pearson" = "Pearson"~"(no"~"clip)", "pearson_analytic" = "Analytic"~"Pearson", "pearson_clip" = "Pearson",
#                            "rand_quantile" = "Random"~"Quantile", "sctransform" = "sctransform",
#                            "pearson_clip_hvg" = "Pearson"%->%"HVG", "pearson_clip_hvg_zscore" = "Pearson"%->%"HVG"%->%"Z", 
#                            "pearson_clip_zscore"= "Pearson"%->%"Z",
#                            "sanity_map"="Sanity"~"MAP", "sanity_dists" = "Sanity"~"Distance", "dino" = "Dino", "normalisr_normvar" = "Normalisr")
trans_labels <- c("acosh" = r"($\textrm{acosh}(2\alpha y/s+1)$)", "logp_alpha" = r"($\log(y/s+1/(4\alpha))$)", "logp_cpm" = r"($\log(\textrm{CPM}+1)$)",
                  "logp1" = r"($\log(y/s+1)$)", "logp1_hvg" = r"($\log(y/s+1)\rightarrow$HVG)", "logp1_hvg_zscore" = r"($\log(y/s+1)\rightarrow$HVG$\rightarrow$Z)", 
                  "logp1_zscore" = r"($\log(y/s+1)\rightarrow$Z)", "logp1_size_normed" = r"($\log(y/s+1)/u$)",
                  "pearson" = "Pearson (no clip)", "pearson_analytic" = "Analytic Pearson", "pearson_clip" = "Pearson",
                  "rand_quantile" = "Random Quantile", "sctransform" = "sctransform",
                  "pearson_clip_hvg" = r"(Pearson$\rightarrow$HVG)", "pearson_clip_hvg_zscore" = r"(Pearson$\rightarrow$HVG$\rightarrow$Z)", 
                  "pearson_clip_zscore"= r"(Pearson$\rightarrow$Z)",
                  "sanity_map"="Sanity MAP", "sanity_dists" = "Sanity Distance", "dino" = "Dino", "normalisr_normvar" = "Normalisr",
                  "glmpca" = "GLM PCA", "newwave" = "NewWave")

trans_labels_plain <- c("acosh" = "acosh(2αy/s+1)", "logp_alpha" = "log(y/s+1/(4α))", "logp_cpm" = "log(CPM + 1)",
                  "logp1" = "log(y/s+1)", "logp1_hvg" = "log(y/s+1)->HVG", "logp1_hvg_zscore" = "log(y/s+1)->HVG->Z", 
                  "logp1_zscore" = "log(y/s+1)->Z", "logp1_size_normed" = "log(y/s+1)/u",
                  "pearson" = "Pearson (no clip)", "pearson_analytic" = "Analytic Pearson", "pearson_clip" = "Pearson",
                  "rand_quantile" = "Random Quantile", "sctransform" = "sctransform",
                  "pearson_clip_hvg" = r"(Pearson->HVG)", "pearson_clip_hvg_zscore" = r"(Pearson->HVG->Z)", 
                  "pearson_clip_zscore"= r"(Pearson->Z)",
                  "sanity_map"="Sanity MAP", "sanity_dists" = "Sanity Distance", "dino" = "Dino", "normalisr_normvar" = "Normalisr",
                  "glmpca" = "GLM PCA", "newwave" = "NewWave")

trans_families$transformation <- factor(trans_families$transformation, levels = trans_families$transformation)

trans_families_labels <- factor(c(delta_method = "Delta Method", glm_residual = "GLM Residuals", latent_expr = "Lat. Expr.", count_model = "Count"), levels = c("Count", "Lat. Expr.", "GLM Residuals", "Delta Method"))
trans_families_labels_long <- factor(c(delta_method = "Delta Method", glm_residual = "GLM Residuals", latent_expr = "Latent Expression", count_model = "Count Model"), levels = c("Count Model", "Latent Expression", "GLM Residuals", "Delta Method"))
trans_families_colors <- c(delta_method = "#66C2A5", glm_residual = "#FC8D62", latent_expr = "#8DA0CB", count_model = "#e78ac3")

dataset_labels <- c(dyngen = "Dyngen", linear_walk = "Linear Walk", muscat = "muscat", random_walk = "Random Walk", scDesign2 = "scDesign2",
                    smartSeq3_fibroblasts = "Fibroblasts (ss3)", smartSeq3_fibroblasts_alt = "Fibroblasts 2 (ss3)", 
                    smartSeq3_hek = "HEK (ss3)", smartSeq3_siRNA_knockdown = "siRNA KD (ss3)", mcSCRB = "mcSCRB",
                    GSE130931 = "GSE130931", GSE142647 = "GSE142647", GSE150068 = "GSE150068", GSE158941 = "GSE158941", GSE163505 = "GSE163505",
                    GSE164017 = "GSE164017", GSE178765 = "GSE178765", GSE179714 = "GSE179714", GSE179831 = "GSE179831", GSE184806 = "GSE184806")

dataset_benchmark <- c(dyngen = "simulation", linear_walk = "simulation", muscat = "simulation", random_walk = "simulation", scDesign2 = "simulation",
                       smartSeq3_fibroblasts = "downsampling", smartSeq3_fibroblasts_alt = "downsampling", 
                       smartSeq3_hek = "downsampling", smartSeq3_siRNA_knockdown = "downsampling", mcSCRB = "downsampling",
                       GSE130931 = "consistency", GSE142647 = "consistency", GSE150068 = "consistency", GSE158941 = "consistency", GSE163505 = "consistency",
                       GSE164017 = "consistency", GSE178765 = "consistency", GSE179714 = "consistency", GSE179831 = "consistency", GSE184806 = "consistency")
