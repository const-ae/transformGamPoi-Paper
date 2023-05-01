source("src/transformations/transformation_helper.R")
file = readRDS("realdataset,rds") #require both count matrix and size factor vector

UMI = file@assays$data$counts
sf = file$total_counts
alpha = 0.05


# example method

out = pearson_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/pearson.rds")

out = pearson_analytic_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/pearson_analytic.rds")

out = pearson_clip_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/pearson_clip.rds")

out = pearson_clip_zscore_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/pearson_clip_zscore.rds")

out = pearson_clip_hvg_zscore_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/pearson_clip_hvg_zscore.rds")

out = pearson_clip_hvg_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/pearson_clip_hvg.rds")

out = sctransform_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/sctransform.rds")

out = sctransformp_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/sctransformp.rds")

out = rand_quantile_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/rand_quantile.rds")

out = logp1_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/logp1.rds")

out = logp_cpm_fnc(UMI, sf, alpha)
saveRDS(out, "real_output/logp_cpm.rds")