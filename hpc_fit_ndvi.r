#!/usr/bin/env Rscript

load("ndvi_df_hpc.Rdata")


ndvi_df$date2 <- as.numeric(temp_df$date)
ndvi_df$bird <- factor(ndvi_df$bird)

fit_ndvi <- brm(bf(ndvi ~ 0 + season + (1|bird) + ar(time = date2, gr = bird, cov = T), 
                   phi ~ 0 + season, family = Beta()), data = ndvi_df, 
                chains = 3, iter = 1000, cores = getOption("mc.cores", 3),
                control = list(adapt_delta = 0.99))
fit_ndvi

inv_logit_scaled(fixef(fit_ndvi))

save(fit_ndvi, file = "fit_ndvi.Rdata")