library(glmnet)
source("/net/xenon/climphys_backedup/maegli/ET_Adj/Scripts/ET_sensitivity/functions.R")

load("/net/xenon/climphys/maegli/data/XAX/seas/anom/CMIP6_hfls_JJA_hist585_anom_2d50_XAX.RData")
load("/net/xenon/climphys/maegli/data/land_mask.RData")


load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/GLEAM/GLEAM_3.6a_JJA_2d50_anom_XAX.RData")
load("/net/xenon/climphys_backedup/maegli/ET_Adj/Obs/ERA5-Land/ERA5Land_ET_JJA_2d50_anom_XAX.RData")


land_mask

lonlatDT <- get_lonlatDT()
lonlatDT[, land := F]
lonlatDT[land_mask, land := T]
lonlatDT[, ix := .I]
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")




NH_mask <- lonlatDT[lat > 30][lat< 60][land == T]$ix
area_weights <- get_area_weights()


area_weights_NH <- area_weights[NH_mask] * (1/sum(area_weights[NH_mask]))
NH_mean_ET_JJA <- CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[, NH_mask] %*% area_weights_NH

ET_JJA_FR_table <- data.table(NH_mean = NH_mean_ET_JJA[,1],
                            year = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$year,
                            mod = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$mod)
ET_JJA_FR_table[, FR := mean(NH_mean), .(mod, year)]




folds <- data.table(mod = unique(ET_JJA_FR_table$mod),
                    fold = c("UKESM", "UKESM", "CanESM", "CESM", "CNRM", "CNRM", "EC-Earth", "EC-Earth", "FGOALS", "GISS", "GISS", "UKESM",
                             "IPSL", "MIROC", "MIROC", "MPI", "MRI", "CESM", "UKESM"))


ET_JJA_FR_table <- ET_JJA_FR_table[folds, on ="mod"]
ET_JJA_FR_table[, ix := .I]
ET_JJA_FR_table[,fold_id := as.numeric(as.factor(fold))]



ET_JJA_FR_table_subset <- ET_JJA_FR_table[year %in% 1950:2100]
weights <- ET_JJA_FR_table_subset[, .N, fold][, weights := 1/(N/sum(N))]
weights[, weights_scaled := weights/sum(weights)]
ET_JJA_FR_table_subset <- ET_JJA_FR_table_subset[weights, on = "fold"]
ET_JJA_FR_table_subset[,fold_id := as.numeric(as.factor(fold))]

cv_fit_NH_ET_JJA <- cv.glmnet(x = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$X[ET_JJA_FR_table_subset$ix, NH_mask],
                            y = ET_JJA_FR_table_subset$FR,
                            foldid = ET_JJA_FR_table_subset$fold_id,
                            keep = T,
                            alpha = 0,
                            weights = ET_JJA_FR_table_subset$weights_scaled)

lambda_1se_ix <-  which(cv_fit_NH_ET_JJA$lambda == cv_fit_NH_ET_JJA$lambda.1se)
lambda_min_ix <-  which(cv_fit_NH_ET_JJA$lambda == cv_fit_NH_ET_JJA$lambda.min)


eval_table_ET_JJA <- data.table( pred_1se = cv_fit_NH_ET_JJA$fit.preval[, lambda_1se_ix],
                          pred_min =cv_fit_NH_ET_JJA$fit.preval[, lambda_min_ix],
                          year = ET_JJA_FR_table_subset$year,
                          mod = ET_JJA_FR_table_subset$mod,
                          true = ET_JJA_FR_table_subset$FR,
                          mem = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$ens.mem[ET_JJA_FR_table_subset$ix])

ET_JJA_FR_table_plot <- data.table(NH_mean = NH_mean_ET_JJA[,1],
                                 year = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$year,
                                 mod = CMIP6_hfls_JJA_hist585_anom_2d50_XAX$M$mod
)[year %in% 1950:2100, .(FR = mean(NH_mean)), .(mod, year)]

eval_table_ET_JJA %>% 
  ggplot()+
  geom_line(aes(year, pred_1se, group = mem), color = "goldenrod", alpha = 0.5) +
  facet_wrap(~mod)+
  geom_line(data = ET_JJA_FR_table_plot, aes(year, FR))



obs_FR_table_ET_JJA <- rbind(data.table(value = predict(cv_fit_NH_ET_JJA, GLEAM_3.6a_JJA_2d50_anom_XAX$X[, NH_mask], s = "lambda.1se")[,1],
                                 year = GLEAM_3.6a_JJA_2d50_anom_XAX$M$year,
                                 obs = "GLEAM"),
                      data.table(value = predict(cv_fit_NH_ET_JJA, ERA5Land_ET_JJA_2d50_anom_XAX$X[, NH_mask], s = "lambda.1se")[,1],
                                 year = ERA5Land_ET_JJA_2d50_anom_XAX$M$year,
                                 obs = "ERA5 Land")
)

obs_FR_table_ET_JJA %>% 
  ggplot()+
  geom_line(aes(year, value, color = obs))+
  ylab("JJA ET anomaly [mm/day]")



CMIP6_ET_JJA_trends <- eval_table_ET_JJA[year %in% 1980:2021, lm(pred_1se~year)$coefficients[2], .(mod, mem)]
obs_trends_JJA_ET <- obs_FR_table_ET_JJA[year %in% 1980:2021, lm(value~year)$coefficients[2], .(obs)]

CMIP6_ET_JJA_trends %>% 
  ggplot()+
  geom_density(aes(V1), fill = "grey", color = NA)+
  geom_vline(data = obs_trends_JJA_ET, aes(xintercept = V1, color = obs)) +
  xlab("Trends")
