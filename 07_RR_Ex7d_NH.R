library(glmnet)
source("/net/xenon/climphys_backedup/maegli/ET_Adj/Scripts/ET_sensitivity/functions.R")

load("/net/xenon/climphys_backedup/maegli/extET/data/CMIP6/Ex7d/Ex7d_hist585_2d50_XAX.RData")

load("/net/xenon/climphys_backedup/maegli/extET/data/obs/GLEAM/Ex7d_GLEAM_XAX_2d50_anom.RData")
load("/net/xenon/climphys_backedup/maegli/extET/data/obs/X-BASE/Ex7d_XBase_XAX_2d50_anom.RData")
load("/net/xenon/climphys_backedup/maegli/extET/data/obs/ERA5/Ex7d_ERA5Land_XAX_2d50_anom.RData")

Ex7d_hist585_2d50_XAX$X <- Ex7d_hist585_2d50_XAX$X *  (3600*24 /(2.5*10^6))

load("/net/xenon/climphys/maegli/data/land_mask.RData")

land_mask <- intersect(which(!is.na(Ex7d_GLEAM_XAX_2d50_anom$X[1,])),
                       which(!is.na(Ex7d_XBase_XAX_2d50_anom$X[1,]))
) %>% 
  intersect(which(!is.na(Ex7d_ERA5Land_XAX_2d50_anom$X[1,])))

lonlatDT <- get_lonlatDT()
lonlatDT[, land := F]
lonlatDT[land_mask, land := T]
lonlatDT[, ix := .I]
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

#=======NH Mean==========


NH_mask <- lonlatDT[lat > 30][lat< 60][land == T]$ix
area_weights <- get_area_weights()


area_weights_NH <- area_weights[NH_mask] * (1/sum(area_weights[NH_mask]))
NH_mean_Ex7d <- Ex7d_hist585_2d50_XAX$X[, NH_mask] %*% area_weights_NH

Ex7d_FR_table <- data.table(NH_mean = NH_mean_Ex7d[,1],
                            year = Ex7d_hist585_2d50_XAX$M$year,
                            mod = Ex7d_hist585_2d50_XAX$M$mod)
Ex7d_FR_table[, FR := mean(NH_mean), .(mod, year)]



CMIP6_sel <-  c("ACCESS-CM2", "ACCESS-ESM1-5", "CESM2",
                "CESM2-WACCM", "CanESM5", "EC-Earth3", 
                "EC-Earth3-Veg","HadGEM3-GC31-LL", "HadGEM3-GC31-MM", 
                "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", 
                "MPI-ESM1-2-HR", "MPI-ESM1-2-LR","MRI-ESM2-0", "NorESM2-LM","NorESM2-MM", "UKESM1-0-LL")


folds <- data.table(mod = unique(Ex7d_FR_table$mod),
           fold = c("UKESM", "UKESM", "BCC", "CAMS", "CanESM", "CESM", "CESM", "CMCC", "CMCC", "FGOALS", "UKESM", "UKESM",
                    "IITM", "MIROC", "MIROC", "MPI", "MPI", "MRI", "NESM", "CESM", "CESM", "CESM", "UKESM"))


Ex7d_FR_table <- Ex7d_FR_table[folds, on ="mod"]
Ex7d_FR_table[, ix := .I]
Ex7d_FR_table[,fold_id := as.numeric(as.factor(fold))]



Ex7d_FR_table_subset <- Ex7d_FR_table[year %in% 1950:2100][mod %in% CMIP6_sel]
weights <- Ex7d_FR_table_subset[, .N, fold][, weights := 1/(N/sum(N))]
weights[, weights_scaled := weights/sum(weights)]
Ex7d_FR_table_subset <- Ex7d_FR_table_subset[weights, on = "fold"]
Ex7d_FR_table_subset[,fold_id := as.numeric(as.factor(fold))]

cv_fit_NH_Ex7d <- cv.glmnet(x = Ex7d_hist585_2d50_XAX$X[Ex7d_FR_table_subset$ix, NH_mask],
                            y = Ex7d_FR_table_subset$FR,
                            foldid = Ex7d_FR_table_subset$fold_id,
                            keep = T,
                            alpha = 0,
                            weights = Ex7d_FR_table_subset$weights_scaled)

lambda_1se_ix <-  which(cv_fit_NH_Ex7d$lambda == cv_fit_NH_Ex7d$lambda.1se)
lambda_min_ix <-  which(cv_fit_NH_Ex7d$lambda == cv_fit_NH_Ex7d$lambda.min)


eval_table <- data.table( pred_1se = cv_fit_NH_Ex7d$fit.preval[, lambda_1se_ix],
                          pred_min =cv_fit_NH_Ex7d$fit.preval[, lambda_min_ix],
            year = Ex7d_FR_table_subset$year,
            mod = Ex7d_FR_table_subset$mod,
            true = Ex7d_FR_table_subset$FR,
            mem = Ex7d_hist585_2d50_XAX$M$mem[Ex7d_FR_table_subset$ix])

Ex7d_FR_table_plot <- data.table(NH_mean = NH_mean_Ex7d[,1],
                            year = Ex7d_hist585_2d50_XAX$M$year,
                            mod = Ex7d_hist585_2d50_XAX$M$mod
                            )[year %in% 1950:2100, .(FR = mean(NH_mean)), .(mod, year)
                              ][mod %in% CMIP6_sel]


eval_table %>% 
  ggplot()+
  geom_line(aes(year, pred_1se, group = mem), color = "goldenrod", alpha = 0.5) +
  facet_wrap(~mod)+
  geom_line(data = Ex7d_FR_table_plot, aes(year, FR))
    

ggsave(filename = "/net/xenon/climphys_backedup/maegli/extET/figures/20241121_RR/lmo_cv_NH_ts.png",
       height = 8,
       width = 11)

eval_table[, .(mse(pred_1se, true),
               mse(pred_min, true)), .(mod)]

beta_table <- lonlatDT[NH_mask]
beta_table[, beta_1se := cv_fit_NH_Ex7d$glmnet.fit$beta[, lambda_1se_ix]]
beta_table[, beta_min := cv_fit_NH_Ex7d$glmnet.fit$beta[, lambda_min_ix]]

beta_table %>% melt(id.vars = c("lat", "lon", "lonlat", "land", "ix")) %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = value))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick")+
  ylim(30,60)+
  facet_wrap(~variable, ncol = 1)+
  #coord_equal()+
  theme(panel.grid.major = element_blank(), axis.text = element_blank()) 

ggsave(filename = "/net/xenon/climphys_backedup/maegli/extET/figures/20241121_RR/betas_NH.png",
       height = 5,
       width = 11)



obs_FR_table <- rbind(data.table(value = predict(cv_fit_NH_Ex7d, Ex7d_GLEAM_XAX_2d50_anom$X[, NH_mask], s = "lambda.1se")[,1],
           year = Ex7d_GLEAM_XAX_2d50_anom$M$year,
           obs = "GLEAM"),
data.table(value = predict(cv_fit_NH_Ex7d, Ex7d_ERA5Land_XAX_2d50_anom$X[, NH_mask], s = "lambda.1se")[,1],
           year = Ex7d_ERA5Land_XAX_2d50_anom$M$year,
           obs = "ERA5 Land"),
data.table(value = predict(cv_fit_NH_Ex7d, Ex7d_XBase_XAX_2d50_anom$X[, NH_mask], s = "lambda.1se")[,1],
           year = Ex7d_XBase_XAX_2d50_anom$M$year,
           obs = "X-Base")
)

obs_FR_table %>% 
  ggplot()+
  geom_line(aes(year, value, color = obs))+
  ylab("Ex7d anomaly [mm/7day]")

ggsave(filename = "/net/xenon/climphys_backedup/maegli/extET/figures/20241121_RR/ts_obs_NH.png",
       height = 5,
       width = 8)


CMIP6_trends <- eval_table[year %in% 1980:2021, lm(pred_1se~year)$coefficients[2], .(mod, mem)]
obs_trends <- obs_FR_table[year %in% 1980:2021, lm(value~year)$coefficients[2], .(obs)]

CMIP6_trends %>% 
  ggplot()+
  geom_density(aes(V1), fill = "grey", color = NA)+
  geom_vline(data = obs_trends, aes(xintercept = V1, color = obs)) +
  xlab("Trends")

ggsave(filename = "/net/xenon/climphys_backedup/maegli/extET/figures/20241121_RR/trend_density_NH.png",
       height = 5,
       width = 8)


#===========global===============


global_mask <- lonlatDT[land == T]$ix
area_weights <- get_area_weights()


area_weights_global <- area_weights[global_mask] * (1/sum(area_weights[global_mask]))
global_mean_Ex7d <- Ex7d_hist585_2d50_XAX$X[, global_mask] %*% area_weights_global

Ex7d_FR_table_global <- data.table(mean = global_mean_Ex7d[,1],
                            year = Ex7d_hist585_2d50_XAX$M$year,
                            mod = Ex7d_hist585_2d50_XAX$M$mod)
Ex7d_FR_table_global[, FR := mean(mean), .(mod, year)]



CMIP6_sel <-  c("ACCESS-CM2", "ACCESS-ESM1-5", "CESM2",
                "CESM2-WACCM", "CanESM5", "EC-Earth3", 
                "EC-Earth3-Veg","HadGEM3-GC31-LL", "HadGEM3-GC31-MM", 
                "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", 
                "MPI-ESM1-2-HR", "MPI-ESM1-2-LR","MRI-ESM2-0", "NorESM2-LM","NorESM2-MM", "UKESM1-0-LL")


folds <- data.table(mod = unique(Ex7d_FR_table_global$mod),
                    fold = c("UKESM", "UKESM", "BCC", "CAMS", "CanESM", "CESM", "CESM", "CMCC", "CMCC", "FGOALS", "UKESM", "UKESM",
                             "IITM", "MIROC", "MIROC", "MPI", "MPI", "MRI", "NESM", "CESM", "CESM", "CESM", "UKESM"))


Ex7d_FR_table_global <- Ex7d_FR_table_global[folds, on ="mod"]
Ex7d_FR_table_global[, ix := .I]
Ex7d_FR_table_global[,fold_id := as.numeric(as.factor(fold))]



Ex7d_FR_table_global_subset <- Ex7d_FR_table_global[year %in% 1950:2100][mod %in% CMIP6_sel]
weights_global <- Ex7d_FR_table_global_subset[, .N, fold][, weights := 1/(N/sum(N))]
weights_global[, weights_scaled := weights/sum(weights)]
Ex7d_FR_table_global_subset <- Ex7d_FR_table_global_subset[weights_global, on = "fold"]
Ex7d_FR_table_global_subset[,fold_id := as.numeric(as.factor(fold))]

cv_fit_NH_Ex7d_global <- cv.glmnet(x = Ex7d_hist585_2d50_XAX$X[Ex7d_FR_table_global_subset$ix, global_mask],
                            y = Ex7d_FR_table_global_subset$FR,
                            foldid = Ex7d_FR_table_global_subset$fold_id,
                            keep = T,
                            alpha = 0,
                            weights = Ex7d_FR_table_global_subset$weights_scaled)



lambda_1se_ix <-  which(cv_fit_NH_Ex7d_global$lambda == cv_fit_NH_Ex7d_global$lambda.1se)
lambda_min_ix <-  which(cv_fit_NH_Ex7d_global$lambda == cv_fit_NH_Ex7d_global$lambda.min)
eval_table_global <- data.table( pred_1se = cv_fit_NH_Ex7d_global$fit.preval[, lambda_1se_ix],
                          pred_min =cv_fit_NH_Ex7d_global$fit.preval[, lambda_min_ix],
                          year = Ex7d_FR_table_global_subset$year,
                          mod = Ex7d_FR_table_global_subset$mod,
                          true = Ex7d_FR_table_global_subset$FR,
                          mem = Ex7d_hist585_2d50_XAX$M$mem[Ex7d_FR_table_global_subset$ix])

Ex7d_FR_table_plot_global <- data.table(NH_mean = global_mean_Ex7d[,1],
                                 year = Ex7d_hist585_2d50_XAX$M$year,
                                 mod = Ex7d_hist585_2d50_XAX$M$mod
                                 )[year %in% 1950:2100, .(FR = mean(NH_mean)), .(mod, year)
                                   ][mod %in% CMIP6_sel]






eval_table_global %>% 
  ggplot()+
  geom_line(aes(year, pred_1se, group = mem, alpha =0.5), color = "goldenrod") +
  facet_wrap(~mod)+
  geom_line(data = Ex7d_FR_table_plot_global, aes(year, FR))


eval_table_global[, .(mse(pred_1se, true),
                      mse(pred_min, true)), .(mod)]

beta_table_global <- lonlatDT[global_mask]
beta_table_global[, beta_1se := cv_fit_NH_Ex7d_global$glmnet.fit$beta[, lambda_1se_ix]]
beta_table_global[, beta_min := cv_fit_NH_Ex7d_global$glmnet.fit$beta[, lambda_min_ix]]

beta_table_global %>% melt(id.vars = c("lat", "lon", "lonlat", "land", "ix")) %>% 
  ggplot()+
  geom_tile(aes(lon, lat, fill = value))+
  geom_sf(data = world, fill=alpha("lightgrey", 0), color="grey50") +
  ylab("") + xlab("") +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "firebrick")+
  ylim(-70,85)+
  facet_wrap(~variable, ncol = 1)+
  #coord_equal()+
  theme(panel.grid.major = element_blank(), axis.text = element_blank()) 

ggsave(filename = "/net/xenon/climphys_backedup/maegli/extET/figures/20241121_RR/betas_global.png",
       height = 5,
       width = 9)


obs_FR_table_global <- rbind(data.table(value = predict(cv_fit_NH_Ex7d_global, Ex7d_GLEAM_XAX_2d50_anom$X[, global_mask], s = "lambda.min")[,1],
                                 year = Ex7d_GLEAM_XAX_2d50_anom$M$year,
                                 obs = "GLEAM"),
                      data.table(value = predict(cv_fit_NH_Ex7d_global, Ex7d_ERA5Land_XAX_2d50_anom$X[, global_mask], s = "lambda.min")[,1],
                                 year = Ex7d_ERA5Land_XAX_2d50_anom$M$year,
                                 obs = "ERA5 Land"),
                      data.table(value = predict(cv_fit_NH_Ex7d_global, Ex7d_XBase_XAX_2d50_anom$X[, global_mask], s = "lambda.min")[,1],
                                 year = Ex7d_XBase_XAX_2d50_anom$M$year,
                                 obs = "X-Base")
)

obs_FR_table_global %>% 
  ggplot()+
  geom_line(aes(year, value, color = obs))+
  ylab("Ex7d anomaly [mm/7day]")

ggsave(filename = "/net/xenon/climphys_backedup/maegli/extET/figures/20241121_RR/betas_global.png",
       height = 5,
       width = 9)

CMIP6_trends_global <- eval_table_global[year %in% 1980:2021, lm(pred_1se~year)$coefficients[2], .(mod, mem)]
obs_trends_global <- obs_FR_table_global[, lm(value~year)$coefficients[2], .(obs)]

CMIP6_trends_global %>% 
  ggplot()+
  geom_density(aes(V1), fill = "grey", color = NA)+
  geom_vline(data = obs_trends_global, aes(xintercept = V1, color = obs)) +
  xlab("Trends")

ggsave(filename = "/net/xenon/climphys_backedup/maegli/extET/figures/20241121_RR/trend_density_global.png",
       height = 5,
       width = 9)



