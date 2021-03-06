
rm(list=ls())

source(here::here("0-config.R"))
source(here::here("src/0-gam-functions.R"))

d<-read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-stress-growth-covariates-stresslab-anthro.csv"))

#Example:

#Fit GAM model with random effects for childid
res_unadj <- fit_RE_gam(d=d, X="t3_cort_z01", Y="laz_t3",  W=NULL)

#Get predictions of differences from the 25th percentile of exposure
preds_unadj <- predict_gam_diff(fit=res_unadj$fit, d=res_unadj$dat, quantile_diff=c(0.25,0.75), Xvar="t3_cort_z01", Yvar="laz_t3")


#Primary parameter we are estimating: difference between 25th and 75th percentile of the exposure
preds_unadj$res

#Plot the difference from the 25th percentile for the full range of the exposure:
#NOTE: not making these plots anymore, just using for diagnostics
p <- plot_gam_diff(preds_unadj$plotdf)
print(p)

#Fit spline with simultaneous confidence intervals
simul_plot <- gam_simul_CI(res_unadj$fit, res_unadj$dat, xlab="t3_cort_z01", ylab="laz_t3", title="unadjusted")
simul_plot$p




#Fit adjusted model
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "elec", "asset_wardrobe",
         "asset_table", "asset_chair", "asset_clock","asset_khat", "asset_chouki", 
         "asset_radio", "asset_tv", "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", 
         "asset_mobile", "n_cattle", "n_goat", "n_chicken", "cesd_sum_t2", "diar7d_t2", "lenhei_med_t2", "weight_med_t2")


#Check if covariates 
Wvars[!(Wvars %in% colnames(d))]


#Fit GAM model -adjusted
res_adj <- fit_RE_gam(d=d, X="t3_cort_z01", Y="laz_t3",  W=Wvars)
preds_adj <- predict_gam_diff(fit=res_adj$fit, d=res_adj$dat, quantile_diff=c(0.25,0.75), Xvar="t3_cort_z01", Yvar="laz_t3")
preds_adj$res
p <- plot_gam_diff(preds_adj$plotdf)
print(p)

simul_plot_adj <- gam_simul_CI(res_adj$fit, res_adj$dat, xlab="t3_cort_z01", ylab="laz_t3", title="adjusted")
simul_plot_adj$p






#Loop over outcomes

Xvars <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i")            
Yvars <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2", 
           "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3",
           "laz_t3", "waz_t3", "whz_t3", "hcz_t3",
           "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3")

#Fit models
H1_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}



#Get primary contrasts
H1_res <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  preds <- predict_gam_diff(fit=H1_models$fit[i][[1]], d=H1_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H1_res <-  bind_rows(H1_res , preds$res)
}
H1_res$adjusted <- 0

#Make list of plots
H1_plot_list <- NULL
H1_plot_data <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  simul_plot <- gam_simul_CI(H1_models$fit[i][[1]], H1_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_plot_list[[i]] <-  simul_plot$p
  H1_plot_data <-  rbind(H1_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}


# #Save models
# saveRDS(H1_models, here("models/H1_models.RDS"))
# 
# #Save results
# saveRDS(H1_res, here("results/unadjusted/H1_res.RDS"))
# 
# 
# #Save plots
# saveRDS(H1_plot_list, here("figure-objects/H1_unadj_splines.RDS"))
# 
# #Save plot data
# saveRDS(H1_plot_data, here("figure-data/H1_unadj_spline_data.RDS"))




#Explore GAM function



res_adj$fit
  res_adj$dat
  xlab="t3_cort_z01"
  ylab="laz_t3"
  title="adjusted"



#fit_RE_gam
  d=d
  X="t3_cort_z01"
  Y="laz_t3"
  W=Wvars
   V=NULL
   id="clusterid"
   family = "gaussian"
   pval = 0.2
   print=TRUE
   
   
#predict_gam_diff
   fit=res_adj$fit
   d=res_adj$dat
   quantile_diff=c(0.25,0.75)
   Xvar="t3_cort_z01"
   Yvar="laz_t3"



