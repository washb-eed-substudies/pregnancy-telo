



rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS("/Users/farheenjamshed/Downloads/bangladesh-cleaned-master-data (2).RDS")
d <- d %>% filter(pregnancy_telo==1)
#d<-readRDS(paste0(dropboxDir, "Data/Cleaned/Audrie/pregnancy_telo_covariates_data.RDS"))
d <- d %>% mutate(vit_A_def = ifelse(RBP_inf_preg < 0.7, 1, 0),
                  vit_A_low = ifelse(RBP_inf_preg < 1.05, 1, 0))

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu","gest_age_weeks", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "HHwealth_scaled",
         "tr", "life_viol_any_t3", "viol_any_preg")

Wvars[!(Wvars %in% colnames(d))]



#Add in time varying covariates:
Wvars2 <- c(Wvars, c("ageday_ht2", "month_blood_t0", "month_ht2"))
Wvars3 <- c(Wvars, c("ageday_ht3", "month_blood_t0", "month_ht3"))
Wvars23 <- c(Wvars, c("ageday_ht2", "ageday_ht3", "month_blood_t0", "month_ht2", "month_ht3"))

pick_covariates <- function(j){
  # j is outcome as string
  # choose correct adjustment set based on outcome
  if(grepl("t2", j)){Wset = Wvars2}
  else if(grepl("t3", j)){Wset = Wvars3}
  else{Wset = Wvars23}
  return(Wset)
}




#Loop over exposure-outcome pairs

##Hypothesis 1
#Maternal nutrition is positively associated with child telomere length
Xvars <- c("vitD_nmol_per_L", "logFERR_inf", "logSTFR_inf", "logRBP_inf", 
           "vit_A_def", "vit_A_low", "iron_def", "vit_D_def")          
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

#Fit models
H1_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}



#Get primary contrasts
H1_adj_res <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  if(grepl("_def|_low", H1_adj_models$X[i])){
    preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y, binary=T)
  }else{
    preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  }
  H1_adj_res <-  bind_rows(H1_adj_res , preds$res)
}

#Make list of plots
H1_adj_plot_list <- NULL
H1_adj_plot_data <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_models$fit[i][[1]], H1_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_plot_list[[i]] <-  simul_plot$p
  H1_adj_plot_data <-  rbind(H1_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H1_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/H1_adj_models.RDS"))

#Save results
saveRDS(H1_adj_res, here("results/adjusted/H1_adj_res.RDS"))


#Save plots
#saveRDS(H1_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H1_adj_splines.RDS"))

#Save plot data
saveRDS(H1_adj_plot_data, here("figure-data/H1_adj_spline_data.RDS"))


## Hypothesis 2
# Maternal stress is negatively associated with child telomere length and postively correlated with
# change in telomere length

Xvars <- c("ln_preg_cort")            
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

#Fit models
H2_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-c(pick_covariates(j), "time_of_day_cort_cont")
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset, forcedW = c("time_of_day_cort_cont"))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_models <- bind_rows(H2_adj_models, res)
  }
}

#Get primary contrasts
H2_adj_res <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H2_adj_models$fit[i][[1]], d=H2_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H2_adj_res <-  bind_rows(H2_adj_res , preds$res)
}

#Make list of plots
H2_adj_plot_list <- NULL
H2_adj_plot_data <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H2_adj_models$fit[i][[1]], H2_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_adj_plot_list[[i]] <-  simul_plot$p
  H2_adj_plot_data <-  rbind(H2_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H2_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_H2_adj_models.RDS"))

#Save results
saveRDS(H2_adj_res, here("results/adjusted/H2_adj_res.RDS"))


#Save plots
#saveRDS(H2_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H2_adj_adj_splines.RDS"))

#Save plot data
saveRDS(H2_adj_plot_data, here("figure-data/H2_adj_spline_data.RDS"))



##Hypothesis 3

Xvars <- c("logCRP", "logAGP", "mom_t0_ln_ifn", "sumscore_t0_mom_Z")            
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

#Fit models
H3_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_models <- bind_rows(H3_models, res)
  }
}

#Get primary contrasts
H3_res <- NULL
for(i in 1:nrow(H3_models)){
  res <- data.frame(X=H3_models$X[i], Y=H3_models$Y[i])
  preds <- predict_gam_diff(fit=H3_models$fit[i][[1]], d=H3_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H3_res <-  bind_rows(H3_res , preds$res)
}

#Make list of plots
H3_plot_list <- NULL
H3_plot_data <- NULL
for(i in 1:nrow(H3_models)){
  res <- data.frame(X=H3_models$X[i], Y=H3_models$Y[i])
  simul_plot <- gam_simul_CI(H3_models$fit[i][[1]], H3_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_plot_list[[i]] <-  simul_plot$p
  H3_plot_data <-  rbind(H3_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H3_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_H3_models.RDS"))

#Save results
saveRDS(H3_res, here("results/adjusted/H3_adj_res.RDS"))


#Save plots
#saveRDS(H3_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H3_adj_splines.RDS"))

#Save plot data
saveRDS(H3_plot_data, here("figure-data/H3_adj_spline_data.RDS"))


##Hypothesis 4
# Maternal estriol is positively associated with child telomere length

Xvars <- c("ln_preg_estri")            
Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")


#Fit models
H4_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset<-pick_covariates(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H4_models <- bind_rows(H4_models, res)
  }
}

#Get primary contrasts
H4_res <- NULL
for(i in 1:nrow(H4_models)){
  res <- data.frame(X=H4_models$X[i], Y=H4_models$Y[i])
  preds <- predict_gam_diff(fit=H4_models$fit[i][[1]], d=H4_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H4_res <-  bind_rows(H4_res , preds$res)
}

#Make list of plots
H4_plot_list <- NULL
H4_plot_data <- NULL
for(i in 1:nrow(H4_models)){
  res <- data.frame(X=H4_models$X[i], Y=H4_models$Y[i])
  simul_plot <- gam_simul_CI(H4_models$fit[i][[1]], H4_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H4_plot_list[[i]] <-  simul_plot$p
  H4_plot_data <-  rbind(H4_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H4_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_H4_models.RDS"))

#Save results
saveRDS(H4_res, here("results/adjusted/H4_adj_res.RDS"))


#Save plots
#saveRDS(H4_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H4_adj_splines.RDS"))

#Save plot data
saveRDS(H4_plot_data, here("figure-data/H4_adj_spline_data.RDS"))

