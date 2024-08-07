
#TO DO!
# Need to merge in full dataset from main trial
# Or possibly just dont subset to the pregnancy_telon group so it is a subset of
# The full EED cohort.
# Or should it be subset to just the EED cohort collected arms (C and N+WSH)?


rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS("C:/Users/andre/Documents//EE/eed-substudy-data/bangladesh-cleaned-master-data.RDS")

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

#-------------------------------------------------------------------------
#set ipcw parameters
#-------------------------------------------------------------------------

Xvars <- c("vitD_nmol_per_L", "logFERR_inf", "logSTFR_inf", "logRBP_inf", 
           "vit_A_def", "vit_A_low", "iron_def", "vit_D_def", #H1
           "ln_preg_cort",#H2 
           "logCRP", "logAGP", "mom_t0_ln_ifn", "sumscore_t0_mom_Z", #H3 
           "ln_preg_estri") #H4
Xvars_bin <- c("vit_A_def", "vit_A_low", "iron_def", "vit_D_def")   
table(d$vit_A_def)

Yvars <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

#Create indicators for missingness
d$TS_t2_Z.miss<-ifelse(is.na(d$TS_t2_Z),0,1)
d$TS_t3_Z.miss<-ifelse(is.na(d$TS_t3_Z),0,1)
d$delta_TS_Z.miss<-ifelse(is.na(d$delta_TS_Z),0,1)

table(d$TS_t2_Z.miss)
table(d$TS_t3_Z.miss)
table(d$delta_TS_Z.miss)

# set missing outcomes to an arbitrary, non-missing value. In this case use 9
d$TS_t2_ZDelta <- d$TS_t2_Z
d$TS_t2_ZDelta[d$TS_t2_Z.miss==0] <- (9)

d$TS_t3_ZDelta <- d$TS_t3_Z
d$TS_t3_ZDelta[d$TS_t3_Z.miss==0] <- (9)

d$delta_TS_ZDelta <- d$delta_TS_Z
d$delta_TS_ZDelta[d$delta_TS_Z.miss==0] <- (9)

#percent missing: (to do: first subset to only those with any exposure and any outcome)
prop.table(table(d$TS_t2_Z.miss)) * 100
prop.table(table(d$TS_t3_Z.miss)) * 100
prop.table(table(d$delta_TS_Z.miss)) * 100

#Order for replication:
d<-d[order(d$block,d$clusterid,d$childid),]


#dataframes of urine biomarkers and missingness:
Y<-d %>% select(TS_t2_ZDelta,TS_t3_ZDelta,delta_TS_ZDelta)
miss<-d %>% select(TS_t2_Z.miss,TS_t3_Z.miss,delta_TS_Z.miss)

#Create empty matrix to hold the ipcw results:
res_adj<-list(neo_t1_adj=matrix(0,5,5), mpo_t1_adj=matrix(0,5,5), aat_t1_adj=matrix(0,5,5))
res_adj<-NULL



for(i in 1:ncol(Y)){
  for(j in 1:length(Xvars_bin)){
    
    Wset <- pick_covariates(names(Y)[i])
    W <- d %>% select(Wset)
    #note the log transformation of the outcome prior to running GLM model:
    temp <- washb_tmle_ipcw(Y=Y[,i], Delta=miss[,i], tr=d[[Xvars_bin[j]]], W=W, id=d$block, pair=NULL, family="gaussian", contrast=c(0,1), Q.SL.library = c("SL.glm"), seed=12345, print=T)
    cat(i," : ",j, "\n")
    temp <- (t(unlist(temp$estimates$ATE)))[,1:5]
    temp$Y <- names(Y)[i]
    temp$X <- Xvars_bin[j]
    res_adj <- bind_rows(res_adj, temp)

  }
}



#Save results
saveRDS(res_adj, here("results/adjusted/ipcw_res.RDS"))


#compare to static TMLE
Y<-d %>% select(TS_t2_Z,TS_t3_Z,delta_TS_Z)
res_adj_tmle=NULL
for(i in 1:ncol(Y)){
  for(j in 1:length(Xvars_bin)){
    
    Wset <- pick_covariates(names(Y)[i])
    W <- d %>% select(Wset)
    #note the log transformation of the outcome prior to running GLM model:
    temp <- washb_tmle(Y=Y[,i], tr=d[[Xvars_bin[j]]], W=W, id=d$block, pair=NULL, family="gaussian", contrast=c(0,1), Q.SL.library = c("SL.glm"), seed=12345, print=T)
    cat(i," : ",j, "\n")
    temp <- (t(unlist(temp$estimates$ATE)))[,1:5]
    temp$Y <- names(Y)[i]
    temp$X <- Xvars_bin[j]
    res_adj_tmle <- bind_rows(res_adj_tmle, temp)
    
  }
}
res_adj_tmle



#-------------------------------------------------------------------------
# make comparison plot with primary results
#-------------------------------------------------------------------------


gam_res <- readRDS(here("results/adjusted/H1_adj_res.RDS"))

res_adj <- res_adj %>% rename(ATE = psi) %>% mutate(Y = gsub("Delta", "", Y))
gam_res <- gam_res %>% rename(ATE = point.diff, CI1=lb.diff, CI2=ub.diff)  


plot_df <- bind_rows(res_adj %>% mutate(Analysis = "IPCW"), 
          gam_res %>% mutate(Analysis = "Primary")) %>%
  filter(!is.na(ATE)) %>% arrange(Y, X) %>%
  group_by(Y, X) %>% filter(n()==2)
head(plot_df)

head(res_adj)
head(gam_res)

table(plot_df$analysis)

plot_df <- plot_df %>% mutate(
  Y_f = factor(Y, levels=c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z"), labels=c("Telomeres Year 1", "Telomeres Year 2", "Telomere Change\nYear 1 and Year 2")),
  X_f = factor(X, levels=c("vit_A_def", "vit_A_low", "iron_def", "vit_D_def"), labels=c("Vit A Deficiency","Low Vit A", "Iron Deficiency", "Vitamin D Deficiency"))
)

p <- ggplot(plot_df, aes(x=X_f, y=ATE, ymin=CI1, ymax=CI2, color=Analysis)) +
  geom_pointrange(position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~Y_f, scales="free") +
  coord_flip() +
  theme_ki() +
  theme(legend.position="bottom") +
  labs(title="Comparison of IPCW and GAM results", x="Exposure", y="Adjusted mean difference") +
  scale_color_manual(values=tableau10)
p


ggsave(p, filename="figures/ipcw_sensitivity.jpg", width=10, height=7)


#Add in tmle for a full comparison
res_adj_tmle <- res_adj_tmle %>% rename(ATE = psi) %>% mutate(Y = gsub("Delta", "", Y))
plot_df <- bind_rows(res_adj %>% mutate(Analysis = "IPCW"), 
                     res_adj_tmle %>% mutate(Analysis = "TMLE"), 
                     gam_res %>% mutate(Analysis = "Primary")) %>%
  filter(!is.na(ATE)) %>% arrange(Y, X) %>%
  group_by(Y, X) %>% filter(n()==3)
head(plot_df)

head(res_adj)
head(gam_res)

table(plot_df$Analysis)

plot_df <- plot_df %>% mutate(
  Y_f = factor(Y, levels=c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z"), labels=c("Telomeres Year 1", "Telomeres Year 2", "Telomere Change\nYear 1 and Year 2")),
  X_f = factor(X, levels=c("vit_A_def", "vit_A_low", "iron_def", "vit_D_def"), labels=c("Vit A Deficiency","Low Vit A", "Iron Deficiency", "Vitamin D Deficiency"))
)

p <- ggplot(plot_df, aes(x=X_f, y=ATE, ymin=CI1, ymax=CI2, color=Analysis)) +
  geom_pointrange(position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~Y_f, scales="free") +
  coord_flip() +
  theme_ki() +
  theme(legend.position="bottom") +
  labs(title="Comparison of IPCW and GAM results", x="Exposure", y="Adjusted mean difference") +
  scale_color_manual(values=tableau10)
p
