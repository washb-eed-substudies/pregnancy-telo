


rm(list=ls())

source(here::here("0-config.R"))


d<-read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/washb-bd-pregnancy-serum-micronutrient-immun-cortisol-covariates-telo.csv"))
names(d)

summary(d$mom_t0_ln_il2)

######### Add combined ratios of cytokines ##########

#drop Z-score, sd, and ratio measures
d <- d[,!(grepl("^(z_)",colnames(d)) | grepl("^(sd_)",colnames(d)))]

x=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0")[1]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))
x=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0")[2]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))
x=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0")[3]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))


#function to create composite score
create_score <- function(d, numerator_vars=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_pro_il10"){
  for(i in numerator_vars){
    if(i==numerator_vars[1]){
      x = as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }else{
      x = x + as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }
  }
  summary(x)
  
  
  for(i in denominator_vars){
    if(i==denominator_vars[1]){
      y = as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }else{
      y = y + as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }
  }
  summary(y)
  
  score=log(x/y)
  summary(score)
  d$score <- score
  colnames(d)[ncol(d)] <- varname
  return(d)
}

# *Pro-inflammatory cytokines / IL-10
# *(IL-1 + IL-6 + TNF-a) / IL-10
d <- create_score(d, numerator_vars=c("il1_mom_t0", "il6_mom_t0", "tnfa_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_pro_il10")
summary(d$mom_t0_ratio_pro_il10)
ggplot(d, aes(x=mom_t0_ratio_pro_il10)) + geom_density()

# *Th1 / IL-10
# *(IL-12 + IFN) / IL-10
# gen mom_t0_ratio_th1_il10 = (il12_mom_t0 + ifng_mom_t0) / il10_mom_t0
d <- create_score(d, numerator_vars=c("il12_mom_t0", "ifng_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_th1_il10")
summary(d$mom_t0_ratio_th1_il10)
ggplot(d, aes(x=mom_t0_ratio_th1_il10)) + geom_density()

# *Th2 / IL-10 
# *(IL-4 + IL-5 + IL-13) / IL-10
# gen mom_t0_ratio_th2_il10 = (il4_mom_t0 + il5_mom_t0 + il13_mom_t0) / il10_mom_t0
d <- create_score(d, numerator_vars=c("il4_mom_t0", "il5_mom_t0", "il13_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_th2_il10")
summary(d$mom_t0_ratio_th2_il10)
ggplot(d, aes(x=mom_t0_ratio_th2_il10)) + geom_density()


# *Th17 / IL-10
# *(IL-17A + IL-21) / IL-10
# gen mom_t0_ratio_th17_il10 = (il17_mom_t0 + il21_mom_t0) / il10_mom_t0
d <- create_score(d, numerator_vars=c("il17_mom_t0", "il21_mom_t0"), denominator_vars="il10_mom_t0", varname="mom_t0_ratio_th17_il10")
summary(d$mom_t0_ratio_th17_il10)
ggplot(d, aes(x=mom_t0_ratio_th17_il10)) + geom_density()


# *Th1 / Th2
# *(IL-12 + IFN) / (IL-4 + IL-5 + IL-13)
# gen mom_t0_ratio_th1_th2 = (il12_mom_t0 + ifng_mom_t0) / (il4_mom_t0 + il5_mom_t0 + il13_mom_t0)
d <- create_score(d, numerator_vars=c("il12_mom_t0", "ifng_mom_t0"), denominator_vars=c("il4_mom_t0", "il5_mom_t0", "il13_mom_t0"), varname="mom_t0_ratio_th1_th2")
summary(d$mom_t0_ratio_th1_th2)
ggplot(d, aes(x=mom_t0_ratio_th1_th2)) + geom_density()


# *Th1 / Th17
# *(IL-12+IFN) / (IL-17A + IL-21)
# gen mom_t0_ratio_th1_th17 = (il12_mom_t0 + ifng_mom_t0) / (il17_mom_t0 + il21_mom_t0)
d <- create_score(d, numerator_vars=c("il12_mom_t0", "ifng_mom_t0"), denominator_vars=c("il17_mom_t0", "il21_mom_t0"), varname="mom_t0_ratio_th1_th17")
summary(d$mom_t0_ratio_th1_th17)
ggplot(d, aes(x=mom_t0_ratio_th1_th17)) + geom_density()


########### Z-score telomere length ############

d <- d %>% 
  mutate(TS_t2_Z = scale(TS_t2, center=TRUE, scale=TRUE)[,1]) %>%
  mutate(TS_t3_Z = scale(TS_t3, center=TRUE, scale=TRUE)[,1]) %>%
  mutate(delta_TS_Z = scale(delta_TS, center=TRUE, scale=TRUE)[,1])



############## Merge in hhwealth ##################
d_hhwealth <- read.csv("C:/Users/Sophia/Documents/ee-secondary/sophia scripts/hhwealth.csv")
d1 <- left_join(d, d_hhwealth, by="dataid")


############## Merge in maternal sum scores ####################
d_sum <- read.csv("C:/Users/Sophia/Documents/immune-growth/results/maternal sum score/maternal inflammation sum score.csv") %>%
  select(-X)
dfull <- left_join(d1, d_sum, by="dataid")


############# Check covariate missingness ###################



saveRDS(dfull, paste0(dropboxDir,"Data/Cleaned/Audrie/pregnancy_telo_covariates_data.RDS"))


