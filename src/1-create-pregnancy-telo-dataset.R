


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


########### Add deficiency cutoff exposures #############
# 1 if deficient, 0 if not deficient
d$vit_A_def <- ifelse(d$RBP_inf_preg < 0.83, 1, 0)
d$vit_D_def <- ifelse(d$vitD_nmol_per_L < 30, 1, 0)
d$iron_def <- ifelse(d$FERR_inf_preg < 12 | d$STFR_inf_preg > 8.3, 1, 0)
  
  
############## Merge in hhwealth ##################

d_hhwealth <- read.csv(paste0(dropboxDir,"Data/Cleaned/Audrie/hhwealth.csv"))
d1 <- left_join(d, d_hhwealth, by="dataid")


############## Merge in maternal sum scores ####################
d_sum <- read.csv(paste0(dropboxDir,"Data/Cleaned/Audrie/maternal inflammation sum score.csv")) %>%
  select(-X)
dfull <- left_join(d1, d_sum, by="dataid")


############# Check covariate missingness ###################
exp <- c("vitD_nmol_per_L", "logFERR_inf", "logSTFR_inf", "logRBP_inf",
         "ln_preg_cort", "logCRP", "logAGP", "ifng_mom_t0", "sumscore_t0_mom_Z", "ln_preg_estri")
out <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")

Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "roof", "HHwealth",
         "tr", "life_viol_any_t3", "viol_any_preg", "ageday_ht2", "ageday_ht3", 
         "month_blood_t0", "month_ht2", "month_ht3") %>% unique()

W <- dfull %>% select(all_of(Wvars))  

miss <- data.frame(name = names(W), missing = colSums(is.na(W))/nrow(W), row.names = c(1:ncol(W)))
for (i in 1:nrow(miss)) {
  miss$class[i] <- class(W[,which(colnames(W) == miss[i, 1])])
}
miss

for (w in Wvars){
  print(w)
  if(is.numeric(W[,w])){
    print(summary(W[,w]))
  }
  else{print(table(W[,w]))}
}

# roof has low variability
mean(W$roof, na.rm=T)
sd(W$roof, na.rm=T)
# remove roof from covariates

# add missingness category to IPV covariates
dfull$life_viol_any_t3<-as.factor(dfull$life_viol_any_t3)
dfull$life_viol_any_t3<-addNA(dfull$life_viol_any_t3)
levels(dfull$life_viol_any_t3)[length(levels(dfull$life_viol_any_t3))]<-"Missing"

dfull$viol_any_preg<-as.factor(dfull$viol_any_preg)
dfull$viol_any_preg<-addNA(dfull$viol_any_preg)
levels(dfull$viol_any_preg)[length(levels(dfull$viol_any_preg))]<-"Missing"

for(outcome in out){
  d_sub <- subset(dfull, !is.na(dfull[,outcome]))
  W_sub <- d_sub %>% select(all_of(Wvars))  
  
  miss_sub <- data.frame(name = names(W_sub), missing = colSums(is.na(W_sub)), row.names = c(1:ncol(W_sub)))
  for (i in 1:nrow(miss_sub)) {
    miss_sub$class[i] <- class(W_sub[,which(colnames(W_sub) == miss_sub[i, 1])])
  }
  print(outcome)
  print(miss_sub)
}



#IPV raw survey for frequencies
ipv_sv <-  haven::read_dta("C:/Users/andre/Dropbox/WASHB-EE-analysis/WBB-EE-analysis/Data/Untouched/MHLE_Mother_clean_data_30May16_deidentified.dta") 
head(ipv_sv)
labels <- sapply(ipv_sv, function(x) attr(x, "label"))
labels <- tibble(name = names(labels),
                 label = labels)
labels$label<- as.character(labels$label)
labels[grepl("ow often",labels$label),]

# NEVER .................................................... 1
# ONCE .......................................................... 2
# SEVERAL TIMES ................................... 3
# MANY TIMES/MOST OF THE TIME.. 4
# DON'T KNOW ................. 8
# REFUSED/NO ANSWER ..... 9

labels[labels$name=="q_1002",]
table(ipv_sv$q_1002)

labels[labels$name=="q_808",]
table(ipv_sv$q_808)

labels[labels$name=="q_902_a",]
table(ipv_sv$q_902_a)

#Has your current/ most recent husband/partner, ever..
labels[labels$name=="q_805_a",]
table(ipv_sv$q_805_a)

labels[labels$name=="q_805_b",]
table(ipv_sv$q_805_b)

labels[labels$name=="q_805_c",]
table(ipv_sv$q_805_c)

labels[labels$name=="q_805_d",]
table(ipv_sv$q_805_d)
#Note! Not enough variation for these variables to be used


#clean and merge into the main dataset
ipv_sv <- ipv_sv %>% select(dataid, q_1002, q_808, q_902_a) %>%
  rename(ipv_child = q_1002, ipv_afraid = q_808, ipv_life_freq = q_902_a) %>%
  mutate(ipv_child=case_when(ipv_child == 1 ~ "Never", ipv_child == 2 ~ "Once", ipv_child == 3 ~ "Several times", ipv_child == 4 ~"Most/all of the times" , TRUE ~ "Missing"),
         ipv_afraid=case_when(ipv_afraid == 1 ~ "Never", ipv_afraid == 2 ~ "Once", ipv_afraid == 3 ~ "Several times", ipv_afraid == 4 ~"Most/all of the times", TRUE ~ "Missing"),
         ipv_life_freq=case_when(ipv_life_freq == 1 ~ "Never", ipv_life_freq == 2 ~ "Once", ipv_life_freq == 3 ~ "Several times", ipv_life_freq == 4 ~ "Most/all of the times", TRUE ~ "Missing"),
         dataid=as.numeric(dataid))

head(ipv_sv)

table(ipv_sv$ipv_life_freq)
dfull <- left_join(dfull, ipv_sv, by="dataid")
dim(dfull)
table(dfull$ipv_life_freq)

saveRDS(dfull, paste0(dropboxDir,"Data/Cleaned/Audrie/pregnancy_telo_covariates_data.RDS"))


