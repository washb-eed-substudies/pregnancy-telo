rm(list=ls())

library('flextable')
library('officer')
library('here')
source(here("table-functions.R"))
source(here::here("0-config.R"))

# load enrollment characteristics and results
# d <- read.csv(paste0(dropboxDir, "Data/Cleaned/Audrie/bangladesh-dm-ee-stress-growth-covariates-stresslab-anthro.csv"))
H1 <- readRDS(here('results/unadjusted/H1_res.RDS'))
H2 <- readRDS(here('results/unadjusted/H2_res.RDS'))
H3 <- readRDS(here('results/unadjusted/H3_res.RDS'))
H4 <- readRDS(here('results/unadjusted/H4_res.RDS'))
H1adj <- readRDS(here('results/adjusted/H1_adj_res.RDS'))
H2adj <- readRDS(here('results/adjusted/H2_adj_res.RDS'))
H3adj <- readRDS(here('results/adjusted/H3_adj_res.RDS'))
H4adj <- readRDS(here('results/adjusted/H4_adj_res.RDS'))


#### MAIN TABLES ####
#### Table 1 ####
# Characteristics of participants
# nperc <- function(vector){
#   n <- sum(vector==1, na.rm=T)
#   perc <- round(n/sum(!is.na(vector))*100)
#   paste(n, " (", perc, "%)", sep="")
# }
# 
# mediqr <- function(vector){
#   quantiles <- round(quantile(vector, na.rm=T), 2)
#   paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
# }
# 
# n_med_col <- c(nperc(d$sex), mediqr(d$t2_f2_8ip), mediqr(d$t2_f2_23d), mediqr(d$t2_f2_VI), mediqr(d$t2_f2_12i),
#                mediqr(d$t3_cort_slope), mediqr(d$t3_residual_cort), mediqr(d$t3_saa_slope), mediqr(d$t3_residual_saa),
#                mediqr(d$t3_map), mediqr(d$t3_hr_mean), mediqr(d$t3_gcr_mean), mediqr(d$t3_gcr_cpg12),
#                mediqr(d$laz_t2), mediqr(d$waz_t2), mediqr(d$whz_t2), mediqr(d$hcz_t2),
#                mediqr(d$laz_t3), mediqr(d$waz_t3), mediqr(d$whz_t3), mediqr(d$hcz_t3),
#                nperc(d$diar7d_t2), nperc(d$diar7d_t3), mediqr(d$momage), mediqr(d$momheight), 
#                mediqr(d$momeduy), mediqr(d$cesd_sum_t2), mediqr(d$cesd_sum_ee_t3), mediqr(d$pss_sum_mom_t3), 
#                nperc(d$life_viol_any_t3))
# 
# tbl1 <- data.table("C1" = c("Child","","","","","","","","","","","","","","","","","","","","","","","Mother","","","","","",""),
#                    "C2" = c("", "Urinary F2-isoprostanes (Year 1)","","","", "Salivary cortisol reactivity (Year 2)","", "sAA reactivity (Year 2)","",
#                            "SAM biomarkers (Year 2)","", "Glucocorticoid receptor","", "Anthropometry (14 months, Year 1)","","","",
#                            "Anthropometry (28 months, Year 2)","","","", "Diarrhea (14 months, Year 1)", "Diarrhea (28 months, Year 2)","",
#                            "Anthropometry at enrollment", "Education", "Depression at Year 1", "Depression at Year 2", "Perceived stress at Year 2", 
#                            "Intimate partner violence"),
#                    "C3" = c("Female", "iPF(2a)-III", "2,3-dinor-iPF(2a)-III", "iPF(2a-VI", "8,12-iso-iPF(2a)-VI", 
#                            "Change in slope between pre- and post-stressor cortisol", "Cortisol residualized gain score", 
#                            "Change in slope between pre- and post-stressor sAA change", "sAA residualized gain score",
#                            "Mean arterial pressure", "Resting heart rate", "NR3C1 exon 1F promoter methylation", "NGFI-A transcription factor binding site methylation",
#                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
#                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
#                            "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Age (years)", "Height (cm)", "Schooling completed (years)",
#                            "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
#                    "C4" = n_med_col)
# 
# tbl1flex <- flextable(tbl1, col_keys=names(tbl1))
# tbl1flex <- set_header_labels(tbl1flex,
#                         values = list("C1" = "", "C2" = "", "C3" = "", "C4" = "n (%) or median (IQR)"))
# tbl1flex <- hline_top(tbl1flex, part="header", border=fp_border(color="black", width = 1))
# tbl1flex <- hline_bottom(tbl1flex, part="all", border=fp_border(color="black", width = 1))
# tbl1flex <- autofit(tbl1flex, part = "all")
# tbl1flex <- align(tbl1flex, j = c(1, 2, 3), align = "left", part="all")
# tbl1flex <- align(tbl1flex, j = 4, align = "center", part="all")
# tbl1flex <- fit_to_width(tbl1flex, max_width=8)
# names(tbl1)<- c("","","","n (%) or median (IQR)")


#### Table 2 ####
exposure <- c("logRBP_inf", "vit_A_low", "vit_A_def", "vitD_nmol_per_L", "vit_D_def", "logFERR_inf", "logSTFR_inf", "iron_def")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Log RBP (umol/L)", "Low Vitamin A", "Vit A Deficiency", "25(OH)D (nmol/L)", "Vitamin D Deficiency", "Log ferritin (ug/L)", "Log sTfR (mg/L)", "Iron Deficiency")
out_var <- c("T/S Ratio Z-score Age 14 months", "T/S Ratio Z-score Age 28 months", "Change in T/S Ratio Z-score between 14 months and 28 months")
results <- H1
results_adj <- H1adj

tbl2 <- growth_tbl("Maternal micronutrients and child telomere length", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl2flex <- growth_tbl_flex("Maternal micronutrients and child telomere length", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl2supp <- growth_tbl("Maternal micronutrients and child telomere length", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl2flexsupp <- growth_tbl_flex("Maternal micronutrients and child telomere length", expo_var, out_var, exposure, outcome, H1, H1adj)

#### Table 3 ####

exposure <- c("ln_preg_cort")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Ln cortisol (ug/dL)")
out_var <- c("T/S Ratio Z-score Age 14 months", "T/S Ratio Z-score Age 28 months", "Change in T/S Ratio Z-score between 14 months and 28 months")
results <- H2
results_adj <- H2adj

tbl3 <- growth_tbl("Maternal cortisol and child telomere length", expo_var, out_var, exposure, outcome, H2, H2adj, T)
tbl3flex <- growth_tbl_flex("Maternal cortisol and child telomere length", expo_var, out_var, exposure, outcome, H2, H2adj, T)
tbl3supp <- growth_tbl("Maternal cortisol and child telomere length", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl3flexsupp <- growth_tbl_flex("Maternal cortisol and child telomere length", expo_var, out_var, exposure, outcome, H2, H2adj)


#### Table 4 ####

exposure <- c("logCRP", "logAGP", "mom_t0_ln_ifn", "sumscore_t0_mom_Z")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Ln CRP (mg/L)", "Ln AGP (g/L)", "Ln IFN-g (pg/ml)", "Inflammation Sum Score of 13 Cytokines")
out_var <- c("T/S Ratio Z-score Age 14 months", "T/S Ratio Z-score Age 28 months", "Change in T/S Ratio Z-score between 14 months and 28 months")
results <- H3
results_adj <- H3adj

tbl4 <- growth_tbl("Maternal inflammation and child telomere length", expo_var, out_var, exposure, outcome, H3, H3adj, T)
tbl4flex <- growth_tbl_flex("Maternal inflammation and child telomere length", expo_var, out_var, exposure, outcome, H3, H3adj, T)
tbl4supp <- growth_tbl("Maternal inflammation and child telomere length", expo_var, out_var, exposure, outcome, H3, H3adj)
tbl4flexsupp <- growth_tbl_flex("Maternal inflammation and child telomere length", expo_var, out_var, exposure, outcome, H3, H3adj)

#### Table 5 ####

exposure <- c("ln_preg_estri")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Ln estriol (ng/mL)")
out_var <- c("T/S Ratio Z-score Age 14 months", "T/S Ratio Z-score Age 28 months", "Change in T/S Ratio Z-score between 14 months and 28 months")
results <- H4
results_adj <- H4adj

tbl5 <- growth_tbl("Maternal estriol and child telomere length", expo_var, out_var, exposure, outcome, H4, H4adj, T)
tbl5flex <- growth_tbl_flex("Maternal estriol and child telomere lengt", expo_var, out_var, exposure, outcome, H4, H4adj, T)
tbl5supp <- growth_tbl("Maternal estriol and child telomere lengt", expo_var, out_var, exposure, outcome, H4, H4adj)
tbl5flexsupp <- growth_tbl_flex("Maternal estriol and child telomere lengt", expo_var, out_var, exposure, outcome, H4, H4adj)

#### Supplementary Tables ####
#### Table S1 ####

# exposure <- c("t2_f2_8ip", "t2_f2_23d", "t2_f2_VI", "t2_f2_12i", "iso.pca")
# outcome <- c("waz_t2", "whz_t2", "hcz_t2", "waz_t3", "whz_t3", "hcz_t3",
#  "wei_velocity_t2_t3", "hc_velocity_t2_t3",
#  "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3")
# expo_var <- c("IPF(2a)-III", "2,3-dinor-iPF(2a)-III", "iPF(2a)-VI", "8,12-iso-iPF(2a)-VI", "Combined urinary oxidative stress biommarkers")
# out_var <- c("WAZ Year 1", "WLZ Year 1", "HCZ Year 1",
#   "WAZ Year 2", "WLZ Year 2", "HCZ Year 2",
#  "Weight velocity (kg/month) Year 1 to Year 2",
# "Head circumference velocity (cm/month) Year 1 to Year 2",
#   "Change in child WAZ from Year 1 to Year 2",
#   "Change in WLZ from Year 1 to Year 2",
#   "Change in HCZ from Year 1 to Year 2")

# tbls1 <- growth_tbl("Urinary oxidative stress biomarker", expo_var, out_var, exposure, outcome, H1, H1adj)
# tbls1flex <- growth_tbl_flex("Urinary oxidative stress biomarker", expo_var, out_var, exposure, outcome, H1, H1adj)


#### Table S2 ####

# exposure <- c("t3_cort_slope", "t3_residual_cort", "t3_saa_slope", "t3_residual_saa")
# outcome <- c("waz_t3", "whz_t3", "hcz_t3")
# expo_var <- c("Pre to post-stress change in slope of cortisol", "Cortisol residualized gain score", "Pre to post-stress change in slope of sAA", "sAA residualized gain score")
# out_var <- c("WAZ Year 2", "WLZ Year 2", "HCZ Year 2")

# tbls2 <- growth_tbl("Salivary stress biomarker", expo_var, out_var, exposure, outcome, H2, H2adj)
# tbls2flex <- growth_tbl_flex("Salivary stress biomarker", expo_var, out_var, exposure, outcome, H2, H2adj)


#### Table S3 ####

# exposure <- c("t3_map", "t3_hr_mean")
# outcome <- c("waz_t3", "whz_t3", "hcz_t3")
# expo_var <- c("Mean arterial pressure", "Mean resting heart rate")
# out_var <- c("WAZ Year 2", "WLZ Year 2", "HCZ Year 2")

# tbls3 <- growth_tbl("Resting SAM biomarker", expo_var, out_var, exposure, outcome, H3, H3adj)
# tbls3flex <- growth_tbl_flex("Resting SAM biomarker", expo_var, out_var, exposure, outcome, H3, H3adj)



#### Table S4 ####

# exposure <- c("t3_gcr_mean", "t3_gcr_cpg12")
# outcome <- c("waz_t3", "whz_t3", "hcz_t3")
# expo_var <- c("Entire promoter region (39 assayed CpG sites)", "NGFI-A transcription factor binding site (CpG site #12)")
# out_var <- c("WAZ Year 2", "WLZ Year 2", "HCZ Year 2")

# tbls4 <- growth_tbl("Methylation site", expo_var, out_var, exposure, outcome, H4, H4adj)
# tbls4flex <- growth_tbl_flex("Methylation site", expo_var, out_var, exposure, outcome, H4, H4adj)


#### SAVE TABLES ####

# write.csv(tbl1, file=here("tables/main/pregnancy-telo-table1.csv"))
write.csv(tbl2, here('tables/pregnancy-telo-table1.csv'))
write.csv(tbl2supp, here('tables/pregnancy-telo-table1-supp.csv'))
write.csv(tbl3, here('tables/pregnancy-telo-table2.csv'))
write.csv(tbl3supp, here('tables/pregnancy-telo-table2-supp.csv'))
write.csv(tbl4, here('tables/pregnancy-telo-table3.csv'))
write.csv(tbl4supp, here('tables/pregnancy-telo-table3-supp.csv'))
write.csv(tbl5, here('tables/pregnancy-telo-table4.csv'))
write.csv(tbl5supp, here('tables/pregnancy-telo-table4-supp.csv'))

save_as_docx("Table 1" = tbl2flex, "Table 2" = tbl3flex, "Table 3" = tbl4flex, "Table 4" = tbl5flex,path='~/Desktop/pregnancy-telo/tables/maintables.docx')

save_as_docx("Table S1" = tbl2flexsupp, "Table S2" = tbl3flexsupp, "Table S3" = tbl4flexsupp, "Table S4" = tbl5flexsupp ,path='~/Desktop/pregnancy-telo/tables/supptables.docx')

# write.csv(tbls1, here('tables/supplementary/stress-growth-tables1.csv'))
# write.csv(tbls2, here('tables/supplementary/stress-growth-tables2.csv'))
# write.csv(tbls3, here('tables/supplementary/stress-growth-tables3.csv'))
# write.csv(tbls4, here('tables/supplementary/stress-growth-tables4.csv'))

# save_as_docx("Table S1" = tbls1flex, "Table S2" = tbls2flex, "Table S3" = tbls3flex, "Table S4" = tbls4flex, path='C:/Users/Sophia/Documents/WASH/WASH Stress and Growth/stress-growth supplementary.docx')

