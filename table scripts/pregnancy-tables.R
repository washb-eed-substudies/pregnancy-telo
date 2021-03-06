rm(list=ls())

library('flextable')
library('officer')
# source(here::here("0-config.R"))

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

#### Functions for growth tables ####

pregnancy_tbl <- function(name, expo_var, out_var, exposure, outcome, H1, H1_adj){
  

  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  
  ### this function produces a table that can be saved as a csv
  
  tbl <- data.table(name = character(), "Outcome" = character(), "N" = character(), "25th Percentile" = character(), "75th Percentile" = character(),
                    " Outcome, 75th Percentile v. 25th Percentile" = character(), " " = character(), " " = character(), " " = character(), " " = character(),
                    " " = character(), " " = character(), " " = character(), " " = character(), " " = character())
  tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", "Unadjusted", " ", " ", " ", " ", "Fully adjusted", " ", " ", " ", " "))
  tbl <- rbind(tbl, list(" ", " ", " ", " ", " ", 
                         "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value", "FDR adjusted P-value", 
                         "Predicted Outcome at 25th Percentile", "Predicted Outcome at 75th Percentile", "Coefficient (95% CI)", "P-value", "FDR adjusted P-value"))
  skipped<-F
  for (i in 1:length(exposure)) {
    for (j in 1:length(outcome)) {
      exp <- exposure[i]
      out <- outcome[j]
      filtered_res <- results[results$Y==out & results$X==exp,]
      filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
      unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
      adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
      if (nrow(filtered_res)==0){
        skipped<-T
        next
      }
      if(j==1|skipped==T){
        tbl <- rbind(tbl, list(expo_var[i], out_var[j], filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
        skipped<-F
      }else {
        tbl <- rbind(tbl, list("", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
      }
    }
    if (i != length(exposure)) {
      tbl <- rbind(tbl, list("","","","","","","","","","","","","","",""))
    }
  }
  tbl
}

pregnancy_tbl_flex <- function(name, expo_var, out_var, exposure, outcome, results, results_adj){
  ### name: string name of group of exposures
  ### expo_var: vector of string exposures to include in table
  ### out_var: vector of string outcomes to include in table
  ### exposure: vector of string exposure variable names
  ### outcome: vector of string outcome variable names
  ### results: data frame with unadjusted results
  ### results_adj: data fram with adjusted results
  
  ### this function produces a table that can be saved as an image or 
  ### directly to a word document!
  
  # build table
  tbl <- data.table(matrix(nrow=0, ncol=15))
  skipped<-F
  for (i in 1:length(exposure)) {
    for (j in 1:length(outcome)) {
      exp <- exposure[i]
      out <- outcome[j]
      filtered_res <- results[results$Y==out & results$X==exp,]
      filtered_adj <- results_adj[results_adj$Y==out & results_adj$X==exp,]
      unadj <- paste(round(filtered_res$`point.diff`, 2), " (", round(filtered_res$`lb.diff`, 2), ", ", round(filtered_res$`ub.diff`, 2), ")", sep="")
      adj <- paste(round(filtered_adj$`point.diff`, 2), " (", round(filtered_adj$`lb.diff`, 2), ", ", round(filtered_adj$`ub.diff`, 2), ")", sep="")
      if (nrow(filtered_res)==0){
        skipped<-T
        next
      }
      if(j==1|skipped==T){
        tbl <- rbind(tbl, list(expo_var[i], out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
        skipped=F
      }else {
        tbl <- rbind(tbl, list(" ", out_var[j],  filtered_res$N, round(filtered_res$q1, 2), round(filtered_res$q3, 2), 
                               round(filtered_res$pred.q1, 2), round(filtered_res$pred.q3, 2), unadj, round(filtered_res$Pval, 2), round(filtered_res$BH.Pval, 2), 
                               round(filtered_adj$pred.q1, 2), round(filtered_adj$pred.q3, 2), adj, round(filtered_adj$Pval, 2), round(filtered_adj$BH.Pval, 2)))
      }
    }
    if (i != length(exposure)) {
      tbl <- rbind(tbl, list("","","","","","","","", "","","","","","",""))
    }
  }
  
  # format for export
  flextbl<-flextable(tbl, col_keys=names(tbl))
  flextbl <- set_header_labels(flextbl,
                               values = list("V1" = " ", "V2" = " ", "V3" = " ", "V4" = " ", "V5" = " ",
                                             "V6" = "Predicted Outcome at 25th Percentile", "V7" = "Predicted Outcome at 75th Percentile", "V8" = "Coefficient (95% CI)", "V9" = "P-value", "V10" = "FDR Corrected P-value",
                                             "V11" = "Predicted Outcome at 25th Percentile", "V12" = "Predicted Outcome at 75th Percentile", "V13" = "Coefficient (95% CI)", "V14" = "P-value", "V15" = "FDR Corrected P-value"))
  flextbl <- add_header_row(flextbl, values = c("","","","","", "Unadjusted", "Fully adjusted"), colwidths=c(1,1,1,1,1,5,5))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- add_header_row(flextbl, values = c(name, "Outcome","N","25th Percentile","75th Percentile", "Outcome, 75th Percentile v. 25th Percentile"), colwidths=c(1,1,1,1,1,10))
  # flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- hline_bottom(flextbl, part="body", border=fp_border(color="black"))
  flextbl <- hline_top(flextbl, part="header", border=fp_border(color="black"))
  flextbl <- align(flextbl, align = "center", part = "all")
  flextbl <- align(flextbl, j = c(1, 2), align = "left", part="all")
  flextbl <- autofit(flextbl, part = "all")
  flextbl <- fit_to_width(flextbl, max_width=8)
  
  flextbl
}


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

exposure <- c("vitD_nmol_per_L", "logFERR_inf", "logSTFR_inf", "logRBP_inf", "vit_A_def", "iron_def", "vit_D_def")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Vitamin D", "Ferritin", "sTfR", "RBP", "Vit A Deficiency", "Iron Deficiency", "Vitamin D Deficiency")
out_var <- c("Telomeres Year 1", "Telomeres Year 2", "Telomere Change Year 1 and Year 2")
results <- H1
results_adj <- H1adj

tbl2 <- pregnancy_tbl("Nutrition Biomarkers", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl2flex <- pregnancy_tbl_flex("Nutrition Biomarkers", expo_var, out_var, exposure, outcome, H1, H1adj)

#### Table 3 ####

exposure <- c("preg_cort")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Cortisol")
out_var <- c("Telomeres Year 1", "Telomeres Year 2", "Telomere Change Year 1 and Year 2")
results <- H2
results_adj <- H2adj

tbl3 <- pregnancy_tbl("Serum Stress Biomarker", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl3flex <- pregnancy_tbl_flex("Serum stress biomarker", expo_var, out_var, exposure, outcome, H2, H2adj)

#### Table 4 ####

exposure <- c("logCRP", "logAGP", "ifng_mom_t0", "sumscore_t0_mom_Z")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("CRP", "AGP", "IFNG", "Sum Score")
out_var <- c("Telomeres Year 1", "Telomeres Year 2", "Telomere Change Year 1 and Year 2")
results <- H3
results_adj <- H3adj

tbl4 <- pregnancy_tbl("Inflammation Biomarkers", expo_var, out_var, exposure, outcome, H3, H3adj)
tbl4flex <- pregnancy_tbl_flex("Inflammation Biomarkers", expo_var, out_var, exposure, outcome, H3, H3adj)


#### Table 5 ####

exposure <- c("preg_estri")
outcome <- c("TS_t2_Z", "TS_t3_Z", "delta_TS_Z")
expo_var <- c("Estriol")
out_var <- c("Telomeres Year 1", "Telomeres Year 2", "Telomere Change Year 1 and Year 2")
results <- H4
results_adj <- H4adj

tbl5 <- pregnancy_tbl("Estriol", expo_var, out_var, exposure, outcome, H4, H4adj)
tbl5flex <- pregnancy_tbl_flex("Estriol", expo_var, out_var, exposure, outcome, H4, H4adj)


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
write.csv(tbl3, here('tables/pregnancy-telo-table2.csv'))
write.csv(tbl4, here('tables/pregnancy-telo-table3.csv'))
write.csv(tbl5, here('tables/pregnancy-telo-table4.csv'))

save_as_docx("Table 2" = tbl2flex, "Table 3" = tbl3flex, "Table 4" = tbl4flex, "Table 5" = tbl5flex, path='~/Desktop/pregnancy-telo/tables/tables.docx')

## "preg-telo tables.docx" is easier to read than "tables.docx"

# write.csv(tbls1, here('tables/supplementary/stress-growth-tables1.csv'))
# write.csv(tbls2, here('tables/supplementary/stress-growth-tables2.csv'))
# write.csv(tbls3, here('tables/supplementary/stress-growth-tables3.csv'))
# write.csv(tbls4, here('tables/supplementary/stress-growth-tables4.csv'))

# save_as_docx("Table S1" = tbls1flex, "Table S2" = tbls2flex, "Table S3" = tbls3flex, "Table S4" = tbls4flex, path='C:/Users/Sophia/Documents/WASH/WASH Stress and Growth/stress-growth supplementary.docx')

