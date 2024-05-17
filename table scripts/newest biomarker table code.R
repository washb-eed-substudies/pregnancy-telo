
library(data.table)
source(here::here("0-config.R"))
library(tidyverse)
library(flextable)
library(officer)

d <- readRDS("/Users/farheenjamshed/Downloads/bangladesh-cleaned-master-data (2).RDS")

d <- d %>% filter(pregnancy_telo==1)

d <- d %>% mutate(vit_A_def = ifelse(RBP_inf_preg < 0.7, 1, 0),
                  vit_A_low = ifelse(RBP_inf_preg < 1.05, 1, 0))

d <- d %>% mutate(crp_high = ifelse(crp > 5, 1, 0),
                  agp_high = ifelse(agp > 1, 1, 0))


writeqntle<-function(vector) {
  quantiles<-round(quantile(vector, na.rm=TRUE), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
}

nperc <- function(vector){
  n <- sum(vector==1, na.rm=T)
  perc <- round(n/sum(!is.na(vector))*100)
  paste(n, " (", perc, "%)", sep="")
}



mom_lab <-c("Maternal Biomarker", 
            "RBP (umol/L)", "Low Vitamin A", "Vitamin A Deficiency", "25(OH)D (nmol/L)", "Vitamin D Deficiency", "Ferritin (ug/L)", "sTfR (mg/L)", "Iron Deficiency",
            "Cortisol (ug/dL)", "Estriol (ng/mL)",
            "IL-1 (pg/ml)", "Il-6 (pg/ml)", "TNF-a (pg/ml)", "IL-12 (pg/ml)", "IFN-g (pg/ml)", "IL-4 (pg/ml)", "IL-5 (pg/ml)", "IL-13 (pg/ml)", "IL-17A (pg/ml)", "IL-21 (pg/ml)", "IL-10 (pg/ml)", "IL-2 (pg/ml)", "GM-CSF (pg/ml)", "AGP (g/L)", "Elevated AGP", "CRP (mg/L)", "Elevated CRP")

view(mom_lab)

# ^WHAT IS ESTRI AND CORT CALLED IN MASTER DATASET
# FARHEEN: ALSO ADD "vit_A_def" "vit_D_def" "iron_def" <- THESE ARE VARIABLE NAMES IN MASTER DATASET
# SHOULD I BE ADDING ALL CYTOKINES LIKE THIS? SHOULD I ADD SOME RATIOS? IS THERE AN ABS SUM SCORE VARIABLE? I ONLY SEE Z SCORE VARIABLE
#        no, dont include ratios; doesnt make sense in context of interpretation; also dont add sum score
# use "il6_mom_t0" not "mean_il6_mom_t0"
#MEAN OR MEDIAN? SHOULD I DO SHAPIRO TEST FOR EVERYTHING? 
# ^^ dont do shapiro test; just use median and iqr

child_lab <-c("T/S Ratio at Age 14 months", "T/S Ratio at Age 28 months", "Change in T/S Ratio between 14 months and 28 months") 

view(child_lab)
# SHOULD I PUT TS RATIO OR BASE PAIRS? "delta_ts_bp" OR "delta_TS"
# ^ use T/S ratio / relative->absolute makes assumptions to get bp


mom <- c("Median (25th, 75th percentile) or n (%)", 
        writeqntle(d$rbp), nperc(d$vit_A_low), nperc(d$vit_A_def), writeqntle(d$vitD_nmol_per_L), nperc(d$vit_D_def), writeqntle(d$ferr), writeqntle(d$stfr), nperc(d$iron_def),
        writeqntle(d$preg_cort), writeqntle(d$preg_estri), 
        writeqntle(d$il1_mom_t0), writeqntle(d$il6_mom_t0), writeqntle(d$tnfa_mom_t0), writeqntle(d$il12_mom_t0), writeqntle(d$ifng_mom_t0), writeqntle(d$il4_mom_t0), writeqntle(d$il5_mom_t0), writeqntle(d$il13_mom_t0), writeqntle(d$il17_mom_t0), writeqntle(d$il21_mom_t0), writeqntle(d$il10_mom_t0), writeqntle(d$il2_mom_t0), writeqntle(d$gmcsf_mom_t0), writeqntle(d$agp), nperc(d$agp_high), writeqntle(d$crp), nperc(d$crp_high))

view(mom)
# make a proportion of moms with deficiencies

child <- c(writeqntle(d$TS_t2),
           writeqntle(d$TS_t3),
           writeqntle(d$delta_TS))

view(child)


mom_tbl<-data.table(" "= mom_lab,
                    "At Enrollment"= mom)

child_tbl<-data.table("Child Outcome"= child_lab,
                      "Median (25th, 75th percentile)" = child)

view(mom_tbl)
view(child_tbl)

sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width=8.5, height=11),
  page_margins = page_mar(bottom=.3, top=.3, right=.3, left=.3, gutter = 0)
)

save_as_docx("Maternal Biomarker Table" = flextable(mom_tbl), path=here("preg-telo maternal biomarkers table 051624.docx"), 
             pr_section = sect_properties) 

save_as_docx("Child Biomarker Table" = flextable(child_tbl), path=here("preg-telo child biomarkers table 051624.docx"), 
             pr_section = sect_properties) 


#write.csv(tbl, file=here('tables/biomarkers.csv'))
#print(xtable(tbl), type="html", file=here("tables/biomarkers.html"))