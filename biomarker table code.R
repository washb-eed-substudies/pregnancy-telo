
library(data.table)
source(here::here("0-config.R"))
library(tidyverse)
library(flextable)
library(officer)
library(boxr)
box_auth()

d <- readRDS("/Users/farheenjamshed/Downloads/bangladesh-cleaned-master-data (1).RDS")
colnames(d)

writeqntle<-function(vector) {
  quantiles<-round(quantile(vector, na.rm=TRUE), 2)
  paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
}

mom_lab <-c("Outcome", 
            "Vitamin D", "Ferritin", "sTfR", "RBP",
            "Vitamin A Deficiency", "Vitamin D Deficiency", "Iron Deficiency",
            "Cortisol", "Estriol",
            "IL-1 (pg/ml)", "Il-6 (pg/ml)", "TNF-a (pg/ml)", "IL-12 (pg/ml)", "IFN-g (pg/ml)", "IL-4 (pg/ml)", "IL-5 (pg/ml)", "IL-13 (pg/ml)", "IL-17A (pg/ml)", "IL-21 (pg/ml)", "IL-10 (pg/ml)", "IL-2 (pg/ml)", "GM-CSF (pg/ml)", "AGP (g/L)", "CRP (mg/L)")

view(mom_lab)

   # ^WHAT IS ESTRI AND CORT CALLED IN MASTER DATASET
   # FARHEEN: ALSO ADD "vit_A_def" "vit_D_def" "iron_def" <- THESE ARE VARIABLE NAMES IN MASTER DATASET
   # SHOULD I BE ADDING ALL CYTOKINES LIKE THIS? SHOULD I ADD SOME RATIOS? IS THERE AN ABS SUM SCORE VARIABLE? I ONLY SEE Z SCORE VARIABLE
   #        no, dont include ratios; doesnt make sense in context of interpretation; also dont add sum score
# use "il6_mom_t0" not "mean_il6_mom_t0"
#MEAN OR MEDIAN? SHOULD I DO SHAPIRO TEST FOR EVERYTHING? 
     # ^^ dont do shapiro test; just use median and iqr

child_lab <-c("Outcome", 
              "TS_t2", "TS_t3", "delta_TS") 

view(child_lab)
   # SHOULD I PUT TS RATIO OR BASE PAIRS? "delta_ts_bp" OR "delta_TS"
       # ^ use T/S ratio / relative->absolute makes assumptions to get bp

mom <-c("Median (25th, 75th percentile)", 
        writeqntle(d$vitD_nmol_per_L), writeqntle(d$ferr), writeqntle(d$stfr), writeqntle(d$rbp),
        writeqntle(d$vit_A_def), writeqntle(d$vit_D_def), writeqntle(d$iron_def),
        writeqntle(d$preg_cort), writeqntle(d$preg_estri), 
        writeqntle(d$il1_mom_t0), writeqntle(d$il6_mom_t0), writeqntle(d$tnfa_mom_t0), writeqntle(d$il12_mom_t0), writeqntle(d$ifng_mom_t0), writeqntle(d$il4_mom_t0), writeqntle(d$il5_mom_t0), writeqntle(d$il13_mom_t0), writeqntle(d$il17_mom_t0), writeqntle(d$il21_mom_t0), writeqntle(d$il10_mom_t0), writeqntle(d$il2_mom_t0), writeqntle(d$gmcsf_mom_t0), writeqntle(d$agp), writeqntle(d$crp))

# make a proportion of moms with deficiencies


child_t2 <- c("Median (25th, 75th percentile)", writeqntle(d$TS_t2))

child_t3<-c("Median (25th, 75th percentile)", writeqntle(d$TS_t3))

child_delta<-c("Median (25th, 75th percentile)", writeqntle(d$delta_TS))

view(child_t2)
view(child_t3)
view(child_delta)


mom_tbl<-data.table(" "= mom_lab,
                    "At Enrollment"=mom)

child_tbl<-data.table(" "= child_lab,
                      "Age 14 Months"=child_t2,
                      "Age 28 Months"=child_t3,
                      "Change in Year 1 and Year 2"=child_delta)

view(mom_tbl)
view(child_tbl)

sect_properties <- prop_section(
  page_size = page_size(orient = "portrait", width=8.5, height=11),
  page_margins = page_mar(bottom=.3, top=.3, right=.3, left=.3, gutter = 0)
)

save_as_docx("Maternal Biomarkers" = flextable(mom_tbl), path=here("preg-telo maternal biomarkers table 042922.docx"), 
             pr_section = sect_properties) 

save_as_docx("Child Biomarkers" = flextable(child_tbl), path=here("preg-telo child biomarkers table 042922.docx"), 
             pr_section = sect_properties) 


#write.csv(tbl, file=here('tables/biomarkers.csv'))
#print(xtable(tbl), type="html", file=here("tables/biomarkers.html"))