data <-  readRDS("/Users/farheenjamshed/Downloads/pregnancy_telo_covariates_data.RDS")

# enrollment characteristics
shapiro.test(data$momage)
shapiro.test(data$momheight)
shapiro.test(data$momeduy)
shapiro.test(data$cesd_sum_t2)
shapiro.test(data$cesd_sum_ee_t3)
shapiro.test(data$pss_sum_mom_t3)
shapiro.test(data$life_viol_any_t3)



# maternal nutrition
shapiro.test(data$iron_def)
shapiro.test(data$vit_D_def)
shapiro.test(data$vit_A_def)
shapiro.test(data$vitD_nmol_per_L)
shapiro.test(data$stfr)
shapiro.test(data$rbp)
shapiro.test(data$ferr)

# maternal hormones
shapiro.test(data$preg_estri)
shapiro.test(data$preg_cort)

# maternal inflammation
shapiro.test(data$agp)
shapiro.test(data$crp)
shapiro.test(data$il1_mom_t0)
shapiro.test(data$il6_mom_t0)
shapiro.test(data$tnfa_mom_t0)
shapiro.test(data$il2_mom_t0)
shapiro.test(data$il12_mom_t0)
shapiro.test(data$ifng_mom_t0)
shapiro.test(data$il4_mom_t0)
shapiro.test(data$il5_mom_t0)
shapiro.test(data$il13_mom_t0)
shapiro.test(data$il17_mom_t0)
shapiro.test(data$il21_mom_t0)
shapiro.test(data$il10_mom_t0)
shapiro.test(data$gmcsf_mom_t0)

# child inflammation
# 14 months
shapiro.test(data$agp_t2)
shapiro.test(data$crp_t2)
shapiro.test(data$il1_t2)
shapiro.test(data$il6_t2)
shapiro.test(data$tnfa_t2)
shapiro.test(data$il2_t2)
shapiro.test(data$il12_t2)
shapiro.test(data$ifng_t2)
shapiro.test(data$il4_t2)
shapiro.test(data$il5_t2)
shapiro.test(data$il13_t2)
shapiro.test(data$il17_t2)
shapiro.test(data$il21_t2)
shapiro.test(data$il10_t2)
shapiro.test(data$gmcsf_t2)

# 28 months
shapiro.test(data$agp_t3)
shapiro.test(data$crp_t3)
shapiro.test(data$il1_t3)
shapiro.test(data$il6_t3)
shapiro.test(data$tnfa_t3)
shapiro.test(data$il2_t3)
shapiro.test(data$il12_t3)
shapiro.test(data$ifng_t3)
shapiro.test(data$il4_t3)
shapiro.test(data$il5_t3)
shapiro.test(data$il13_t3)
shapiro.test(data$il17_t3)
shapiro.test(data$il21_t3)
shapiro.test(data$il10_t3)
shapiro.test(data$gmcsf_t3)