



rm(list=ls())

#------------------------------------------------------------------------------
# update washbglm function to drop near zero variance covariates to avoid error
#------------------------------------------------------------------------------


washb_glm2 <- function (Y, tr, pair = NULL, W = NULL, forcedW = NULL, V = NULL, 
                        id, contrast, family = "gaussian", pval = 0.2, print = TRUE, 
                        verbose = FALSE, FECR = NULL){
  #browser()
  require(sandwich)
  require(lmtest)
  options(scipen = 20)
  Subgroups = NULL
  if (!is.null(FECR)) {
    if (FECR != "arithmetic" & FECR != "geometric") {
      stop(paste("You specified FECR=", fnargs$FECR[[length(fnargs$FECR)]], 
                 "to estimate the fecal egg count reduction %\nYou need to supply either 'arithmetic' or 'geometric' as an argument to the FECR option."))
    }
    if ((FECR == "arithmetic" | FECR == "geometric") & family != 
        "gaussian") {
      stop(paste("You specified FECR=", fnargs$FECR[[length(fnargs$FECR)]], 
                 "to estimate the fecal egg count reduction %\nThis parameter is a ratio of means: FECR=(EY1/EY0)-1\nso you need to specify family='gaussian' to estimate it properly."))
    }
  }
  if (!is.null(W)) {
    W <- data.frame(W)
    if (sum("tr" %in% colnames(W)) > 0) {
      colnames(W)[which(colnames(W) == "tr")] <- "trW"
    }
  }
  if (!is.null(pair)) {
    if (!is.null(W)) {
      glmdat <- data.frame(id, Y, tr, pair, W)
    }else {
      glmdat <- data.frame(id, Y, tr, pair)
    }
    glmdat$tr <- factor(glmdat$tr, levels = contrast[1:2])
    glmdat$pair <- factor(glmdat$pair)
  }else {
    if (!is.null(W)) {
      glmdat <- data.frame(id, Y, tr, W)
    }else {
      glmdat <- data.frame(id, Y, tr)
    }
    glmdat$tr <- factor(glmdat$tr, levels = contrast[1:2])
  }
  glmdat <- subset(glmdat, tr == contrast[1] | tr == contrast[2])
  glmdat$tr <- factor(glmdat$tr, levels = contrast[1:2])
  if (!is.null(pair)) {
    n.orig <- dim(glmdat)[1]
    miss <- NULL
    activeOnly <- ((subset(glmdat, tr == contrast[1])))
    nomiss <- sort(unique(activeOnly$pair))
    miss1 <- (unique(pair)[which(!(unique(pair) %in% (nomiss)))])
    activeOnly2 <- ((subset(glmdat, tr == contrast[2])))
    nomiss2 <- sort(unique(activeOnly2$pair))
    miss2 <- (unique(pair)[which(!(unique(pair) %in% (nomiss2)))])
    miss <- append(miss1, miss2)
    glmdat <- subset(glmdat, !(pair %in% miss))
    n.sub <- dim(glmdat)[1]
    if (print == TRUE) 
      if (n.orig > n.sub) 
        cat("\n-----------------------------------------\n", 
            "Starting N:  ", n.orig, "\nN after block dropping: ", 
            n.sub)
    if (print == TRUE) 
      if (n.orig > n.sub) 
        cat("\n-----------------------------------------\n", 
            "Pairs/blocks dropped due to missingness in at least one treatment level:\n", 
            sort(unique(miss)), "\n\nDropping", n.orig - 
              n.sub, "observations due to missing pairs.", 
            "\n-----------------------------------------\n")
  }
  n.orig <- dim(glmdat)[1]
  rowdropped <- rep(1, nrow(glmdat))
  rowdropped[which(complete.cases(glmdat))] <- 0
  glmdat <- glmdat[complete.cases(glmdat), ]
  n.sub <- dim(glmdat)[1]
  if (print == TRUE) 
    if (n.orig > n.sub) 
      cat("\n-----------------------------------------\nDropping", 
          n.orig - n.sub, "observations due to missing values in 1 or more variables\n", 
          "Final sample size:", n.sub, "\n-----------------------------------------\n")
  if (!is.null(W)) {
    colnamesW <- names(W)
  }
  if (!is.null(W)) {
    if (!is.null(V)) {
      forcedW = c(V, forcedW)
    }
    if (!is.null(forcedW)) {
      screenW <- subset(glmdat, select = colnamesW)
      toexclude <- names(screenW) %in% forcedW
      if (length(which(toexclude == TRUE)) != length(forcedW)) 
        stop("A forcedW variable name is not a variable within the W data frame.")
      screenW = screenW[!toexclude]
      if (ncol(screenW) == 0) {
        screenW <- NULL
      }
      if (print == TRUE) {
        cat("\n-----------------------------------------\nInclude the following adjustment covariates without screening:\n-----------------------------------------\n")
        print(forcedW, sep = "\n")
      }
    }else {
      screenW <- subset(glmdat, select = colnamesW)
    }
  }else {
    screenW <- NULL
  }
  if (!is.null(screenW)) {
    if (print == TRUE) 
      cat("\n-----------------------------------------\nPre-screening the adjustment covariates:\n-----------------------------------------\n")
    if(length(nearZeroVar(screenW))>0){
      screenW<- screenW[,-nearZeroVar(screenW)]
    }
    suppressWarnings(Wscreen <- washb_prescreen(Y = glmdat$Y, 
                                                Ws = screenW, family = family, pval = pval, print = print))
  }
  else {
    Wscreen = NULL
  }
  if (!is.null(pair)) {
    if (!is.null(forcedW)) {
      if (!is.null(Wscreen)) {
        dmat <- subset(glmdat, select = c("Y", "tr", 
                                          forcedW, Wscreen, "pair"))
      }
      else {
        dmat <- subset(glmdat, select = c("Y", "tr", 
                                          forcedW, "pair"))
      }
    }
    else {
      if (!is.null(Wscreen)) {
        dmat <- subset(glmdat, select = c("Y", "tr", 
                                          Wscreen, "pair"))
      }
      else {
        dmat <- subset(glmdat, select = c("Y", "tr", 
                                          "pair"))
      }
    }
  }
  else {
    if (!is.null(forcedW)) {
      if (!is.null(Wscreen)) {
        dmat <- subset(glmdat, select = c("Y", "tr", 
                                          forcedW, Wscreen))
      }
      else {
        dmat <- subset(glmdat, select = c("Y", "tr", 
                                          forcedW))
      }
    }
    else {
      if (!is.null(Wscreen)) {
        dmat <- subset(glmdat, select = c("Y", "tr", 
                                          Wscreen))
      }
      else {
        dmat <- subset(glmdat, select = c("Y", "tr"))
      }
    }
  }
  if (family[1] == "binomial" | family[1] == "poisson" | family[1] == 
      "gaussian") {
    if (!is.null(FECR)) {
      if (print == TRUE) {
        cat(paste("\n-----------------------------------------\nEstimating the fecal egg count reduction\n(FECR) proportion = (EY1/EY0) - 1\nfrom GLM results using", 
                  FECR, "means\nand the delta method (for a ratio of means)\n-----------------------------------------\n"))
      }
      require(msm)
      if (!is.null(V)) {
        colnames(dmat)[which(colnames(dmat) == V)] <- "V"
        if (class(dmat$V) == "factor") 
          Subgroups <- levels(dmat$tr:dmat$V)
        if (class(dmat$V) != "factor") 
          warning("V is not a factor variable within the W covariate data frame. An interaction term will be added to the model but not linear combination of coefficients will be calculated.")
        suppressWarnings(fit <- glm(Y ~ tr * V + ., family = family, 
                                    data = dmat))
      }
      else {
        suppressWarnings(fit <- glm(Y ~ ., family = family, 
                                    data = dmat))
      }
      vcovCL <- sandwichSE(dmat, fm = fit, cluster = glmdat$id)
      rfit <- coeftest(fit, vcovCL)
      df1 <- df0 <- dmat
      df1$tr <- contrast[2]
      df0$tr <- contrast[1]
      Qst1 <- predict(fit, type = "response", newdata = df1)
      Qst0 <- predict(fit, type = "response", newdata = df0)
      Ey1 <- mean(Qst1, na.rm = T)
      Ey0 <- mean(Qst0, na.rm = T)
      modelfit <- washb_glmFormat(glmModel = fit, rfit = rfit, 
                                  dmat = dmat, rowdropped = rowdropped, contrast = contrast, 
                                  pair = pair, vcovCL = vcovCL, family = family, 
                                  V = V, Subgroups = Subgroups, print = print, 
                                  verbose = verbose)
      if (FECR == "arithmetic") {
        n_coef <- length(fit$coefficients)
        fecr <- (Ey1/Ey0) - 1
        vars <- paste0("x", 1:n_coef)
        numerator <- paste(vars, collapse = "+")
        denominator <- gsub("\\+x2\\+", "+", numerator)
        denominator <- gsub("x1\\+x2", "x1", denominator)
        delta_formula <- as.formula(paste0("~(", numerator, 
                                           ")/(", denominator, ")-1"))
        fecr_se <- deltamethod(g = delta_formula, mean = coef(fit), 
                               cov = vcovCL(fit, glmdat$id), ses = TRUE)
      }
      if (FECR == "geometric") {
        fecr <- (exp(Ey1)/exp(Ey0)) - 1
        fecr_se <- deltamethod(g = ~exp(x2) - 1, mean = coef(fit), 
                               cov = vcovCL(fit, glmdat$id), ses = TRUE)
      }
      fecr_lb <- fecr - 1.96 * fecr_se
      fecr_ub <- fecr + 1.96 * fecr_se
      fecr_p <- 2 * (1 - pnorm(abs(fecr/fecr_se)))
      if (print == TRUE) {
        cat(paste("\n-----------------------------------------\nFecal egg count reduction (EY1/EY0)-1,\nestimated using ", 
                  FECR, " means", "\n-----------------------------------------\n", 
                  sep = ""))
        cat(paste("FECR (95% CI) : ", sprintf("%1.3f", 
                                              fecr), " (", sprintf("%1.3f", fecr_lb), ", ", 
                  sprintf("%1.3f", fecr_ub), ")", sep = ""))
        cat("\n     SE(FECR) :", sprintf("%1.4f", fecr_se))
        cat("\n      p-value :", sprintf("%1.4f", fecr_p))
        cat("\n-----------------------------------------\n")
      }
      modelfit$TR <- data.frame(psi = fecr, var.psi = fecr_se^2, 
                                ci.lb = fecr_lb, ci.ub = fecr_ub, pvalue = fecr_p, 
                                method = FECR)
      return(modelfit)
    }
    else {
      if (!is.null(V)) {
        colnames(dmat)[which(colnames(dmat) == V)] <- "V"
        if (class(dmat$V) == "factor") 
          Subgroups <- levels(dmat$tr:dmat$V)
        if (class(dmat$V) != "factor") 
          warning("V is not a factor variable within the W covariate data frame. An interaction term will be added to the model but not linear combination of coefficients will be calculated.")
        suppressWarnings(fit <- glm(Y ~ tr * V + ., family = family, 
                                    data = dmat))
        vcovCL <- sandwichSE(dmat, fm = fit, cluster = glmdat$id)
        rfit <- coeftest(fit, vcovCL)
      }
      else {
        suppressWarnings(fit <- glm(Y ~ ., family = family, 
                                    data = dmat))
        vcovCL <- sandwichSE(dmat, fm = fit, cluster = glmdat$id)
        rfit <- coeftest(fit, vcovCL)
      }
      modelfit <- washb_glmFormat(glmModel = fit, rfit = rfit, 
                                  dmat = dmat, rowdropped = rowdropped, contrast = contrast, 
                                  pair = pair, vcovCL = vcovCL, family = family, 
                                  V = V, Subgroups = Subgroups, print = print, 
                                  verbose = verbose)
      return(modelfit)
    }
  }
  else {
    if (family[1] == "neg.binom") {
      require(MASS)
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("MASS needed for this function to work. Please install it.", 
             call. = FALSE)
      }
      if (!is.null(V)) {
        colnames(dmat)[which(colnames(dmat) == V)] <- "V"
        Subgroups <- levels(dmat$tr:dmat$V)
        if (class(dmat$V) != "factor") 
          warning("V is not a factor variable within the W covariate data frame. An interaction term will be added to the model but not linear combination of coefficients will be calculated.")
        suppressWarnings(fit <- glm.nb(Y ~ tr * V + ., 
                                       data = dmat))
        vcovCL <- sandwichSE(dmat, fm = fit, cluster = glmdat$id)
        rfit <- coeftest(fit, vcovCL)
      }
      else {
        suppressWarnings(fit <- glm.nb(Y ~ ., data = dmat))
        vcovCL <- sandwichSE(dmat, fm = fit, cluster = glmdat$id)
        rfit <- coeftest(fit, vcovCL)
      }
      modelfit <- washb_glmFormat(glmModel = fit, rfit = rfit, 
                                  dmat = dmat, rowdropped = rowdropped, contrast = contrast, 
                                  pair = pair, vcovCL = vcovCL, family = family, 
                                  V = V, Subgroups = Subgroups, print = print, 
                                  verbose = verbose)
      if (print == TRUE) 
        cat("\n-----------------------------------------\nAssess whether conditional mean is equal to conditional variance:\n-----------------------------------------\n")
      if (!is.null(V)) {
        pois <- glm(Y ~ tr * V + ., family = "poisson", 
                    data = dmat)
      }
      else {
        pois <- glm(Y ~ ., family = "poisson", data = dmat)
      }
      X2 <- 2 * (logLik(fit) - logLik(pois))
      Pois_LRtest <- pchisq(X2, df = 1, lower.tail = FALSE)
      if (print == TRUE) {
        cat("\nLog-likelihood ratio test P-value:\n")
        cat("\nIf <0.05, negative binomial model is more appropriate than a Poisson model.\n\n")
        print(Pois_LRtest)
      }
      modelfit <- c(modelfit, Pois_LRtest)
      return(modelfit)
    }
    else {
      stop("Error in family specified. Must choose Gaussian, Poisson, Binomial, Binomial(link-log), or neg.binom.")
    }
  }
}

#------------------------------------------------------------------------------
# load data and set up analysis
#------------------------------------------------------------------------------


source(here::here("0-config.R"))
library(caret)

dfull <- readRDS("C:/Users/andre/Documents//EE/eed-substudy-data/bangladesh-cleaned-master-data.RDS")
dfull <- dfull %>% filter(pregnancy_telo==1)
summary(dfull$gest_age_weeks)

dfull <- dfull %>% select(all_of(c("childid", "gest_age_weeks", "HHwealth_scaled")))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/pregnancy_telo_covariates_data.RDS"))
dim(d)
d <- d %>% left_join(dfull, by="childid")
dim(d)

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

Xvars <- c("ipv_child", "ipv_afraid", "ipv_life_freq")          
Yvars <- c("TS_t3_Z", "delta_TS_Z")




#------------------------------------------------------------------------------
# run analysis
#------------------------------------------------------------------------------


#Fit models
freq_res <- NULL
for(i in Xvars){
  for(j in Yvars){
    cat(i,"\n")
    cat(j,"\n")
    
    if(grepl("_t2",j)){
      df <- d %>% select(all_of(c("childid",i,j,Wvars2)))
      df <- df[complete.cases(df),]
      Ws = df %>% select(all_of(Wvars2))
      if(length(nearZeroVar(Ws))>0){
        Ws<- Ws[,-nearZeroVar(Ws)]
      }
    }
    if(grepl("_t3",j)){
      df <- d %>% select(all_of(c("childid",i,j,Wvars3)))
      df <- df[complete.cases(df),]
      Ws = df %>% select(all_of(Wvars3))
      if(length(nearZeroVar(Ws))>0){
        Ws<- Ws[,-nearZeroVar(Ws)]
      }
    }
    if(grepl("delta",j)){
      df <- d %>% select(all_of(c("childid",i,j,Wvars23)))
      df <- df[complete.cases(df),]
      Ws = df %>% select(all_of(Wvars23))
      if(length(nearZeroVar(Ws))>0){
        Ws<- Ws[,-nearZeroVar(Ws)]
      }
    }

    res=NULL
    for(k in c("Once","Several times","Most/all of the times")){
      temp<-NULL
      if(sum(df[[i]]==k,na.rm = T)>0){
        temp <- washb_glm2(Y=df[[j]], tr=df[[i]], 
                          W = Ws, forcedW = NULL, V = NULL, 
                          id=df$childid, contrast=c("Never",k), family = "gaussian", print = F)
      }
      if(!is.null(temp)){
        temp <- temp$TR
        temp$exposure <- j
        temp$outcome <- i
        temp$contrast <- k
      }
      res=bind_rows(res, temp)
    }
    
    freq_res <- bind_rows(freq_res, res)
  }
  rownames(freq_res) <- NULL
}


#format for plot data
freq_res <- freq_res %>% mutate(contrast=factor(contrast, levels=c("Never","Once","Several times","Most/all of the times")),
                                outcome_f=case_when(outcome=="ipv_child" ~ "Child present/overhear physical violence",
                                                    outcome=="ipv_afraid" ~ "Fear of husband",
                                                    outcome=="ipv_life_freq" ~ "Injury due to acts by husband"),
                                exposure_f=case_when(exposure=="TS_t3_Z" ~ "Telomeres Year 2",
                                                     exposure=="delta_TS_Z" ~ "Telomere Change Year 1 and Year 2"))
labels[labels$name=="q_1002",]
labels[labels$name=="q_808",]
labels[labels$name=="q_902_a",]

#make plot
p <- ggplot(aes(x=contrast, y=`Coef.`, ymin=`2.5%`, ymax=`97.5%`, color= outcome_f), data=freq_res) +
  geom_point() +
  geom_errorbar(width=0.1) +
  geom_hline(yintercept = 0) +
  facet_wrap(~exposure_f + outcome_f, scales="free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  theme(legend.position = "none") +
  labs(x="Contrast (ref: never)", y="Adjusted mean difference") 

p
ggsave(p, file = here("figures/ipv-freq-telo-figure.jpg"), height=8, width=8)

