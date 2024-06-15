
#-------------------------------------
# EE substudies analysis 

# configure data directories
# source base functions
# load libraries
#-------------------------------------

library(tidyverse)
library(haven)
library(washb)
library(foreign)
library(data.table)
library(tmle)
library(SuperLearner)
library(devtools)
library(kableExtra)
library(here)
library(cowplot)
library(mgcv)
library(psych)

if(!require(washbgam)){
  devtools::install_github("washb-eed-substudies/washbgam")
  library(washbgam)
}

dropboxDir <- NULL
if(dir.exists("C:/Users/andre/Dropbox/WASHB-EE-analysis/WBB-EE-analysis/")){ 
  dropboxDir <- "C:/Users/andre/Dropbox/WASHB-EE-analysis/WBB-EE-analysis/"
}
if(dir.exists("/Users/audrielin/Dropbox/WBB-EE-analysis/")){ 
  dropboxDir <- "/Users/audrielin/Dropbox/WBB-EE-analysis/"
}
if(dir.exists("C:/Users/Sophia/Dropbox/WASH/")){ 
  dropboxDir <- "C:/Users/Sophia/Dropbox/WASH/"
}
if(dir.exists("/Users/lisa/Dropbox/WASH/")){ 
  dropboxDir <- "/Users/lisa/Dropbox/WASH/"
}



theme_ki<-function(){
  theme_bw() %+replace%
    theme(
      strip.background = element_blank(),
      legend.position="none",
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size=14),
      axis.title = element_text(size=12),
      axis.text.y = element_text(size=10),
      axis.text.x = element_text(size=10, angle = 0, hjust = 0.5, vjust=.1)
    )
}

theme_set(theme_ki())

tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F",
               "#BCBD22","#17BECF")


#save R package versions

# # Only run thise lines once when project is initialized 
# #Call renv::init() to initialize a new project-local environment with a private R library,
# renv::init(project=here()) 
# 
# # Only run thise line when packages are updated
# #Call renv::snapshot() to save the state of the project library to the lockfile (called renv.lock),
# renv::snapshot()
# 
# # Only run these lines when needed (upon initialization and then when package versions need to be restored)
# #call renv::restore() to  revert to the previous state as encoded in the lockfile 
# renv::restore()


washb_tmle_ipcw = function(Y, tr, W = NULL, id = 1:length(Y), pair = NULL, Delta = rep(1, length(Y)), family = "gaussian", contrast, Q.SL.library = "SL.glm", 
                           g.SL.library = "SL.glm", pval = 0.2, FECR = NULL, seed = NULL, print = TRUE){
  
  require(tmle)
  require(SuperLearner)
  fnargs <- as.list(match.call())
  if (!is.null(W) & is.null(names(W))) {
    W <- data.frame(W)
  }
  Wvars = colnames(W)
  if (!is.null(pair)) {
    cat("\n-----------------------------------------", 
        "\nBy specifying the pair argument,\nyou have indicated that this is a matched pair analysis.\n\nNote: the analysis will only include pairs\nthat have a contrast in the treatment variable.", 
        "\n-----------------------------------------")
    if (length(fnargs$id) == 0) 
      fnargs$id <- "(unspecified), which defaults to 1:length(Y)"
    if (fnargs$pair[[length(fnargs$pair)]] != fnargs$id[[length(fnargs$id)]]) {
      stop(paste("\nBy specifying the pair argument, you have indicated this is a pair-matched analysis.\n\nTo get correct variance, the variable you pass to id\nmust be the same as the variable you pass to pair. \n\nIf you are not doing a pair-matched analysis, then specify pair=NULL\n(or leave it unspecified).\n\n   You specified pair=", 
                 fnargs$pair[[length(fnargs$pair)]], "and id=", 
                 fnargs$id[[length(fnargs$id)]]))
    }
  }
  if (!is.null(FECR)) {
    if (FECR != "arithmetic" & FECR != "geometric") {
      stop(paste("You specified FECR=", fnargs$FECR[[length(fnargs$FECR)]], 
                 "to estimate the fecal egg count reduction %\nYou need to supply either 'arithmetic' or 'geometric' as an argument to the FECR option."))
    }
    if ((FECR == "arithmetic" | FECR == "geometric") & 
        family != "gaussian") {
      stop(paste("You specified FECR=", fnargs$FECR[[length(fnargs$FECR)]], 
                 "to estimate the fecal egg count reduction %\nThis parameter is a ratio of means: FECR=(EY1/EY0)-1\nso you need to specify family='gaussian' to estimate it properly."))
    }
  }
  if (is.null(W)) {
    if (is.null(pair)) {
      tmledat <- data.frame(id, Y, Delta, tr)
    }
    else {
      tmledat <- data.frame(id, pair, Y, Delta, tr)
    }
  }
  else {
    if (is.null(pair)) {
      tmledat <- data.frame(id, Y, Delta, tr, W)
    }
    else {
      tmledat <- data.frame(id, pair, Y, Delta, tr, W)
    }
  }
  tmledat <- subset(tmledat, tr == contrast[1] | tr == contrast[2])
  tmledat$tr <- factor(tmledat$tr, levels = contrast[1:2])
  if (!is.null(pair)) {
    n.orig <- dim(tmledat)[1]
    miss <- NULL
    activeOnly <- ((subset(tmledat, tr == contrast[1])))
    nomiss <- sort(unique(activeOnly$pair))
    miss1 <- (unique(pair)[which(!(unique(pair) %in% (nomiss)))])
    activeOnly2 <- ((subset(tmledat, tr == contrast[2])))
    nomiss2 <- sort(unique(activeOnly2$pair))
    miss2 <- (unique(pair)[which(!(unique(pair) %in% (nomiss2)))])
    miss <- append(miss1, miss2)
    tmledat <- subset(tmledat, !(pair %in% miss))
    n.sub <- dim(tmledat)[1]
    if ((print == TRUE) & (n.orig > n.sub)) {
      cat("\n-----------------------------------------", 
          "\nThere were", length(unique(miss)), "pairs dropped because they were\nmissing at least one treatment level.\nThis is the list of their IDs:\n", 
          sort(unique(miss)))
      cat("\n-----------------------------------------", 
          "\nStarting N:  ", n.orig, "\nN after dropping incomplete blocks: ", 
          n.sub, "\n\nTotal of", n.orig - n.sub, 
          "observations dropped\n because of unmatched pairs.", 
          "\n-----------------------------------------\n")
    }
  }
  n_orig <- dim(tmledat)[1]
  tmledat <- tmledat[complete.cases(tmledat), ]
  n_sub <- dim(tmledat)[1]
  if ((print == TRUE) & (n_orig > n_sub)) {
    cat("\n-----------------------------------------\nTotal of", 
        n_orig - n_sub, "observations dropped due to missing\nvalues in one or more variables\n", 
        " Final sample size:", n_sub, "\n-----------------------------------------\n")
  }
  if (!is.null(W)) {
    if (print == TRUE) {
      cat("\n-----------------------------------------\nPre-screening the adjustment covariates\nusing a univariate liklihood ratio test:\n-----------------------------------------\n")
    }
    Wscreen <- washb_prescreen(Y = tmledat$Y[tmledat$Delta == 
                                               1], Ws = subset(tmledat, tmledat$Delta == 1, select = Wvars), 
                               family = family, pval = pval, print = print)
    if (print == TRUE) {
      cat("\n-----------------------------------------\n")
    }
    if (length(Wscreen) > 0) {
      Wselect <- subset(tmledat, select = Wscreen)
      Wselect <- design_matrix(Wselect)
      if (ncol(Wselect) <= 1 & (length(grep("SL.glmnet", 
                                            Q.SL.library)) + length(grep("SL.glmnet", 
                                                                         g.SL.library)) > 0)) {
        cat("\nNOTE: Dropping SL.glmnet from the library\nbecause there is only 1 covariate selected and glmnet\nrequires 2+ covariates to run\n")
        Q.SL.library <- Q.SL.library[-grep("SL.glmnet", 
                                           Q.SL.library)]
        g.SL.library <- g.SL.library[-grep("SL.glmnet", 
                                           g.SL.library)]
      }
    }
    else {
      cat("\n\nSince no covariates were associated with the outcome,\nthe estimates below are unadjusted...")
      if (n_orig > n_sub) {
        cat("\n\nIn this case, since", n_orig - 
              n_sub, "observations were dropped\ndue to missing covariates,\nthose observations were not included in the analysis.\nIt would be best to re-estimate the effect without covariates (W=NULL)\nto avoid unnecessarily dropping these observations.")
      }
      W <- NULL
    }
  }
  if (is.null(W)) {
    Wselect <- data.frame(w1 = runif(n = length(tmledat$Y)))
    Q.SL.library <- c("SL.glm")
    g.SL.library <- c("SL.glm")
  }
  tmle_Y <- tmledat$Y
  tmle_A <- ifelse(tmledat$tr == contrast[2], 1, 0)
  tmle_Delta <- tmledat$Delta
  tmle_id <- as.numeric(tmledat$id)
  if (!is.null(seed)) 
    set.seed(seed)
  tmle_fit <- tmle(Y = tmle_Y, A = tmle_A, W = Wselect, Delta = tmle_Delta, 
                   id = tmle_id, family = family, 
                   Q.SL.library = Q.SL.library, 
                   g.Delta.SL.library =  g.SL.library,
                   g.SL.library = g.SL.library)
  
  if (print == TRUE) {
    cat("\n-----------------------------------------\nEstimation Results:\n-----------------------------------------\n")
    print(summary(tmle_fit))
    cat("\n-----------------------------------------\n")
  }
  if (!is.null(FECR)) {
    if (print == TRUE) {
      cat(paste("\n-----------------------------------------\nEstimating the fecal egg count reduction\n(FECR) proportion = (EY1/EY0) - 1\nfrom TMLE results using", 
                FECR, "means\nand the delta method (for a ratio of means)\n-----------------------------------------\n"))
    }
    n_id <- length(unique(tmle_id))
    pDelta0 <- tmle_fit$g.Delta$g1W[, 1]
    pDelta1 <- tmle_fit$g.Delta$g1W[, 2]
    g1 <- tmle_fit$g$g1W
    Qst0 <- tmle_fit$Qstar[, 1]
    Qst1 <- tmle_fit$Qstar[, 2]
    Ey0 <- mean(Qst0)
    Ey1 <- mean(Qst1)
    IC0 <- (tmle_Delta/pDelta0) * ((1 - tmle_A)/(1 - g1)) * 
      (tmle_Y - Qst0) + Qst0 - Ey0
    IC1 <- (tmle_Delta/pDelta1) * (tmle_A/g1) * (tmle_Y - 
                                                   Qst1) + Qst1 - Ey1
    if (n_id < length(id)) {
      IC0 <- as.vector(by(IC0, tmle_id, mean))
      IC1 <- as.vector(by(IC1, tmle_id, mean))
    }
    vc <- (1/n_id) * var(cbind(IC0, IC1))
    if (FECR == "arithmetic") {
      fecr <- (Ey1/Ey0) - 1
      fderiv <- c(-Ey1/(Ey0^2), 1/Ey0)
    }
    if (FECR == "geometric") {
      fecr <- (exp(Ey1)/exp(Ey0)) - 1
      fderiv <- c(-exp(Ey1)/exp(Ey0), exp(Ey1)/exp(Ey0))
    }
    fecr_se <- as.vector(sqrt(t(fderiv) %*% vc %*% fderiv))
    fecr_lb <- fecr - 1.96 * fecr_se
    fecr_ub <- fecr + 1.96 * fecr_se
    fecr_p <- 2 * (1 - pnorm(abs(fecr/fecr_se)))
    tmle_fit$estimates$FECR <- list(psi = fecr, var.psi = fecr_se^2, 
                                    CI = c(fecr_lb, fecr_ub), pvalue = fecr_p, method = FECR)
    if (print == TRUE) {
      cat(paste("\n-----------------------------------------\nFecal egg count reduction (EY1/EY0)-1,\nestimated using ", 
                FECR, " means", "\n-----------------------------------------\n", 
                sep = ""))
      cat(paste("FECR (95% CI) : ", sprintf("%1.3f", 
                                            fecr), " (", sprintf("%1.3f", fecr_lb), 
                ", ", sprintf("%1.3f", fecr_ub), 
                ")", sep = ""))
      cat("\n     SE(FECR) :", sprintf("%1.4f", 
                                       fecr_se))
      cat("\n      p-value :", sprintf("%1.4f", 
                                       fecr_p))
      cat("\n-----------------------------------------\n")
    }
  }
  return(tmle_fit)
}