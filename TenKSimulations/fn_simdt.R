##
## script that simulates 10,000 pairs of ipd's
##
require(dplyr, quietly = TRUE)
require(tidyverse, quietly = TRUE)
## for maicChecks usage see section on Exact Matching on https://clyau.github.io/maicChecks/
require(maicChecks) ## version 0.2.0
## 
source("funs.R") 
## common parameter values
n1_obs <- 300 ## number of observations in ipd1
n2_obs <- 300 ## number of observations in ipd2
# Set the number of variables and the block size
n_variables <- 15
block_size <- 5
# correlation within each block
rho <- c(0.3, 0.5, 0.7)
## correlation matrix for the 15 covariates
corr_mat <- diag(1, n_variables, n_variables)
for (i in 1:(n_variables/block_size)) {
  start_index <- (i-1) * block_size + 1
  end_index <- i * block_size
  corr_mat[start_index:end_index, 
           start_index:end_index] <- rho[i]
}
diag(corr_mat) <- 1 ## correct diagonal values
##
## random choice of some variables to categorize
vars_to_categorize <- c("X1", "X2", "X3", ## in block 1
                        "X6", "X7", "X8", ## in block 2
                        "X11", "X12", "X13", "X14" ## in block 3
)
## some random choices of quantiles from standard normal for 
## ... categorization of selected variables
vc_co <- list(
  q.x1 = qnorm(0.238),                    ## x1, binary
  q.x2 = qnorm(0.312),                    ## x2, binary
  q.x3 = qnorm(c(0.12, 0.335, 0.68)),     ## x3, 4 levels
  q.x6 = qnorm(0.439),                    ## x6, binary
  q.x7 = qnorm(0.581),                    ## x7, binary
  q.x8 = qnorm(c(0.23, 0.56)),            ## x8, 3 levels
  q.x11 = qnorm(0.607),                   ## x11, binary
  q.x12 = qnorm(0.712),                   ## x12, binary
  q.x13 = qnorm(0.842),                   ## x13, binary
  q.x14 = qnorm(c(0.18, 0.3, 0.56, 0.72)) ## x14, 5 levels
)
##
## function that simulates and match one pair of ipd data ::::::::::::::: 
##
sim.dt.fn <- function(r_and_s, n1_obs, n2_obs) {
  
  ##
  ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
  ## :::::::::::::::::::: simulate data :::::::::::::::::::: ##
  ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
  ## 
  ## IPD A ::::::::::::::::::::::::::::::::::::::::::::::::: ##
  ## 
  
  set.seed(r_and_s)
  
  simdt1 <- mvtnorm::rmvnorm(
    n = n1_obs, 
    mean = rep(0, n_variables), 
    sigma = corr_mat
  )
  
  colnames(simdt1) <- paste0('X', 1:15)
  simdt1 <- as.data.frame(simdt1) 
  
  ## simulate response variable Y
  set.seed(r_and_s + 100)
  
  ## simulate Y to depend on X1, X3, X8, X11 (all categorical)
  ## ... and X9 and X15 (numerical)
  simdt1 <- simdt1 %>%
    mutate(Y = 0.3 * (X1 > vc_co$q.x1[1]) + 
             0.2 * (X3 > vc_co$q.x3[1]) + 
             0.2 * (X3 > vc_co$q.x3[2]) + 
             0.2 * (X3 > vc_co$q.x3[3]) +
             0.15 * (X8 > vc_co$q.x8[1]) + 
             0.15 * (X8 > vc_co$q.x8[2]) +
             0.1 * X9 + 
             0.2 * (X11 > vc_co$q.x11[1]) + 
             0.1 * X15 ,
           Y = Y + rnorm(n = n1_obs, mean = 0, sd = 1)
    ) 
  
  ## simulate another version of the response variable with the same seed
  set.seed(r_and_s + 100)
  
  ## simulate Yc to depend on X1, X3, X8, X9, X11, and X15 (all numerical)
  simdt1 <- simdt1 %>%
    mutate(Yc = 0.3 * X1 + 
             0.2 * X3 +
             0.3 * X8 +
             0.1 * X9 + 
             0.2 * X11 + 
             0.1 * X15 ,
           Yc = Yc + rnorm(n = n1_obs, mean = 0, sd = 1)
    ) %>%
    ## convert selected variables to categorical scale 
    num2cat(df = ., 
            vars_to_cat = vars_to_categorize, 
            thresholds = vc_co)
  
  ## 
  ## IPD B :::::::::::::::::::::::::::::::::::::::::::::::::: ##
  ## 
  
  set.seed(r_and_s + 10)
  
  simdt2 <- mvtnorm::rmvnorm(
    n = n2_obs, 
    ## shift means for 3 variables
    mean = c(1, 0, 0, 0, 0, ## block 1, x1 shifted
             0, 0, 1, 0, 0, ## block 2, x8 shifted
             0, 0, 0, 0, 1  ## block 3, x15 shifted
    ),
    sigma = corr_mat
  )
  
  colnames(simdt2) <- paste0('X', 1:15)
  simdt2 <- as.data.frame(simdt2) ## convert to dataframe
  
  ## simulate response variable
  set.seed(r_and_s + 10 + 100)
  
  ## same dependency on X's as simdt1 
  simdt2 <- simdt2 %>%
    mutate(Y = 0.3 * (X1 > vc_co$q.x1[1]) + 
             0.2 * (X3 > vc_co$q.x3[1]) + 
             0.2 * (X3 > vc_co$q.x3[2]) + 
             0.2 * (X3 > vc_co$q.x3[3]) +
             0.15 * (X8 > vc_co$q.x8[1]) + 
             0.15 * (X8 > vc_co$q.x8[2]) +
             0.1 * X9 + 
             0.2 * (X11 > vc_co$q.x11[1]) + 
             0.1 * X15 ,
           Y = Y + rnorm(n = n2_obs, mean = 0, sd = 1)
    ) 
  
  ## simulate another version of the response variable with the same seed
  set.seed(r_and_s + 10 + 100)
  
  ## simulate Yc to depend on X1, X3, X8, X9, X11, and X15 (all numerical)
  simdt2 <- simdt2 %>%
    mutate(Yc = 0.3 * X1 + 
             0.2 * X3 +
             0.3 * X8 +
             0.1 * X9 + 
             0.2 * X11 + 
             0.1 * X15 ,
           Yc = Yc + rnorm(n = n2_obs, mean = 0, sd = 1)
    ) %>%
    ## convert selected variables to categorical scale 
    num2cat(df = ., 
            vars_to_cat = vars_to_categorize, 
            thresholds = vc_co)
  
  ##
  ## combine datasets :::::::::::::::::::::::::::::::::::::: ##
  ## 
  
  simdt <- data.frame(rbind(simdt1, simdt2)) %>%
    mutate(study = c(rep('IPD A', n1_obs), rep('IPD B', n2_obs)),
           id = sprintf("%03d", (1 : (n1_obs + n2_obs))),
           ipd1 = (study == 'IPD A')) ## for logistic regression for PS weights
  
  ##
  ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
  ## ::::::::::: matching & calculating weights :::::::::::: ##
  ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
  ##
  ## exact matching :::::::::::::::::::::::::::::::::::::::: ##
  ##
  ##
  ## check if matching can be done for unconstrained weights
  ##
  checks <- maicChecks::exmLP.2ipd(simdt[simdt$study == 'IPD A',], 
                             simdt[simdt$study == 'IPD B',],
                             vars_to_match = paste0('X', 1:15),
                             cat_vars_to_01 = vars_to_categorize,
                             mean.constrained = FALSE
  )
  if(checks$lp.check != 0) {
    ## this checks whether the *unconstrained* solution exists
    print(paste0("LP checks for unconstrained weights = ", 
                 checks$lp.check, "rand.seed = ", r_and_s))
  }
  start_time_proc <- proc.time()
  x.unc <- tryCatch({ ## in case the linear check does not catch it
    maicChecks::exmWt.2ipd(simdt[simdt$study == 'IPD A',], 
                     simdt[simdt$study == 'IPD B',],
                     vars_to_match = paste0('X', 1:15),
                     cat_vars_to_01 = vars_to_categorize,
                     mean.constrained = FALSE)
  }, error = function(e) {
    ## manually create dummy data if no solution
    ipd1 <- cat201(simdt[simdt$ipd1, paste0('X', 1:15)], vars_to_categorize)
    ipd1$study <- 'IPD A'
    ipd1$exm.wts <- NA
    ipd2 <- cat201(simdt[!simdt$ipd1, paste0('X', 1:15)], vars_to_categorize)
    ipd2$study <- 'IPD B'
    ipd2$exm.wts <- NA
    list(ipd1 = ipd1,
         ipd2 = ipd2, 
         wtd.summ = NA
    )
  })
  stop_time_proc <- proc.time()
  run_time_unc <- stop_time_proc - start_time_proc
  ##
  ## check if matching can be done for constrained weights
  ##
  checks <- maicChecks::exmLP.2ipd(simdt[simdt$study == 'IPD A',], 
                             simdt[simdt$study == 'IPD B',],
                             vars_to_match = paste0('X', 1:15),
                             cat_vars_to_01 = vars_to_categorize,
                             mean.constrained = TRUE
  )
  ##
  if(checks$lp.check != 0) {
    ## this checks whether the *constrained* solution exists
    print(paste0("LP checks for constrained weights = ", 
                 checks$lp.check, ", rand.seed = ", r_and_s))
  }
  ##
  start_time_proc <- proc.time()
  x.con <- tryCatch({ ## in case the linear check does not catch it
    ## this checks if *constrained* solution exists
    maicChecks::exmWt.2ipd(simdt[simdt$study == 'IPD A',], 
                     simdt[simdt$study == 'IPD B',],
                     vars_to_match = paste0('X', 1:15),
                     cat_vars_to_01 = vars_to_categorize,
                     mean.constrained = TRUE)
  }, error = function(e) {
    ## manually create dummy data if no solution
    ipd1 <- cat201(simdt[simdt$ipd1, paste0('X', 1:15)], vars_to_categorize)
    ipd1$study <- 'IPD A'
    ipd1$ipd1 <- TRUE
    ipd1$exm.wts <- NA
    ipd2 <- cat201(simdt[!simdt$ipd1, paste0('X', 1:15)], vars_to_categorize)
    ipd2$study <- 'IPD B'
    ipd2$ipd1 <- FALSE
    ipd2$exm.wts <- NA
    list(ipd1 = ipd1,
         ipd2 = ipd2, 
         wtd.summ = NA
    )
  })
  stop_time_proc <- proc.time()
  run_time_con <- stop_time_proc - start_time_proc
  
  ##
  ## propensity score matching :::::::::::::::::::::::::::::: ##
  ##
 
  ## ps matching
  psmatch <- glm(ipd1 ~ X4 + X5 + X9 + X10 + X15 +
                   factor(X1) + factor(X2) + factor(X3) +
                   factor(X6) + factor(X7) + factor(X8) + 
                   factor(X11) + factor(X12) + factor(X13) + factor(X14),
                 data = simdt,
                 family = binomial())
  
  ## psvalue is the propensity score
  simdt$psvalue <- predict(psmatch, type = "response")
  
  ##
  ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
  ## :::::::: putting pieces together for summary :::::::::: ##
  ## ::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
  ##
  ## adding weights to original dataset :::::::::::::::::::: ##
  ##  
  
  simdt <- simdt %>%
    mutate(
      ## ps weights
      wt.ps = ifelse(study == 'IPD A', 
                     1 / psvalue,      ## ipd1
                     1 / (1 - psvalue) ## ipd2
      ),
      ## unconstrained exact matching weights
      wt.unc = c(x.unc$ipd1$exm.wts,  ## ipd1
                 x.unc$ipd2$exm.wts   ## ipd2
      ),
      ## constrained exact matching weights
      wt.con = c(x.con$ipd1$exm.wts,  ## ipd1
                 x.con$ipd2$exm.wts   ## ipd2
      )
    ) %>%
    ## standardize all the weights so they can be compared
    group_by(study) %>%
    mutate(
      wt.unc.sdz = ifelse(is.na(wt.unc), NA, 100 * wt.unc / sum(wt.unc)),
      wt.con.sdz = ifelse(is.na(wt.con), NA, 100 * wt.con / sum(wt.con)),
      wt.ps.sdz = 100 * wt.ps / sum(wt.ps)
    ) %>%
    ungroup()
  
  ## 
  ## observed and weighted means of the variables used in matching ::::: ##
  ##
  
  ipd <- data.frame(rbind(x.con$ipd1, x.con$ipd2))
  ## v_names contains variables (including indicator names ...
  ## ... for categorical variables) used in exact matching
  v_names <- names(ipd)[!(names(ipd) %in% 
                                   c('study', 'ipd1', 'id', 'exm.wts', 'Y', 'Yc'))]
  ## weighted mean and variance for ps weighted covariates
  ## ... variance (SD) needed to calculate SMD
  ipd1_ps_wtmean <- sapply(ipd[ipd$study == 'IPD A', v_names],
                           function(x) 
                             Hmisc::wtd.mean(
                               x, 
                               simdt[which(ipd$study == 'IPD A'), ]$wt.ps.sdz,
                               normwt = TRUE
                             ))
  ipd2_ps_wtmean <- sapply(ipd[ipd$study == 'IPD B', v_names], 
                           function(x) 
                             Hmisc::wtd.mean(
                               x, 
                               simdt[ipd$study == 'IPD B', ]$wt.ps.sdz,
                               normwt = TRUE
                             ))
  ipd1_ps_wtvar <- sapply(ipd[ipd$study == 'IPD A', v_names], 
                          function(x)
                            Hmisc::wtd.var(
                              x,
                              simdt[ipd$study == 'IPD A', ]$wt.ps.sdz,
                              normwt = TRUE,
                              method = 'ML'
                            ))
  ipd2_ps_wtvar <- sapply(ipd[ipd$study == 'IPD B', v_names], 
                          function(x)
                            Hmisc::wtd.var(
                              x,
                              simdt[ipd$study == 'IPD B', ]$wt.ps.sdz,
                              normwt = TRUE,
                              method = 'ML'
                            ))
  ## pooled standard deviation of the two ipd's
  ps_wtsd_pooled <- sqrt((n1_obs * ipd1_ps_wtvar + 
                            n2_obs * ipd2_ps_wtvar) / (n1_obs + n2_obs))
  ## weighted sample size, or ESS of ps matching
  ipd1_ps_ess <- round(sum(simdt[ipd$study == 'IPD A', ]$wt.ps.sdz) ^ 2 / 
                         sum(simdt[ipd$study == 'IPD A', ]$wt.ps.sdz ^ 2), 1)
  ipd2_ps_ess <- round(sum(simdt[ipd$study == 'IPD B', ]$wt.ps.sdz) ^ 2 / 
                         sum(simdt[ipd$study == 'IPD B', ]$wt.ps.sdz ^ 2), 1)
  
  ##
  ## putting all the means (observed and weighted) into a dataframe for reporting
  ##
  
  ## unconstrained matching means (same for ipd1 and ipd2)
  ## ... unconstrained matching uses full-rank design matrix for
  ## ... categorical variables, different from constrained matching
  ## ... therefore, need to be handled separately
  unc_means <- data.frame(ipd12.unc = t(x.unc[['wtd.summ']])) %>%
    mutate(v.names = rownames(.)) %>%
    filter(!(v.names %in% c('ipd1.ess', 'ipd2.ess')))

  all_means <- data.frame(
    ## observed means
    ipd1.obs = apply(ipd[ipd$study == 'IPD A', v_names], 2, mean),
    ipd2.obs = apply(ipd[ipd$study == 'IPD B', v_names], 2, mean),
    ## constraint exact weighted means of the covariates (same for ipd1 and ipd2)
    ipd12.con = x.con$wtd.summ[3:(length(x.con$wtd.summ))],
    ## ps weighted means
    ipd1.ps = ipd1_ps_wtmean, ## ps weighted means
    ipd2.ps = ipd2_ps_wtmean  ## ps weighted means
  ) %>%
    mutate(v.names = rownames(.)) %>%
    left_join(unc_means, by = 'v.names') %>%
    select(v.names, ipd1.obs, ipd2.obs, ipd12.unc, ipd12.con, ipd1.ps, ipd2.ps)
  
  ## smd (standardized mean difference) for ps matching
  all_means$abs.smd <- abs(all_means$ipd1.ps - all_means$ipd2.ps) / ps_wtsd_pooled
  
  ## ess
  ess <- data.frame(
    ess.ipd1.unc = x.unc$wtd.summ[1],
    ess.ipd2.unc = x.unc$wtd.summ[2],
    ess.ipd1.con = x.con$wtd.summ[1],
    ess.ipd2.con = x.con$wtd.summ[2],
    ess.ipd1.ps = ipd1_ps_ess,
    ess.ipd2.ps = ipd2_ps_ess
  )
  
  ## difference in response variable Y
  ybar <- data.frame(ipd1.obs = mean(simdt1$Y),
                     ipd2.obs = mean(simdt2$Y),
                     
                     ipd1.unc = Hmisc::wtd.mean(simdt1$Y, 
                                                simdt[ipd$study == 'IPD A', ]$wt.unc.sdz),
                     ipd2.unc = Hmisc::wtd.mean(simdt2$Y, 
                                                simdt[ipd$study == 'IPD B', ]$wt.unc.sdz),
                     
                     ipd1.con = Hmisc::wtd.mean(simdt1$Y, 
                                                simdt[ipd$study == 'IPD A', ]$wt.con.sdz),
                     ipd2.con = Hmisc::wtd.mean(simdt2$Y, 
                                                simdt[ipd$study == 'IPD B', ]$wt.con.sdz),
                     
                     ipd1.ps = Hmisc::wtd.mean(simdt1$Y, 
                                               simdt[ipd$study == 'IPD A', ]$wt.ps.sdz),
                     ipd2.ps = Hmisc::wtd.mean(simdt2$Y, 
                                               simdt[ipd$study == 'IPD B', ]$wt.ps.sdz)
  ) %>%
    mutate(diff.obs = ipd1.obs - ipd2.obs,
           diff.unc = ipd1.unc - ipd2.unc,
           diff.con = ipd1.con - ipd2.con,
           diff.ps = ipd1.ps - ipd2.ps)
  
  ## difference in response variable Yc
  ycbar <- data.frame(ipd1.obs = mean(simdt1$Yc),
                     ipd2.obs = mean(simdt2$Yc),
                     
                     ipd1.unc = Hmisc::wtd.mean(simdt1$Yc, 
                                                simdt[ipd$study == 'IPD A', ]$wt.unc.sdz),
                     ipd2.unc = Hmisc::wtd.mean(simdt2$Yc, 
                                                simdt[ipd$study == 'IPD B', ]$wt.unc.sdz),
                     
                     ipd1.con = Hmisc::wtd.mean(simdt1$Yc, 
                                                simdt[ipd$study == 'IPD A', ]$wt.con.sdz),
                     ipd2.con = Hmisc::wtd.mean(simdt2$Yc, 
                                                simdt[ipd$study == 'IPD B', ]$wt.con.sdz),
                     
                     ipd1.ps = Hmisc::wtd.mean(simdt1$Yc, 
                                               simdt[ipd$study == 'IPD A', ]$wt.ps.sdz),
                     ipd2.ps = Hmisc::wtd.mean(simdt2$Yc, 
                                               simdt[ipd$study == 'IPD B', ]$wt.ps.sdz)
  ) %>%
    mutate(diff.obs = ipd1.obs - ipd2.obs,
           diff.unc = ipd1.unc - ipd2.unc,
           diff.con = ipd1.con - ipd2.con,
           diff.ps = ipd1.ps - ipd2.ps)
  
  return(list(simdt = data.frame(simdt), 
              ess = data.frame(ess),
              all_means = data.frame(all_means), 
              ybar = ybar,
              ycbar = ycbar,
              run_time = data.frame(cbind(run_time_unc = run_time_unc[1:3], 
                                          run_time_con = run_time_con[1:3]))))
}
##
## :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ##
