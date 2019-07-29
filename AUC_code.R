#--- Setting the current working directory to source file location ---#

if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

compute_AUC = function(cohort.data, case.control.data, full.model.formula, 
                       reduced.model.formula, risk.factors.logRR, 
                       no.riskscore.cat){
  
  N1 = sum(cohort.data$event)
  N0 = dim(cohort.data)[1] - N1
  
  n1 = sum(case.control.data$event)
  n0 = sum(case.control.data$event == 0)
  
  case.status = cohort.data$event
  
  cohort.covariate.data = model.matrix(reduced.model.formula, 
                                       cohort.data[,all.vars(reduced.model.formula)])[,-1]
  
  cohort.logRR = risk.factors.logRR[colnames(cohort.covariate.data)]
  
  obs.riskscore = as.numeric(cohort.covariate.data %*% cohort.logRR)
  
  obs.riskscore.cat.cases = quantcut(obs.riskscore[case.status == 1],
                                     q = no.riskscore.cat,na.rm = TRUE)
  obs.riskscore.cat.controls = quantcut(obs.riskscore[case.status == 0],
                                        q = no.riskscore.cat,na.rm = TRUE)
  levels(obs.riskscore.cat.cases) = 
    levels(obs.riskscore.cat.controls) = 1:no.riskscore.cat
  
  cohort.data$obs.riskscore = obs.riskscore
  
  cohort.data$obs.riskscore.cat.cases[cohort.data$event == 1] = obs.riskscore.cat.cases
  cohort.data$obs.riskscore.cat.controls[cohort.data$event == 0] = obs.riskscore.cat.controls
  
  case.control.covariate.data = model.matrix(full.model.formula, 
                                             case.control.data[,all.vars(full.model.formula)])[,-1]
  
  case.control.data$riskscore = as.numeric(case.control.covariate.data 
                                           %*% risk.factors.logRR)
  
  case.control.data$obs.riskscore.cat.cases = 
    cohort.data$obs.riskscore.cat.cases[which(is.element(cohort.data$id,case.control.data$id))]
  case.control.data$obs.riskscore.cat.controls = 
    cohort.data$obs.riskscore.cat.controls[which(is.element(cohort.data$id,case.control.data$id))]
  
  cc.cases = case.control.data[case.control.data$event == 1,]
  cc.controls = case.control.data[case.control.data$event == 0,]
  
  cohort.cases = cohort.data[cohort.data$event == 1,]
  cohort.controls = cohort.data[cohort.data$event == 0,]
  
  obs.riskscore.cases = cohort.cases$obs.riskscore
  obs.riskscore.controls = cohort.controls$obs.riskscore
  
  indicator = vapply(obs.riskscore.cases, 
                     function(x) x > obs.riskscore.controls,logical(length(obs.riskscore.controls))) * 1
  
  auc_cohort = mean(indicator)
  
  #compute p_0(S_1) = E_S_0[I(S_1 > S_0)]
  mean_S0_indicator = apply(indicator,2,mean)
  
  #compute E_S_1[I(S_1 > S_0)]
  mean_S1_indicator = apply(indicator,1,mean)
  
  var_auc_cohort = (mean((mean_S0_indicator - auc_cohort)^2))/N1 + 
    (mean((mean_S1_indicator - auc_cohort)^2))/N0
  
  #Compute the conditional expectation matrix for the categories
  
  #Break up the cohort into 4 groups: (1) case-control pair both included, 
  #(2) case included control excluded, (3) case excluded control included, 
  #(4) case-control pair excluded
  
  #Computing \psi_1(S1) and \psi_0(S0)
  
  psi1_by_wt = rep(NA,n1)
  psi0_by_wt = rep(NA,n0)
  
  #Compute E[I(S1>S0)|S1,S0obs]
  cond.exp2 = matrix(0,no.riskscore.cat,n1)
  
  #Compute E[I(S1>S0) | S1obs,S0]
  cond.exp3 = matrix(0,n0,no.riskscore.cat)
  
  for(i in 1:no.riskscore.cat){
    
    riskscore.cases = cc.cases$riskscore[cc.cases$obs.riskscore.cat.cases == i]
    riskscore.controls = cc.controls$riskscore
    
    indicator = vapply(riskscore.cases, 
                       function(x) x > riskscore.controls,logical(length(riskscore.controls))) * 1
    
    sampling.weights.cases = cc.cases$sampling.weights[cc.cases$obs.riskscore.cat.cases == i]
    sampling.weights.controls = cc.controls$sampling.weights
    
    freq.cases = 1/sampling.weights.cases
    freq.controls = 1/sampling.weights.controls
    
    weight.mat = matrix(kronecker(freq.controls,freq.cases), 
                        nrow = length(freq.controls), byrow = TRUE)
    
    indicator.times.one_minus_casewt = t(indicator) * (1 - sampling.weights.cases)
    
    first_term = weighted.mean((1 - sampling.weights.cases),w = freq.cases) * 
      apply(indicator,2,weighted.mean,w = freq.controls)
    
    second_term = sum(t(indicator.times.one_minus_casewt) * weight.mat)/sum(weight.mat)
    
    psi1_by_wt[cc.cases$obs.riskscore.cat.cases == i] = 
      (first_term - second_term)/sampling.weights.cases
    
    
    #Compute E[I(S1>S0) | S1obs,S0]
    cond.exp3[,i] = apply(indicator,1,weighted.mean,w = freq.cases)
    
    riskscore.cases = cc.cases$riskscore
    riskscore.controls = cc.controls$riskscore[cc.controls$obs.riskscore.cat.controls == i]
    
    indicator = vapply(riskscore.cases, function(x) x > riskscore.controls,logical(length(riskscore.controls))) * 1
    
    sampling.weights.cases = cc.cases$sampling.weights
    sampling.weights.controls = cc.controls$sampling.weights[cc.controls$obs.riskscore.cat.controls == i]
    
    freq.cases = 1/sampling.weights.cases
    freq.controls = 1/sampling.weights.controls
    
    weight.mat = matrix(kronecker(freq.controls,freq.cases), nrow = length(freq.controls), byrow = TRUE)
    
    indicator.times.one_minus_controlwt = indicator * (1 - sampling.weights.controls)
    
    first_term = weighted.mean((1 - sampling.weights.controls),w = freq.controls) * 
       apply(indicator,1,weighted.mean,w = freq.cases)
    
    second_term = sum(indicator.times.one_minus_controlwt * weight.mat)/sum(weight.mat)
    
    psi0_by_wt[cc.controls$obs.riskscore.cat.controls == i] = (first_term - second_term)/sampling.weights.controls
    
    #Compute E[I(S1>S0)|S1,S0obs]
    cond.exp2[i,] = apply(indicator,2,weighted.mean,w = freq.controls)
    
  }
  
  #Group 1: 
  
  riskscore.cases = cc.cases$riskscore
  riskscore.controls = cc.controls$riskscore
  
  indicator = vapply(riskscore.cases, 
                     function(x) x > riskscore.controls,logical(length(riskscore.controls))) * 1
  
  auc_term1 = sum(indicator)
  
  sampling.weights.cases = cc.cases$sampling.weights
  sampling.weights.controls = cc.controls$sampling.weights
  
  freq.cases = 1/sampling.weights.cases
  freq.controls = 1/sampling.weights.controls
  
  cases.obs.riskscore.cat = cc.cases$obs.riskscore.cat.cases
  controls.obs.riskscore.cat = cc.controls$obs.riskscore.cat.controls
  
  #compute p_0(S_1) = E_S_0[I(S_1 > S_0)]
  mean_S0_indicator = apply(indicator,2,weighted.mean,w = freq.controls)
  
  #compute E_S_1[I(S_1 > S_0)]
  mean_S1_indicator = apply(indicator,1,weighted.mean,w = freq.cases)
  
  ##########[This piece of code also computes the IPW estimator for AUC and its variance]##############
  
  weight.mat = matrix(kronecker(freq.controls,freq.cases), 
                      nrow = length(freq.controls), byrow = TRUE)
  auc_ipw = sum(indicator * weight.mat)/sum(weight.mat)
  
  var_auc_ipw = (weighted.mean(((mean_S0_indicator - auc_ipw)^2/sampling.weights.cases), 
                               w = freq.cases))/sum(freq.cases) + 
    (weighted.mean(((mean_S1_indicator - auc_ipw)^2/sampling.weights.controls), 
                   w = freq.controls))/sum(freq.controls)
  
  
  ####################################################################################################
  
  
  #compute p(S_1_obs) and p(S_0_obs)
  
  mean_S0_indicator.wtcases = mean_S0_indicator * freq.cases
  
  mean_S0_indicator.wtcases.cat = aggregate(mean_S0_indicator.wtcases ~ cases.obs.riskscore.cat, 
                                            FUN = sum)
  wtcases.cat = aggregate(freq.cases ~ cases.obs.riskscore.cat, FUN = sum)
  
  mean_S1S0.S1obs_indicator = mean_S0_indicator.wtcases.cat$mean_S0_indicator.wtcases/wtcases.cat$freq.cases
  
  mean_S1_indicator.wtcontrols = mean_S1_indicator * freq.controls
  
  mean_S1_indicator.wtcontrols.cat = aggregate(mean_S1_indicator.wtcontrols ~ controls.obs.riskscore.cat, FUN = sum)
  wtcontrols.cat = aggregate(freq.controls ~ controls.obs.riskscore.cat, FUN = sum)
  
  mean_S1S0.S0obs_indicator = mean_S1_indicator.wtcontrols.cat$mean_S1_indicator.wtcontrols/wtcontrols.cat$freq.controls
  
  p_S1obs = mean_S1S0.S1obs_indicator[cohort.data$obs.riskscore.cat.cases[cohort.data$event == 1]]
  
  p_S0obs = mean_S1S0.S0obs_indicator[cohort.data$obs.riskscore.cat.controls[cohort.data$event == 0]]
  
  p0_S1 = psi1.by.wt = rep(0,N1)
  p0_S1[which(is.element(cohort.cases$id,cc.cases$id))] = mean_S0_indicator
  
  psi1.by.wt[which(is.element(cohort.cases$id,cc.cases$id))] = psi1_by_wt
  
  p1_S0 = psi0.by.wt = rep(0,N0)
  p1_S0[which(is.element(cohort.controls$id,cc.controls$id))] = mean_S1_indicator
  
  psi0.by.wt[which(is.element(cohort.controls$id,cc.controls$id))] = psi0_by_wt
  
  case.included = cohort.data$include[cohort.data$event == 1]
  control.included = cohort.data$include[cohort.data$event == 0]
  
  cohort.case.weights = cohort.data$sampling.weights[cohort.data$event == 1]
  cohort.control.weights = cohort.data$sampling.weights[cohort.data$event == 0]
  
  cohort2 = cohort.data[(cohort.data$event == 1 & cohort.data$include == 1) | 
                          (cohort.data$event == 0 & cohort.data$include == 0),]
  
  #cond.exp2 = matrix(0,no.riskscore.cat,sum(cohort2$event))
  
  samp.size2 = tabulate(cohort2$obs.riskscore.cat.controls)
  
  auc_term2 = sum(sweep(cond.exp2,1,samp.size2,"*"))
  
  cohort3 = cohort.data[(cohort.data$event == 1 & cohort.data$include == 0) | 
                          (cohort.data$event == 0 & cohort.data$include == 1),]
  
  #cond.exp3 = matrix(0,sum(cohort3$event==0),no.riskscore.cat)
  
  samp.size3 = tabulate(cohort3$obs.riskscore.cat.cases)
  
  auc_term3 = sum(sweep(cond.exp3,2,samp.size3,"*"))
  
  #Group 4
  cohort4 = cohort.data[cohort.data$include == 0,]
  
  #E[I(S1>S0)|S1obs,S0obs]
  cond.exp4 = matrix(0,nrow = no.riskscore.cat,ncol = no.riskscore.cat)
  
  cohort4.cases = cohort4[cohort4$event == 1,]
  cohort4.controls = cohort4[cohort4$event == 0,]
  
  for(i in 1:no.riskscore.cat){
    for(j in 1:no.riskscore.cat){
      
      riskscore.cases = cc.cases$riskscore[cc.cases$obs.riskscore.cat.cases == i]
      riskscore.controls = cc.controls$riskscore[cc.controls$obs.riskscore.cat.controls == j]
      
      freq.cases = 1/(cc.cases$sampling.weights[cc.cases$obs.riskscore.cat.cases == i])
      freq.controls = 1/(cc.controls$sampling.weights[cc.controls$obs.riskscore.cat.controls == j])
      
      weight.mat = matrix(kronecker(freq.controls,freq.cases), nrow = length(freq.controls), 
                          byrow = TRUE)
      
      indicator = vapply(riskscore.cases, function(x) x > riskscore.controls,
                         logical(length(riskscore.controls))) * 1
      
      cond.exp4[i,j] = sum(indicator * weight.mat)/sum(weight.mat)
    }
  }
  
  samp.size4 = matrix(kronecker(tabulate(cohort4.cases$obs.riskscore.cat.cases, 
                                         nbins = no.riskscore.cat),
                                tabulate(cohort4.controls$obs.riskscore.cat.controls, 
                                         nbins = no.riskscore.cat)), 
                      nrow = no.riskscore.cat, byrow = TRUE)
  
  auc_term4 = sum(cond.exp4 * samp.size4)
  
  auc_mis = (auc_term1 + auc_term2 + auc_term3 + auc_term4)/(N1 * N0)
  
  case.infl = case.included * (p0_S1 + psi1.by.wt) + 
    (1 - case.included) * (p_S1obs) - auc_mis
  control.infl = control.included * (p1_S0 + psi0.by.wt) + 
    (1 - control.included) * (p_S0obs) - auc_mis
  
  var_auc_mis = mean(case.infl^2)/N1 + mean(control.infl^2)/N0
  
  return(list(Full_Cohort_Estimator = auc_cohort, 
              Variance_Full_Cohort_Estimator = var_auc_cohort,
              IPW_Estimator = auc_ipw,
              Variance_IPW_Estimator = var_auc_ipw,
              Two_Phase_Estimator = auc_mis,
              Variance_Two_Phase_Estimator = var_auc_mis))
  
}


compute_AUC_adj = function(cohort.data, case.control.data, full.model.formula, 
                           reduced.model.formula, risk.factors.logRR, 
                           no.riskscore.cat){
  
  N1 = sum(cohort.data$event)
  N0 = dim(cohort.data)[1] - N1
  
  n1 = sum(case.control.data$event)
  n0 = sum(case.control.data$event == 0)
  
  case.status = cohort.data$event
  
  cohort.covariate.data = model.matrix(reduced.model.formula, 
                                       cohort.data[,all.vars(reduced.model.formula)])[,-1]
  
  cohort.logRR = risk.factors.logRR[colnames(cohort.covariate.data)]
  
  obs.riskscore = as.numeric(cohort.covariate.data %*% cohort.logRR)
  
  obs.riskscore.cat.cases = quantcut(obs.riskscore[case.status == 1],
                                     q = no.riskscore.cat,na.rm = TRUE)
  obs.riskscore.cat.controls = quantcut(obs.riskscore[case.status == 0],
                                        q = no.riskscore.cat,na.rm = TRUE)
  levels(obs.riskscore.cat.cases) = 
    levels(obs.riskscore.cat.controls) = 1:no.riskscore.cat
  
  cohort.data$obs.riskscore = obs.riskscore
  
  cohort.data$obs.riskscore.cat.cases[cohort.data$event == 1] = obs.riskscore.cat.cases
  cohort.data$obs.riskscore.cat.controls[cohort.data$event == 0] = obs.riskscore.cat.controls
  
  case.control.covariate.data = model.matrix(full.model.formula, 
                                             case.control.data[,all.vars(full.model.formula)])[,-1]
  
  case.control.data$riskscore = as.numeric(case.control.covariate.data 
                                           %*% risk.factors.logRR)
  
  case.control.data$obs.riskscore.cat.cases = 
    cohort.data$obs.riskscore.cat.cases[which(is.element(cohort.data$id,case.control.data$id))]
  case.control.data$obs.riskscore.cat.controls = 
    cohort.data$obs.riskscore.cat.controls[which(is.element(cohort.data$id,case.control.data$id))]
  
  cc.cases = case.control.data[case.control.data$event == 1,]
  cc.controls = case.control.data[case.control.data$event == 0,]
  
  cohort.cases = cohort.data[cohort.data$event == 1,]
  cohort.controls = cohort.data[cohort.data$event == 0,]
  
  obs.riskscore.cases = cohort.cases$obs.riskscore
  obs.riskscore.controls = cohort.controls$obs.riskscore
  
  indicator = vapply(obs.riskscore.cases, 
                     function(x) x > obs.riskscore.controls,logical(length(obs.riskscore.controls))) * 1
  
  sampling.weights.phase1.cases = cohort.cases$sampling.weights.phase1
  sampling.weights.phase1.controls = cohort.controls$sampling.weights.phase1
  
  freq.phase1.cases = 1/sampling.weights.phase1.cases
  freq.phase1.controls = 1/sampling.weights.phase1.controls
  
  #compute p_0(S_1) = E_S_0[I(S_1 > S_0)]
  mean_S0_indicator = apply(indicator,2,weighted.mean,w = freq.phase1.controls)
  
  #compute E_S_1[I(S_1 > S_0)]
  mean_S1_indicator = apply(indicator,1,weighted.mean,w = freq.phase1.cases)
  
  weight.mat = matrix(kronecker(freq.phase1.controls,freq.phase1.cases), 
                      nrow = length(freq.phase1.controls), byrow = TRUE)
  
  auc_cohort = sum(indicator * weight.mat)/sum(weight.mat)
  
  var_auc_cohort = (weighted.mean(((mean_S0_indicator - auc_cohort)^2/sampling.weights.phase1.cases), 
                                  w = freq.phase1.cases))/sum(freq.phase1.cases) + 
    (weighted.mean(((mean_S1_indicator - auc_cohort)^2/sampling.weights.phase1.controls), 
                   w = freq.phase1.controls))/sum(freq.phase1.controls)
  
  #Compute the conditional expectation matrix for the categories
  
  #Break up the cohort into 4 groups: (1) case-control pair both included, (2) case included control excluded, (3) case excluded control included, (4) case-control pair excluded
  
  #Computing \psi_1(S1) and \psi_0(S0)
  
  psi1_by_wt = rep(NA,n1)
  psi0_by_wt = rep(NA,n0)
  
  #Compute E[I(S1>S0)|S1,S0obs]
  cond.exp2 = matrix(0,no.riskscore.cat,n1)
  
  #Compute E[I(S1>S0) | S1obs,S0]
  cond.exp3 = matrix(0,n0,no.riskscore.cat)
  
  for(i in 1:no.riskscore.cat){
    
    riskscore.cases = cc.cases$riskscore[cc.cases$obs.riskscore.cat.cases == i]
    riskscore.controls = cc.controls$riskscore
    
    indicator = vapply(riskscore.cases, 
                       function(x) x > riskscore.controls,logical(length(riskscore.controls))) * 1
    
    sampling.weights.cases = cc.cases$sampling.weights[cc.cases$obs.riskscore.cat.cases == i]
    sampling.weights.controls = cc.controls$sampling.weights
    
    freq.cases = 1/sampling.weights.cases
    freq.controls = 1/sampling.weights.controls
    
    weight.mat = matrix(kronecker(freq.controls,freq.cases), 
                        nrow = length(freq.controls), byrow = TRUE)
    
    indicator.times.one_minus_casewt = t(indicator) * (1 - sampling.weights.cases)
    
    first_term = weighted.mean((1 - sampling.weights.cases),w = freq.cases) * 
      apply(indicator,2,weighted.mean,w = freq.controls)
    
    second_term = sum(t(indicator.times.one_minus_casewt) * weight.mat)/sum(weight.mat)
    
    psi1_by_wt[cc.cases$obs.riskscore.cat.cases == i] = 
      (first_term - second_term)/sampling.weights.cases
    
    
    #Compute E[I(S1>S0) | S1obs,S0]
    cond.exp3[,i] = apply(indicator,1,weighted.mean,w = freq.cases)
    
    riskscore.cases = cc.cases$riskscore
    riskscore.controls = cc.controls$riskscore[cc.controls$obs.riskscore.cat.controls == i]
    
    indicator = vapply(riskscore.cases, function(x) x > riskscore.controls,logical(length(riskscore.controls))) * 1
    
    sampling.weights.cases = cc.cases$sampling.weights
    sampling.weights.controls = cc.controls$sampling.weights[cc.controls$obs.riskscore.cat.controls == i]
    
    freq.cases = 1/sampling.weights.cases
    freq.controls = 1/sampling.weights.controls
    
    weight.mat = matrix(kronecker(freq.controls,freq.cases), nrow = length(freq.controls), byrow = TRUE)
    
    indicator.times.one_minus_controlwt = indicator * (1 - sampling.weights.controls)
    
    first_term = weighted.mean((1 - sampling.weights.controls),w = freq.controls) * 
      apply(indicator,1,weighted.mean,w = freq.cases)
    
    second_term = sum(indicator.times.one_minus_controlwt * weight.mat)/sum(weight.mat)
    
    psi0_by_wt[cc.controls$obs.riskscore.cat.controls == i] = (first_term - second_term)/sampling.weights.controls
    
    #Compute E[I(S1>S0)|S1,S0obs]
    cond.exp2[i,] = apply(indicator,2,weighted.mean,w = freq.controls)
    
  }
  
  #Group 1: 
  
  riskscore.cases = cc.cases$riskscore
  riskscore.controls = cc.controls$riskscore
  
  indicator = vapply(riskscore.cases, 
                     function(x) x > riskscore.controls,logical(length(riskscore.controls))) * 1
  
  sampling.weights.phase1.cases = cc.cases$sampling.weights.phase1
  sampling.weights.phase1.controls = cc.controls$sampling.weights.phase1
  
  freq.phase1.cases = 1/sampling.weights.phase1.cases
  freq.phase1.controls = 1/sampling.weights.phase1.controls
  
  weight.mat = matrix(kronecker(freq.phase1.controls,freq.phase1.cases), 
                      nrow = length(freq.phase1.controls), byrow = TRUE)
  
  auc_term1 = sum(indicator * weight.mat)
  
  sampling.weights.cases = cc.cases$sampling.weights
  sampling.weights.controls = cc.controls$sampling.weights
  
  freq.cases = 1/sampling.weights.cases
  freq.controls = 1/sampling.weights.controls
  
  sampling.weights.cases.total = sampling.weights.cases * sampling.weights.phase1.cases
  sampling.weights.controls.total = sampling.weights.controls * sampling.weights.phase1.controls
  
  freq.cases.total = 1/sampling.weights.cases.total
  freq.controls.total = 1/sampling.weights.controls.total
  
  cases.obs.riskscore.cat = cc.cases$obs.riskscore.cat.cases
  controls.obs.riskscore.cat = cc.controls$obs.riskscore.cat.controls
  
  #compute p_0(S_1) = E_S_0[I(S_1 > S_0)]
  mean_S0_indicator = apply(indicator,2,weighted.mean,w = freq.controls)
  
  #compute E_S_1[I(S_1 > S_0)]
  mean_S1_indicator = apply(indicator,1,weighted.mean,w = freq.cases)
  
  ##########[This piece of code also computes the IPW estimator for AUC and its variance]##############
  
  weight.mat = matrix(kronecker(freq.controls.total,freq.cases.total), 
                      nrow = length(freq.controls.total), byrow = TRUE)
  auc_ipw = sum(indicator * weight.mat)/sum(weight.mat)
  
  var_auc_ipw = (weighted.mean(((mean_S0_indicator - auc_ipw)^2/sampling.weights.cases.total), 
                               w = freq.cases.total))/sum(freq.cases.total) + 
    (weighted.mean(((mean_S1_indicator - auc_ipw)^2/sampling.weights.controls.total), 
                   w = freq.controls.total))/sum(freq.controls.total)
  
  
  ####################################################################################################
  
  
  #compute p(S_1_obs) and p(S_0_obs)
  
  mean_S0_indicator.wtcases = mean_S0_indicator * freq.cases
  
  mean_S0_indicator.wtcases.cat = aggregate(mean_S0_indicator.wtcases ~ cases.obs.riskscore.cat, 
                                            FUN = sum)
  wtcases.cat = aggregate(freq.cases ~ cases.obs.riskscore.cat, FUN = sum)
  
  mean_S1S0.S1obs_indicator = mean_S0_indicator.wtcases.cat$mean_S0_indicator.wtcases/wtcases.cat$freq.cases
  
  mean_S1_indicator.wtcontrols = mean_S1_indicator * freq.controls
  
  mean_S1_indicator.wtcontrols.cat = aggregate(mean_S1_indicator.wtcontrols ~ controls.obs.riskscore.cat, FUN = sum)
  wtcontrols.cat = aggregate(freq.controls ~ controls.obs.riskscore.cat, FUN = sum)
  
  mean_S1S0.S0obs_indicator = mean_S1_indicator.wtcontrols.cat$mean_S1_indicator.wtcontrols/wtcontrols.cat$freq.controls
  
  p_S1obs = mean_S1S0.S1obs_indicator[cohort.data$obs.riskscore.cat.cases[cohort.data$event == 1]]
  
  p_S0obs = mean_S1S0.S0obs_indicator[cohort.data$obs.riskscore.cat.controls[cohort.data$event == 0]]
  
  p0_S1 = psi1.by.wt = rep(0,N1)
  p0_S1[which(is.element(cohort.cases$id,cc.cases$id))] = mean_S0_indicator
  
  psi1.by.wt[which(is.element(cohort.cases$id,cc.cases$id))] = psi1_by_wt
  
  p1_S0 = psi0.by.wt = rep(0,N0)
  p1_S0[which(is.element(cohort.controls$id,cc.controls$id))] = mean_S1_indicator
  
  psi0.by.wt[which(is.element(cohort.controls$id,cc.controls$id))] = psi0_by_wt
  
  case.included = cohort.data$include[cohort.data$event == 1]
  control.included = cohort.data$include[cohort.data$event == 0]
  
  cohort.case.weights = cohort.data$sampling.weights[cohort.data$event == 1]
  cohort.control.weights = cohort.data$sampling.weights[cohort.data$event == 0]
  
  
  cohort2 = cohort.data[(cohort.data$event == 1 & cohort.data$include == 1) | 
                          (cohort.data$event == 0 & cohort.data$include == 0),]
  
  
  sampling.weights.phase1.cases = cohort2$sampling.weights.phase1[cohort2$event == 1]
  sampling.weights.phase1.controls = cohort2$sampling.weights.phase1[cohort2$event == 0]
  
  freq.phase1.cases = 1/sampling.weights.phase1.cases
  freq.phase1.controls = 1/sampling.weights.phase1.controls
  
  samp.size2 = as.numeric(wtd.table(x = cohort2$obs.riskscore.cat.controls[cohort2$event == 0], 
                                    weights = freq.phase1.controls, 
                                    type = "table"))
  
  
  auc_term2 = sum(sweep(sweep(cond.exp2,1,samp.size2,"*"),2,freq.phase1.cases,"*"))
  
  cohort3 = cohort.data[(cohort.data$event == 1 & cohort.data$include == 0) | 
                          (cohort.data$event == 0 & cohort.data$include == 1),]
  
  sampling.weights.phase1.cases = cohort3$sampling.weights.phase1[cohort3$event == 1]
  sampling.weights.phase1.controls = cohort3$sampling.weights.phase1[cohort3$event == 0]
  
  freq.phase1.cases = 1/sampling.weights.phase1.cases
  freq.phase1.controls = 1/sampling.weights.phase1.controls
  
  samp.size3 = as.numeric(wtd.table(x = cohort3$obs.riskscore.cat.cases[cohort3$event == 1], 
                                    weights = freq.phase1.cases, 
                                    type = "table"))
  
  
  auc_term3 = sum(sweep(sweep(cond.exp3,2,samp.size3,"*"),1,freq.phase1.controls,"*"))
  
  #Group 4
  cohort4 = cohort.data[cohort.data$include == 0,]
  
  #E[I(S1>S0)|S1obs,S0obs]
  cond.exp4 = matrix(0,nrow = no.riskscore.cat,ncol = no.riskscore.cat)
  
  cohort4.cases = cohort4[cohort4$event == 1,]
  cohort4.controls = cohort4[cohort4$event == 0,]
  
  for(i in 1:no.riskscore.cat){
    for(j in 1:no.riskscore.cat){
      
      riskscore.cases = cc.cases$riskscore[cc.cases$obs.riskscore.cat.cases == i]
      riskscore.controls = cc.controls$riskscore[cc.controls$obs.riskscore.cat.controls == j]
      
      freq.cases = 1/(cc.cases$sampling.weights[cc.cases$obs.riskscore.cat.cases == i])
      freq.controls = 1/(cc.controls$sampling.weights[cc.controls$obs.riskscore.cat.controls == j])
      
      weight.mat = matrix(kronecker(freq.controls,freq.cases), nrow = length(freq.controls), 
                          byrow = TRUE)
      
      indicator = vapply(riskscore.cases, function(x) x > riskscore.controls,
                         logical(length(riskscore.controls))) * 1
      
      cond.exp4[i,j] = sum(indicator * weight.mat)/sum(weight.mat)
    }
  }
  
  
  samp.size4 = matrix(kronecker(as.numeric(wtd.table(x = cohort4$obs.riskscore.cat.cases, 
                                                     weights = 1/cohort4$sampling.weights.phase1,
                                                     type = "table")),
                                as.numeric(wtd.table(x = cohort4$obs.riskscore.cat.controls, 
                                                     weights = 1/cohort4$sampling.weights.phase1,
                                                     type = "table"))), 
                      nrow = no.riskscore.cat, byrow = TRUE)
  
  auc_term4 = sum(cond.exp4 * samp.size4)
  
  sampling.weights.phase1.cases = cohort.cases$sampling.weights.phase1
  sampling.weights.phase1.controls = cohort.controls$sampling.weights.phase1
  
  freq.phase1.cases = 1/sampling.weights.phase1.cases
  freq.phase1.controls = 1/sampling.weights.phase1.controls
  
  auc_mis = (auc_term1 + auc_term2 + auc_term3 + auc_term4)/(sum(freq.phase1.cases) * sum(freq.phase1.controls))
  
  case.infl = case.included * (p0_S1 + psi1.by.wt) + 
    (1 - case.included) * (p_S1obs) - auc_mis
  control.infl = control.included * (p1_S0 + psi0.by.wt) + 
    (1 - control.included) * (p_S0obs) - auc_mis
  
  var_auc_mis = (weighted.mean((case.infl^2/sampling.weights.phase1.cases), 
                               w = freq.phase1.cases))/sum(freq.phase1.cases) + 
    (weighted.mean((control.infl^2/sampling.weights.phase1.controls), 
                   w = freq.phase1.controls))/sum(freq.phase1.controls)
  
  return(list(Adjsuted_Full_Cohort_Estimator = auc_cohort, 
              Variance_Adjsuted_Full_Cohort_Estimator = var_auc_cohort,
              Adjusted_IPW_Estimator = auc_ipw,
              Variance_AdjustedIPW_Estimator = var_auc_ipw,
              Adjusted_Two_Phase_Estimator = auc_mis,
              Variance_Adjusted_Two_Phase_Estimator = var_auc_mis))
  
}

#Illustration to use the compute_AUC function

#Install and load the required libraries if they are not installed
if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}

if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}

if(!require(Hmisc)){
  install.packages("Hmisc")
  library(Hmisc)
}
if(!require(sas7bdat)){
  install.packages("sas7bdat")
  library(sas7bdat)
}


#Load the datasets
load("cohort.data1.txt")
load("case.control.data1.txt")

full.model.formula = ~X1 + X2 + X3 + X4 + factor(X5) + factor(X6) + factor(X7) + factor(X8)
reduced.model.formula =  ~X1 + X2 + X3 + X4 + factor(X5) + factor(X6)

risk.factors.logRR = c(rep(log(1.1),4),rep(log(1.2),2),rep(-log(1.2),2))
names(risk.factors.logRR) = colnames(model.matrix(full.model.formula,
                                                  case.control.data1[,all.vars(full.model.formula)])[,-1])

no.riskscore.cat = 10


output = compute_AUC(cohort.data1, case.control.data1, full.model.formula, 
                     reduced.model.formula, risk.factors.logRR, 
                     no.riskscore.cat = 10)


#Illustration to use the compute_AUC_adj function

#Install and load the required libraries if they are not installed
if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}

if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}

if(!require(Hmisc)){
  install.packages("Hmisc")
  library(Hmisc)
}
if(!require(sas7bdat)){
  install.packages("sas7bdat")
  library(sas7bdat)
}


#Load the datasets
load("cohort.data2.txt")
load("case.control.data2.txt")

full.model.formula = ~X1 + X2 + X3 + X4 + factor(X5) + factor(X6) + factor(X7) + factor(X8)
reduced.model.formula =  ~X1 + X2 + X3 + X4 + factor(X5) + factor(X6)

risk.factors.logRR = c(rep(log(1.1),4),rep(log(1.2),2),rep(-log(1.2),2))
names(risk.factors.logRR) = colnames(model.matrix(full.model.formula,
                                                  case.control.data2[,all.vars(full.model.formula)])[,-1])

no.riskscore.cat = 10


output_adj = compute_AUC_adj(cohort.data2, case.control.data2, full.model.formula, 
                     reduced.model.formula, risk.factors.logRR, 
                     no.riskscore.cat = 10)




