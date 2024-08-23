library(R.matlab)
library(splitTools)
library(ranger)
library(doParallel)
library(foreach)
library(gbm)
library(glmnet)
library(longterm)
library(lemon)
library(mgcv)
library(np)
library(tidyverse)
select = dplyr::select


#####################
#   ATE Estimation  #
#####################
est_truth <- function(data, nuisances){
  # difference-in-mean ATE estimator using the full data 
  Y = data$Y; e = data$e
  mean(Y[e == 1]) - mean(Y[e==0])
}
est_dr_X_label <- function(data, nuisances){
  # doubly robust ATE estimator using only labelled data ("DR + Labelled" in Figure 1)
  data_temp = data[data$R == 1, ]
  mu_x1 = nuisances$mu_x1[data$R == 1]
  mu_x0 = nuisances$mu_x0[data$R == 1]
  ehat = nuisances$ehat
  Y = data_temp$Y; e = data_temp$e
  
  temp = mu_x1 - mu_x0 +
    (Y - mu_x1) * e/ehat -
    (Y - mu_x0) * (1 - e)/(1 - ehat)
  mean(temp)
}
est_dr_X_full <- function(data, nuisances){
  # doubly robust ATE estimator using both labelled and unlabelled data ("DR + s=0" in Figure 1)
  mu_x1 = nuisances$mu_x1
  mu_x0 = nuisances$mu_x0
  ehat = nuisances$ehat
  rhat = nuisances$rhat
  Y = data$Y; e = data$e; R = data$R
  
  temp = mu_x1 - mu_x0 + 
    (Y - mu_x1) * (e * R)/(ehat * rhat) - 
    (Y - mu_x0) * ((1 - e) * R)/((1 - ehat) * rhat)
  mean(temp)
}
est_XS_proposal <- function(data, nuisances){
  # our proposed estimator ("Proposal" in Figure 1)
  mu_x1 = nuisances$mu_x1
  mu_x0 = nuisances$mu_x0
  ehat = nuisances$ehat
  rhat = nuisances$rhat
  Y = data$Y; e = data$e; R = data$R
  
  estimates = rep(NA, length(nuisances$mu_xs1))
  for (i in seq_along(nuisances$mu_xs1)){
    mu_xs1 = nuisances$mu_xs1[[i]]
    mu_xs0 = nuisances$mu_xs0[[i]]
    temp = mu_x1 - mu_x0 + 
      (Y - mu_xs1) * (e * R)/(ehat *rhat) - 
      (Y - mu_xs0) * ((1 - e) * R)/((1 - ehat) *rhat) +
      (mu_xs1 - mu_x1) * e/ehat -
      (mu_xs0 - mu_x0) * (1 - e)/(1 - ehat)
    estimates[i] = mean(temp)
  }
  names(estimates) = paste0("propDR_", names(nuisances$mu_xs))
  estimates
}
est_XS_athey <- function(data, nuisances){
  # the surrogate index estimator based on regression imputation in 
  # Athey et al. [2019] (denoted as "SIndex REG" in Figure 1)
  
  ehat = nuisances$ehat
  rhat = nuisances$rhat
  Y = data$Y; e = data$e; R = data$R
  estimates = rep(NA, length(nuisances$mu_xs))
  for (i in seq_along(nuisances$mu_xs)){
    mu_xs = nuisances$mu_xs[[i]]
    estimates[i] = mean(mu_xs[e==1]) - mean(mu_xs[e==0])
  }
  names(estimates) = paste0("athey_", names(nuisances$mu_xs))
  estimates
}
est_athey_dr_package <- function(data, nuisances, s_length_list, x_list, type){
  #  the semiparametrically efficient surrogate index estimator proposed in 
  # Chen and Ritzwoller [2023] (denoted as "SIndex DR" in Figure 1)
  
  data$treatment = data$e; data$observe = data$R; 
  X_vars = x_list
  estimates = rep(NA, length(s_length_list))
  for (i in seq_along(s_length_list)){
    S_vars = paste0("tcedd", 1:s_length_list[i])
    temp = longterm(data, S_vars= S_vars, X_vars = X_vars, Y_var = "Y", obs = F, estimand = T, type = type)
    estimates[i] = temp$hat_tau
  }
  names(estimates) = paste0("atheyDRPackage_", names(nuisances$mu_xs))
  estimates
}
########################
#  Nuisance Estimation #
########################
compute_nuisances_rf <- function(data, folds, x_list, s_length_list, obs_ratio){
  ehat = sum(data$e)/nrow(data) 
    # treatment is randomized so the propensity score can be easily estimated 
  rhat = sum(data$R)/nrow(data)   
    # this is a missing-completely-at-random setting so labeling propensity score 
    # can be also easily estimated 
  
  # fit short regressions mu(t, x) 
  mu_x1 = rep(NA, nrow(data))
  mu_x0 = rep(NA, nrow(data))
  ehat_fit = rep(NA, nrow(data))
  for (i in 1:length(folds)){
    fold=folds[[i]]
    data_train = data[fold, ]
    data_eval =  data[-fold, ]
    data_train$e = factor(data_train$e, levels = c(1, 0))
    data_eval$e = factor(data_eval$e, levels = c(1, 0))
    
    form_x = formula(paste("Y ~", reduce(c(x_list), .f = function(x, y) paste(x, "+", y))))
    learner_x1 = ranger(formula = form_x,
                        data = data_train[(data_train$R==1)&(data_train$e==1), ])
    learner_x0 = ranger(formula = form_x,
                        data = data_train[(data_train$R==1)&(data_train$e==0), ])
    mu_x1[-fold] = predict(learner_x1, data = data_eval)$predictions
    mu_x0[-fold] = predict(learner_x0, data = data_eval)$predictions
  }
  
  # fit long regressions: \tilde{mu}(t, x, s)
  mu_xs1 = vector("list", length=length(s_length_list))
  mu_xs0 = vector("list", length=length(s_length_list))
  mu_xs = vector("list", length=length(s_length_list))
  for (j in seq_along(s_length_list)){
    form_xs = formula(paste("Y ~", reduce(c(x_list, paste0("tcedd", 1:s_length_list[j])),
                                          .f = function(x, y) paste(x, "+", y))))
    mu_xs1[[j]] = rep(NA, nrow(data))
    mu_xs0[[j]] = rep(NA, nrow(data))
    mu_xs[[j]] = rep(NA, nrow(data))
    for (i in 1:length(folds)){
      fold=folds[[i]]
      data_train = data[fold, ]
      data_eval =  data[-fold, ]
      
      learner_xs1 = ranger(formula = form_xs,
                           data = data_train[(data_train$R==1)&(data_train$e==1), ])
      learner_xs0 = ranger(formula = form_xs,
                           data = data_train[(data_train$R==1)&(data_train$e==0), ])
      learner_xs = ranger(formula = form_xs,
                          data = data_train[(data_train$R==1), ])
      mu_xs1[[j]][-fold] = predict(learner_xs1, data = data_eval)$predictions
      mu_xs0[[j]][-fold] = predict(learner_xs0, data = data_eval)$predictions
      mu_xs[[j]][-fold] = predict(learner_xs, data = data_eval)$predictions
    }
  }
  names(mu_xs1) = paste0("s", s_length_list)
  names(mu_xs0) = paste0("s", s_length_list)
  names(mu_xs) = paste0("s", s_length_list)
  
  list(mu_xs1 = mu_xs1, mu_xs0 = mu_xs0, mu_xs = mu_xs,
       mu_x1 = mu_x1, mu_x0 = mu_x0,
       ehat = ehat, rhat = rhat)
}
compute_nuisances_gbm <- function(data, folds, x_list, s_length_list, obs_ratio){
  ehat = sum(data$e)/nrow(data) 
  # treatment is randomized so the propensity score can be easily estimated 
  rhat = sum(data$R)/nrow(data)   
  # this is a missing-completely-at-random setting so labeling propensity score 
  # can be also easily estimated 
  
  # fit short regressions mu(t, x)
  mu_x1 = rep(NA, nrow(data))
  mu_x0 = rep(NA, nrow(data))
  ehat_fit = rep(NA, nrow(data))
  for (i in 1:length(folds)){
    fold=folds[[i]]
    data_train = data[fold, ]
    data_eval =  data[-fold, ]
    
    form_x = formula(paste("Y ~", reduce(c(x_list), .f = function(x, y) paste(x, "+", y))))
    learner_x1 = gbm(formula = form_x, distribution = "gaussian",
                     n.trees = 500,
                     data = data_train[(data_train$R==1)&(data_train$e==1), ], n.minobsinnode = 5)
    learner_x0 = gbm(formula = form_x, distribution = "gaussian",
                     n.trees = 500,
                     data = data_train[(data_train$R==1)&(data_train$e==0), ], n.minobsinnode = 5)
    mu_x1[-fold] = predict(learner_x1, data_eval)
    mu_x0[-fold] = predict(learner_x0, data_eval)
  }
  
  # fit long regressions: \tilde{mu}(t, x, s)
  mu_xs1 = vector("list", length=length(s_length_list))
  mu_xs0 = vector("list", length=length(s_length_list))
  mu_xs = vector("list", length=length(s_length_list))
  for (j in seq_along(s_length_list)){
    form_xs = formula(paste("Y ~", reduce(c(x_list, paste0("tcedd", 1:s_length_list[j])),
                                          .f = function(x, y) paste(x, "+", y))))
    mu_xs1[[j]] = rep(NA, nrow(data))
    mu_xs0[[j]] = rep(NA, nrow(data))
    mu_xs[[j]] = rep(NA, nrow(data))
    for (i in 1:length(folds)){
      fold=folds[[i]]
      data_train = data[fold, ]
      data_eval =  data[-fold, ]
      
      learner_xs1 = gbm(formula = form_xs,  distribution = "gaussian",
                        n.trees = 500,
                        data = data_train[(data_train$R==1)&(data_train$e==1), ], n.minobsinnode = 5)
      learner_xs0 = gbm(formula = form_xs,  distribution = "gaussian",
                        n.trees = 500,
                        data = data_train[(data_train$R==1)&(data_train$e==0), ], n.minobsinnode = 5)
      learner_xs = gbm(formula = form_xs,  distribution = "gaussian",
                       n.trees = 500,
                       data = data_train[(data_train$R==1), ], n.minobsinnode = 5)
      mu_xs1[[j]][-fold] = predict(learner_xs1, data_eval)
      mu_xs0[[j]][-fold] = predict(learner_xs0, data_eval)
      mu_xs[[j]][-fold] = predict(learner_xs, data_eval)
    }
  }
  names(mu_xs1) = paste0("s", s_length_list)
  names(mu_xs0) = paste0("s", s_length_list)
  names(mu_xs) = paste0("s", s_length_list)
  
  list(mu_xs1 = mu_xs1, mu_xs0 = mu_xs0, mu_xs = mu_xs,
       mu_x1 = mu_x1, mu_x0 = mu_x0,
       ehat = ehat, rhat = rhat)
}
compute_nuisances_lasso <- function(data, folds, x_list, s_length_list, obs_ratio){
  ehat = sum(data$e)/nrow(data) 
  # treatment is randomized so the propensity score can be easily estimated 
  rhat = sum(data$R)/nrow(data)   
  # this is a missing-completely-at-random setting so labeling propensity score 
  # can be also easily estimated 
  
  # fit short regressions mu(t, x)
  mu_x1 = rep(NA, nrow(data))
  mu_x0 = rep(NA, nrow(data))
  ehat_fit = rep(NA, nrow(data))
  for (i in 1:length(folds)){
    fold=folds[[i]]
    data_train = data[fold, ]
    data_eval =  data[-fold, ]
    
    learner_x1 = cv.glmnet(as.matrix(data_train[(data_train$R==1)&(data_train$e==1), x_list]),
                           data_train$Y[(data_train$R==1)&(data_train$e==1)], family = "gaussian", alpha = 1)
    learner_x0 = cv.glmnet(as.matrix(data_train[(data_train$R==1)&(data_train$e==0), x_list]),
                           data_train$Y[(data_train$R==1)&(data_train$e==0)], family = "gaussian", alpha = 1)
    mu_x1[-fold] = c(predict(learner_x1, newx = as.matrix(data_eval[, x_list]), s = "lambda.min"))
    mu_x0[-fold] = c(predict(learner_x0, newx = as.matrix(data_eval[, x_list]), s = "lambda.min"))
  }
  
  # fit long regressions: \tilde{mu}(t, x, s)
  mu_xs1 = vector("list", length=length(s_length_list))
  mu_xs0 = vector("list", length=length(s_length_list))
  mu_xs = vector("list", length=length(s_length_list))
  for (j in seq_along(s_length_list)){
    xs_list = c(x_list, paste0("tcedd", 1:s_length_list[j]))
    
    mu_xs1[[j]] = rep(NA, nrow(data))
    mu_xs0[[j]] = rep(NA, nrow(data))
    mu_xs[[j]] = rep(NA, nrow(data))
    for (i in 1:length(folds)){
      fold=folds[[i]]
      data_train = data[fold, ]
      data_eval =  data[-fold, ]
      
      learner_xs1 = cv.glmnet(as.matrix(data_train[(data_train$R==1)&(data_train$e==1), xs_list]),
                              data_train$Y[(data_train$R==1)&(data_train$e==1)], family = "gaussian", alpha = 1)
      learner_xs0 = cv.glmnet(as.matrix(data_train[(data_train$R==1)&(data_train$e==0), xs_list]),
                              data_train$Y[(data_train$R==1)&(data_train$e==0)], family = "gaussian", alpha = 1)
      learner_xs = cv.glmnet(as.matrix(data_train[(data_train$R==1), xs_list]),
                             data_train$Y[(data_train$R==1)], family = "gaussian", alpha = 1)
      mu_xs1[[j]][-fold] = c(predict(learner_xs1, newx = as.matrix(data_eval[, xs_list]), s = "lambda.min"))
      mu_xs0[[j]][-fold] = c(predict(learner_xs0, newx = as.matrix(data_eval[, xs_list]), s = "lambda.min"))
      mu_xs[[j]][-fold] = c(predict(learner_xs0, newx = as.matrix(data_eval[, xs_list]), s = "lambda.min"))
    }
  }
  names(mu_xs1) = paste0("s", s_length_list)
  names(mu_xs0) = paste0("s", s_length_list)
  names(mu_xs) = paste0("s", s_length_list)
  
  list(mu_xs1 = mu_xs1, mu_xs0 = mu_xs0, mu_xs = mu_xs,
       mu_x1 = mu_x1, mu_x0 = mu_x0,
       ehat = ehat, rhat = rhat)
}
########################
#      Evaluation      #
########################
compute_S_power <- function(data, nuisances){
  # calculate the cross-validated predictiveness of surrogates according to our 
  # nuisance estimates of \tilde{mu}(t, x, s) and mu(t, x)
  total1 = mean((data$Y[data$e == 1] - nuisances$mu_x1[data$e == 1])^2)
  total0 = mean((data$Y[data$e == 0] - nuisances$mu_x0[data$e == 0])^2)
  sub1 = rep(NA, length(nuisances$mu_xs1)); sub0 = rep(NA, length(nuisances$mu_xs0))
  for (i in 1:length(nuisances$mu_xs1)){
    sub1[i] = mean((data$Y[data$e == 1] - nuisances$mu_xs1[[i]][data$e == 1])^2)
    sub0[i] = mean((data$Y[data$e == 0] - nuisances$mu_xs0[[i]][data$e == 0])^2)
  }
  R1 = 1 - sub1/total1
  R0 = 1 - sub0/total0
  list(R1 = R1, R0 = R0)
}
convert_res_to_df_obs_ratio <- function(result, obs_ratio_list, nuisance){
  # convert the .rds results to a data frame
  temp = vector("list", length=length(obs_ratio_list))
  for (i in seq_along(obs_ratio_list)){
    temp_inner = transpose(result)[[i]]
    temp[[i]] = as.data.frame(do.call(rbind, lapply(temp_inner, function(x) x[["est"]][[nuisance]])))
    temp[[i]] = temp[[i]] %>% pivot_longer(everything(), names_to = "method", 
                                           values_to = "estimate") %>% 
      mutate(obs_ratio = paste0("r = ", obs_ratio_list[i]))
  }
  do.call(rbind, temp)
}
extract_R2 <- function(result, obs_ratio_list, nuisance, data, s_length_list){
  # extract the R square from the results 
  pi = mean(data$e)
  
  R1_temp = vector("list", length=length(obs_ratio_list))
  R0_temp = vector("list", length=length(obs_ratio_list))
  R_temp = vector("list", length=length(obs_ratio_list))
  for (i in seq_along(obs_ratio_list)){
    temp_inner = transpose(result)[[i]]
    R1_temp[[i]] = as.data.frame(do.call(rbind, lapply(temp_inner, function(x) x[["Rsquare"]][[nuisance]][["R1"]]))) %>% 
      colMeans()
    R0_temp[[i]] = as.data.frame(do.call(rbind, lapply(temp_inner, function(x) x[["Rsquare"]][[nuisance]][["R0"]]))) %>% 
      colMeans()
    R_temp[[i]] = pi * R1_temp[[i]] + (1 - pi) * R0_temp[[i]]
  }
  
  list_to_df <- function(temp){
    temp = as.data.frame(do.call(rbind, temp))
    temp$obs_ratio = paste0("r = ", obs_ratio_list)
    colnames(temp)[-ncol(temp)] = paste0("s=", 4*(1:ncol(temp[, -1])))
    temp = temp %>% pivot_longer(starts_with("s="), names_to = "type", values_to = "Rsquare") %>% filter(type != "S=36")
    temp$type = factor(temp$type, levels =  paste0("s=", s_length_list))
    temp
  }
  
  list(R1 = list_to_df(R1_temp), R0 = list_to_df(R0_temp), R = list_to_df(R_temp))
}
modify_notations <- function(res_df, s_length_list, s_list_plot_error){
  # modify a few notations 
  res_df$S = 0
  ind = grepl("_s", res_df$method) 
  res_df$S[ind] = sapply(str_split(res_df$method[ind], "_s"), function(x) as.numeric(x[2]))
  res_df$type = NA
  res_df$type[grepl("label", res_df$method)] = "Labelled"
  res_df$type[grepl("full", res_df$method)] = "s=0"
  res_df$type[grepl("truth", res_df$method)] = "truth"
  res_df$type[is.na(res_df$type)] = paste0("s=", res_df$S[is.na(res_df$type)])
  res_df = res_df %>% filter(type %in% c("Labelled", "s=0", paste0("s=", s_list_plot_error)))
  res_df$type = factor(res_df$type, levels = c("Labelled", "s=0", paste0("s=", s_list_plot_error)))
  res_df$method = sapply(str_split(res_df$method, "_"), function(x) x[1])
  res_df
}
