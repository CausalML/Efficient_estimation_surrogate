require(foreach);
require(doParallel)
library(tidyverse)
library(MASS)
library(gbm)
library(splitTools)
source("utility.R")

experiment_once <- function(N_list, r_N_list, p, K,
                            mu1, mu0, sigma0, sigma1, ps, omegas, mus1, mus0, eta, alpha, beta, thres,
                            est_r_functions, extract_r_functions, est_e_functions, 
                            est_mux_functions, est_muxs_functions){
  
  errors_r = matrix(0, length(N_list), length(est_r_functions))
  errors_log_lambda = matrix(0, length(N_list), length(est_r_functions))
  errors_e = matrix(0, length(N_list), length(est_e_functions))
  errors_mux1 = matrix(0, length(N_list), length(est_mux_functions))
  errors_mux0 = matrix(0, length(N_list), length(est_mux_functions))
  errors_muxs1 = matrix(0, length(N_list), length(est_muxs_functions))
  errors_muxs0 = matrix(0, length(N_list), length(est_muxs_functions))
  
  dr_estimate_list = vector("list", length = length(N_list))
  
  for (i in 1:length(N_list)){
    N = N_list[i]; r_N = r_N_list[i]
    df_XRT = simuate_T(simulate_XR(N, p, r_N, mu1, mu0, sigma1, sigma0), eta, thres = thres)
    df = simulate_SY(df_XRT, ps, mus1, mus0, omegas, beta, alpha)
    
    folds = create_folds(df$R, k = K, type = "stratified")
    
    # calculate values of true nuisances 
    gamma_true = calculate_true_coef(mu1, mu0, sigma1, sigma0)
    true_muxs = calculate_true_muxs(df, p, ps, omegas, mus1, mus0, beta, alpha)
    true_mux = calculate_true_mux(df, p, ps, omegas, mus1, mus0, beta, alpha)
    true_r_lambda = calculate_true_r_lambda(df, gamma_true, df$r_N[1])
    df$muxs1 = true_muxs$muxs1; df$muxs0 = true_muxs$muxs0
    df$mux1 = true_mux$mux1; df$mux0 = true_mux$mux0
    df$r=true_r_lambda$r; df$log_lambda=true_r_lambda$log_lambda
    
    # matrices to store estimated values of nuisances 
    r_hat = matrix(NA, N, length(est_r_functions))
    log_lambda_hat = matrix(NA, N, length(est_r_functions))
    muxs1_hat = matrix(NA, N, length(est_muxs_functions)); muxs0_hat = matrix(NA, N, length(est_muxs_functions))
    mux1_hat = matrix(NA, N, length(est_mux_functions)); mux0_hat = matrix(NA, N, length(est_mux_functions))
    e_hat = matrix(NA, N, length(est_e_functions))
    
    # cross-fitting for nuisance estimation 
    for (j in 1:K){
      ind = folds[[j]]
      train = df[ind, ]
      valid = df[-ind, ]
      
      # fitting r 
      for (k in 1:length(est_r_functions)){
        fit = est_r_functions[[k]](train, p)
        temp = extract_r_functions[[k]](fit, valid)
        r_hat[-ind, k] = temp$pi_hat; log_lambda_hat[-ind, k] = temp$log_lambda_hat
      }
      
      # fitting e
      for (k in 1:length(est_e_functions)){
        fit = est_e_functions[[k]](train, p)
        e_hat[-ind, k] = truncate(predict(fit,newdata=valid,type="response"), thres = thres)
      }
      
      # fitting muxs
      for (k in 1:length(est_muxs_functions)){
        fit1 = est_muxs_functions[[k]](train[train$Treat==1, ], p, ps)
        fit0 = est_muxs_functions[[k]](train[train$Treat==0, ], p, ps)
        
        muxs1_hat[-ind, k] = predict(fit1,newdata=valid) 
        muxs0_hat[-ind, k] = predict(fit0,newdata=valid) 
      }
      
      # fitting mux
      for (k in 1:length(est_mux_functions)){
        fit1 = est_mux_functions[[k]](train[train$Treat==1, ], p, ps)
        fit0 = est_mux_functions[[k]](train[train$Treat==0, ], p, ps)
        
        mux1_hat[-ind, k] = predict(fit1,newdata=valid) 
        mux0_hat[-ind, k] = predict(fit0,newdata=valid)
      }
    }
    
    # calucate the errors of nuisance estimators 
    errors_r[i, ] = apply(1-r_hat/c(df$r), 2, mse)
    errors_log_lambda[i, ] = apply(log_lambda_hat - c(df$log_lambda), 2, mse)
    errors_e[i, ] = apply(e_hat - c(df$e), 2, mse)
    
    errors_mux1[i, ] = apply(mux1_hat - c(df$mux1), 2, mse)
    errors_mux0[i, ] = apply(mux0_hat - c(df$mux0), 2, mse)
    
    errors_muxs1[i, ] = apply(muxs1_hat - c(df$muxs1), 2, mse)
    errors_muxs0[i, ] = apply(muxs0_hat - c(df$muxs0), 2, mse)
    
    # calculate ATE estimators based on three different types of nuisance estimators
    #     method = 1: nuisances are estimated by linear models 
    #     method = 2: nuisances are estimated by GBM
    #     method = 0: nuisances true values are used 
    dr_estimates = matrix(0, length(est_r_functions)+1, 4)
    colnames(dr_estimates) = c("method", "N", "est", "std")
    dr_estimates[, "N"] = N
    for (j in seq_along(est_r_functions)){
      r_hat0 = r_hat[, j]; e_hat0 = e_hat[, j]
      muxs1_hat0 = muxs1_hat[, j]; muxs0_hat0 = muxs0_hat[, j]
      mux1_hat0 = mux1_hat[, j]; mux0_hat0 = mux0_hat[, j]
      
      est_sn = estimate_dr_sn(r_hat0, e_hat0, muxs1_hat0, muxs0_hat0, mux1_hat0, mux0_hat0, df) 
      dr_estimates[j, "method"] = j
      # method = 1 indicates linear, method = 2 indicates GBM 
      dr_estimates[j, "est"] = est_sn;   # point estimation 
      dr_estimates[j, "std"] = estimate_dr_std_limit_sn(r_hat0, e_hat0, muxs1_hat0, muxs0_hat0, mux1_hat0, mux0_hat0, df, est_sn)
      # standard error estimation 
    }
    muxs1_hat0 = c(true_muxs$muxs1); muxs0_hat0 = c(true_muxs$muxs0); 
    mux1_hat0 = c(true_mux$mux1); mux0_hat0 = c(true_mux$mux0)
    r_hat0=c(true_r_lambda$r); e_hat0 = c(df$e)
    dr_estimates[length(est_r_functions)+1, "method"] = 0  # method = 0 indicates oracle 
    est_sn = estimate_dr_sn(r_hat0, e_hat0, muxs1_hat0, muxs0_hat0, mux1_hat0, mux0_hat0, df)
    dr_estimates[length(est_r_functions)+1, "est"] = est_sn;    # point estimation 
    dr_estimates[length(est_r_functions)+1, "std"] = estimate_dr_std_limit_sn(r_hat0, e_hat0, muxs1_hat0, muxs0_hat0, mux1_hat0, mux0_hat0, df, est_sn)
    # standard error estimation 
    
    dr_estimate_list[[i]] = dr_estimates
  }
  
  dr_estimates_df = do.call(rbind, dr_estimate_list)
  errors_mux1_df = data.frame(error = c(errors_mux1), N = rep(N_list, ncol(errors_mux1)), 
                              mu_id = rep(1:ncol(errors_mux1), each = length(N_list)))
  errors_mux0_df = data.frame(error = c(errors_mux0), N = rep(N_list, ncol(errors_mux0)), 
                              mu_id = rep(1:ncol(errors_mux0), each = length(N_list)))
  errors_muxs1_df = data.frame(error = c(errors_muxs1), N = rep(N_list, ncol(errors_muxs1)), 
                               mu_id = rep(1:ncol(errors_muxs1), each = length(N_list)))
  errors_muxs0_df = data.frame(error = c(errors_muxs0), N = rep(N_list, ncol(errors_muxs0)), 
                               mu_id = rep(1:ncol(errors_muxs0), each = length(N_list)))
  errors_r_df = data.frame(error = c(errors_r), N = rep(N_list, ncol(errors_r)), 
                           r_id = rep(1:ncol(errors_r), each = length(N_list)))
  errors_e_df = data.frame(error = c(errors_e), N = rep(N_list, ncol(errors_e)), 
                           e_id = rep(1:ncol(errors_e), each = length(N_list)))
  errors_log_lambda_df = data.frame(error = c(errors_log_lambda), N = rep(N_list, ncol(errors_log_lambda)), 
                                    r_id = rep(1:ncol(errors_log_lambda), each = length(N_list)))
  
  list(dr_estimates_df = dr_estimates_df,
       errors_mux1_df = errors_mux1_df, errors_mux0_df = errors_mux0_df,
       errors_muxs1_df = errors_muxs1_df, errors_muxs0_df = errors_muxs0_df,
       errors_r_df = errors_r_df, errors_e_df = errors_e_df, errors_log_lambda_df = errors_log_lambda_df)
}


est_r_functions = c(est_r_logistic_square_offset, est_r_gbm_offset)  
extract_r_functions = c(extract_r_logistic_offset_fit, extract_r_gbm_offset_fit)
est_e_functions = c(est_e_logistic, est_e_gbm)
est_mux_functions = c(est_mux_square, est_mux_gbm)
est_muxs_functions = c(est_muxs_square, est_muxs_gbm)



eta = c(1, -0.5, -0.5, -0.5, -0.5, 1)
ps = 5; omegas = 0.5; mus1 = 1; mus0 = -1

p = 6
mu1 = rep(1, p); mu0 = c(rep(0.5, 3), rep(1.5, p-3)); sigma0 = 1/2; sigma1 = 1
alpha = c(1, 0, 1, 0, 1, 0); beta = c(0, 1, 0, 1, 0, 1)
K = 5; thres = 0.05


r_const = 1
N_list = 500 * 2^{c(2:7)}

seed = 2024
cluster_num = 10; rep_num = 1000   # use 10 cores to run 1000 replications of experiments 

r_rate_list = c(1/4, 1/3, 1/2)   
  # size of labelled data = N^{-r_rate} for r_rate in {1/4, 1/3, 1/2}
r_N_list_all = vector("list", length(r_rate_list))
for (j in 1:length(r_rate_list)){
  r_rate = r_rate_list[j]
  r_N_list_all[[j]] = generate_rN(N_list, rate = r_rate, const = r_const)
}
r_N_list_all[[3]] = 2.5 * r_N_list_all[[3]]

for (j in c(1:3)){
  r_rate = r_rate_list[j]
  r_N_list = r_N_list_all[[j]]
  
  cat("The ", j, "th iteration", "\n")
  
  time0 = Sys.time()
  set.seed(seed)
  cluster = makeCluster(cluster_num)
  registerDoParallel(cluster)
  res = foreach(i = 1:rep_num, .packages = c("MASS", "gbm", "splitTools")) %dopar% {
    tryCatch({
      experiment_once(N_list, r_N_list, p, K,
                      mu1, mu0, sigma0, sigma1, ps, omegas, mus1, mus0, eta, alpha, beta, thres,
                      est_r_functions, extract_r_functions, est_e_functions, 
                      est_mux_functions, est_muxs_functions)
    },  error = function(e) return(paste0("'", e, "'")))
  }
  stopCluster(cluster)
  time1 = Sys.time()
  print(time1 - time0)
  
  saveRDS(res, paste0("res_rate", round(r_rate, 2), "_nmax", max(N_list), "_rep", rep_num, ".rds"))
}



