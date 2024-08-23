source("utility.R")

rep_num = 120
core_num = 60
s_length_list = 4*(1:9)
obs_ratio_list = c(0.1, 0.3, 0.5)

# list of covariates 
x_list = c("age", "xsexf",
           "grde911", "grade12", "grade13above",
           "white", "black", "hisp")  


est_functions_no_S = list(est_truth, est_dr_X_label, est_dr_X_full)
# baseline ATE estimators without using surrogates  
#   est_truth: the oracle difference-in-mean ATE estimator based on the 
#              whole dataset, whose result is viewed as the "truth" for 
#              calculating the bias of other estimators 
#  est_dr_X_label: doubly robust ATE estimator using only labelled data 
#                   ("DR + Labelled" in Figure 1)
#   est_dr_X_full: doubly robust ATE estimator using both labelled and unlabelled 
#                  data, but ignore the surrogates ("DR + s=0" in Figure 1)
est_functions_S = list(est_XS_proposal, est_XS_athey)
# ATE estimators that use the surrogate information 
#   est_XS_proposal: our proposed estimator ("Proposal" in Figure 1)
#   est_XS_athey: the surrogate index estimator based on regression imputation in 
#       Athey et al. [2019] (denoted as "SIndex REG" in Figure 1)
#    the semiparametrically efficient surrogate index estimator proposed in 
#     Chen and Ritzwoller [2023] (denoted as "SIndex DR" in Figure 1) is 
#     implemented by calling the est_athey_dr_package() function below
compute_nuisance_functions = c("grf" = compute_nuisances_rf,
                               "xgboost" = compute_nuisances_gbm,
                               "glmnet" = compute_nuisances_lasso)
# list of estimators used to estimate nuisances: 
#   random forests (grf), gradient boosting (xgboost), lasso (glmnet)
#   here the names grf, xgboost and glmnet are used to keep consistent
#   with the "type" argument of the longterm() function for implementing the 
#   estimator in Chen and Ritzwoller [2023]

experiment_obs_ratio_once <- function(obs_ratio, s_length_list, 
                                      data, x_list, compute_nuisance_functions, 
                                      est_functions_no_S, est_functions_S){
  
  n1 = sum(data$e == 1)
  n0 = sum(data$e == 0)
  ind = partition(data$e, p = c(train = obs_ratio, test = 1-obs_ratio)) 
    # subsample 1-obs_ratio of the data to have missing labels 
  ind_label = ind$train
  ind_miss = ind$test
  data$R = 1
  data$R[ind_miss] = 0
  
  data$RE_type = 4
  data$RE_type[(data$R==0)&(data$e==0)] = 1
  data$RE_type[(data$R==0)&(data$e==1)] = 2
  data$RE_type[(data$R==1)&(data$e==0)] = 3
  folds = create_folds(data$RE_type, k = 5)   
    # split the data into 5 even folds for cross-fitting; the splitting is 
    # stratified on the missingness and the treatment indicators  
  
  est = vector("list", length(compute_nuisance_functions))
  Rsquare = vector("list", length(compute_nuisance_functions))
  for (i in seq_along(compute_nuisance_functions)){
    nuisances = compute_nuisance_functions[[i]](data, folds, x_list, s_length_list, obs_ratio)
      # fit the nuisances using the ith method (i is grf, xgboost, or lasso)
    est_res_no_S = sapply(est_functions_no_S, function(f) f(data, nuisances))
    names(est_res_no_S) = c("truth", "DR_label", "DR_X_full")
      # estimate ATE using estimators that do not need surrogate observations 
    est_res_S = lapply(est_functions_S, function(f) f(data, nuisances))
    est_res_S = unlist(est_res_S)
    est_res_athey_dr = est_athey_dr_package(data, nuisances, s_length_list, x_list, names(compute_nuisance_functions)[i])
      # estimate ATE using estimators that need surrogate observations 
    est[[i]] = c(est_res_no_S, est_res_S, est_res_athey_dr)
    Rsquare[[i]] = compute_S_power(data, nuisances)
      # calculate the cross-validated predictiveness of the surrogates 
  }
  names(est) = names(compute_nuisance_functions)
  names(Rsquare) = names(compute_nuisance_functions)
  list(est = est, Rsquare = Rsquare)
}

set.seed(2024)

city = "river"
data = read_csv(paste0("data_", city, "_income.csv"))
data$Y = data$tcedd36

time1 = proc.time()
cl <- makeCluster(core_num)
registerDoParallel(cl)
result = foreach(i = 1:rep_num, .packages = c("tidyverse", "splitTools", "ranger", "gbm", "glmnet", "np", "mgcv", "longterm")) %dopar% {
  cat(paste("The", i, "th repetition \n"))
  res_temp = vector("list", length = length(obs_ratio_list))
  for (i in seq_along(obs_ratio_list)){
    res_temp[[i]] =  tryCatch({
      experiment_obs_ratio_once(obs_ratio_list[i],
                                s_length_list, data, x_list, compute_nuisance_functions,
                                est_functions_no_S, est_functions_S)
    },  error = function(e) return(paste0("'", e, "'")))
  }
  res_temp
}
stopCluster(cl)
time2 = proc.time()
time2 - time1

saveRDS(result, paste0("result_", city,
                       "_rep", rep_num,
                       ".rds"))

city = "sd"
data = read_csv(paste0("data_", city, "_income.csv"))
data$Y = data$tcedd36

time1 = proc.time()
cl <- makeCluster(core_num)
registerDoParallel(cl)
result = foreach(i = 1:rep_num, .packages = c("tidyverse", "splitTools", "ranger", "gbm", "glmnet", "np", "mgcv", "longterm")) %dopar% {
  cat(paste("The", i, "th repetition \n"))
  res_temp = vector("list", length = length(obs_ratio_list))
  for (i in seq_along(obs_ratio_list)){
    res_temp[[i]] =  tryCatch({
      experiment_obs_ratio_once(obs_ratio_list[i],
                                s_length_list, data, x_list, compute_nuisance_functions,
                                est_functions_no_S, est_functions_S)
    },  error = function(e) return(paste0("'", e, "'")))
  }
  res_temp
}
stopCluster(cl)
time2 = proc.time()
time2 - time1

saveRDS(result, paste0("result_", city,
                       "_rep", rep_num,
                       ".rds"))

city = "la"
data = read_csv(paste0("data_", city, "_income.csv"))
data$Y = data$tcedd36

time1 = proc.time()
cl <- makeCluster(core_num)
registerDoParallel(cl)
result = foreach(i = 1:rep_num, .packages = c("tidyverse", "splitTools", "ranger", "gbm", "glmnet", "np", "mgcv", "longterm")) %dopar% {
  cat(paste("The", i, "th repetition \n"))
  res_temp = vector("list", length = length(obs_ratio_list))
  for (i in seq_along(obs_ratio_list)){
    res_temp[[i]] =  tryCatch({
      experiment_obs_ratio_once(obs_ratio_list[i],
                                s_length_list, data, x_list, compute_nuisance_functions,
                                est_functions_no_S, est_functions_S)
    },  error = function(e) return(paste0("'", e, "'")))
  }
  res_temp
}
stopCluster(cl)
time2 = proc.time()
time2 - time1

saveRDS(result, paste0("result_", city,
                       "_rep", rep_num,
                       ".rds"))
