require(foreach);
require(doParallel)
library(tidyverse)
library(MASS)
library(gbm)
library(splitTools)
source("utility.R")

eta = c(1, -0.5, -0.5, -0.5, -0.5, 1)
ps = 5; omegas = 0.5; mus1 = 1; mus0 = -1

p = 6
mu1 = rep(1, p); mu0 = c(rep(0.5, 3), rep(1.5, p-3)); sigma0 = 1/2; sigma1 = 1
alpha = c(1, 0, 1, 0, 1, 0); beta = c(0, 1, 0, 1, 0, 1)
K = 5; thres = 0.05


r_const = 1
N_list = 500 * 2^{c(2:7)}

seed = 2024
cluster_num = 20; rep_num = 1000

r_rate_list = c(1/4, 1/3, 1/2)
r_N_list_all = vector("list", length(r_rate_list))
for (j in 1:length(r_rate_list)){
  r_rate = r_rate_list[j]
  r_N_list_all[[j]] = generate_rN(N_list, rate = r_rate, const = r_const)
}
r_N_list_all[[3]] = 2.5 * r_N_list_all[[3]]



### Plot nuisance estimation errors
j = 1
r_rate = r_rate_list[j]    
res = readRDS(paste0("res_rate", round(r_rate, 2), "_nmax", max(N_list), "_rep", rep_num, ".rds"))

true_mean = calculate_true_mean(omegas, mus1, mus0)

dr_estimates = do.call(rbind,
                       lapply(1:rep_num, function(i) res[[i]]$dr_estimates))
errors_mux1 = do.call(rbind,
                      lapply(1:rep_num, function(i) res[[i]]$errors_mux1_df))
errors_mux0 = do.call(rbind,
                      lapply(1:rep_num, function(i) res[[i]]$errors_mux0_df))
errors_muxs1 = do.call(rbind,
                       lapply(1:rep_num, function(i) res[[i]]$errors_muxs1_df))
errors_muxs0 = do.call(rbind,
                       lapply(1:rep_num, function(i) res[[i]]$errors_muxs0_df))
errors_r = do.call(rbind,
                   lapply(1:rep_num, function(i) res[[i]]$errors_r_df))
errors_e = do.call(rbind,
                   lapply(1:rep_num, function(i) res[[i]]$errors_e_df))
errors_log_lambda = do.call(rbind,
                            lapply(1:rep_num, function(i) res[[i]]$errors_log_lambda_df))

r_plt = as_tibble(errors_r) %>%
  mutate(method = case_when(r_id == 1 ~ "offset-logistic",
                            r_id == 2 ~ "offset-GB")) %>% 
  mutate(N = factor(N)) %>%
  ggplot(aes(x = N, y = error, fill = method)) + ggtitle("Labeling Propensity Score Error") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_boxplot()
lambda_plt = as_tibble(errors_log_lambda) %>%
  mutate(method = case_when(r_id == 1 ~ "offset-logistic",
                            r_id == 2 ~ "offset-GB")) %>%
  mutate(N = factor(N)) %>%
  ggplot(aes(x = N, y = error, fill = method)) + ggtitle("Log Density Ratio Error")  + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_boxplot()
grid_arrange_shared_legend(r_plt, lambda_plt)

### Calculate ATE estimation errors and CI coverage
j = 1
r_rate = r_rate_list[j]
res = readRDS(paste0("res_rate", round(r_rate, 2), "_nmax", max(N_list), "_rep", rep_num, ".rds"))

true_mean = calculate_true_mean(omegas, mus1, mus0)

dr_estimates = as_tibble(dr_estimates) %>%
  mutate(method = case_when(method == 1 ~ "linear",
                            method == 2 ~ "gbm",
                            method == 0 ~ "oracle")) %>%
  mutate(r_N = generate_rN(N, rate = r_rate, const = r_const)) %>%
  mutate(true_mean = true_mean) %>%
  mutate(error = est - true_mean) %>%
  mutate(CI_lower = est - 1.96 * std, CI_upper = est + 1.96 * std) %>% 
  mutate(coverage = (CI_lower <= true_mean) & (true_mean <= CI_upper))

coverage = dr_estimates %>% group_by(method, N) %>%
  summarise(p_coverage = mean(coverage), bias = abs(mean(error)), 
            sd = sd(est), rmse = sqrt(mean(error^2)), length = mean(2 * 1.96 * std),
            se_mean = mean(std))
coverage = as.data.frame(coverage)

