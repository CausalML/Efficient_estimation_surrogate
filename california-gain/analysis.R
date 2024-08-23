source("utility.R")


plot_errors_new <- function(result, nuisance, s_list_plot_error){
  # result = result[sapply(result, length) > 1]
  res_df = convert_res_to_df_obs_ratio(result, obs_ratio_list, nuisance) 
  truth = res_df$estimate[res_df$method == "truth"][1]
  res_df = res_df %>% modify_notations(s_length_list, s_list_plot_error)
  res_df[res_df$method == "athey", "method"] = "SIndex_REG"
  
  temp = res_df %>% 
    filter(!(method %in% c("truth"))) %>% 
    filter(type != "S=36") 
  temp$method[temp$method == "atheyDRPackage"] = "SIndex_DR"
  temp$method[temp$method == "propDR"] = "Proposal"
  
  temp = temp %>% group_by(method, type, obs_ratio) %>%  
    summarise(mean = mean(estimate), std = sd(estimate)) %>% 
    mutate(bias = abs(mean - truth)) %>% 
    mutate(rmse = sqrt(bias^2 + std^2))
  plts = vector("list", length = 3)
  plts[[1]] = temp %>% ggplot(aes(x=type, y=bias, fill=method)) + 
    geom_col(position = "dodge", color = "black") +
    facet_wrap(~obs_ratio) + 
    ylab("Bias")
  plts[[2]] = temp %>% ggplot(aes(x=type, y=std, fill=method)) + 
    geom_col(position = "dodge", color = "black") +
    facet_wrap(~obs_ratio) + ylab("Standard Error")
  plts[[3]] = temp %>% ggplot(aes(x=type, y=rmse, fill=method)) + 
    geom_col(position = "dodge", color = "black") +
    facet_wrap(~obs_ratio) + ylab("RMSE")
  plts
}


###################
#    Riverside
###################
city = "river"
rep_num = 60
core_num = 120
s_length_list = 4*(1:9)
obs_ratio_list = c(0.1, 0.3, 0.5)


s_list_plot_error = c(8, 16, 24, 32)

data = read_csv(paste0("data_", city, "_income.csv"))
result = readRDS(paste0("result_", city, "_rep", rep_num, ".rds"))

nuisance = "grf"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")

nuisance = "xgboost"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")

nuisance = "glmnet"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")


###################
#    San Diego 
###################
city = "sd"
rep_num = 60
core_num = 120
s_length_list = 4*(1:9)
obs_ratio_list = c(0.1, 0.3, 0.5)


s_list_plot_error = c(8, 16, 24, 32)

data = read_csv(paste0("data_", city, "_income.csv"))
result = readRDS(paste0("result_", city, "_rep", rep_num, ".rds"))

nuisance = "grf"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")

nuisance = "xgboost"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")

nuisance = "glmnet"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")


###################
#    Los Angeles  
###################
city = "la"
rep_num = 60
core_num = 120
s_length_list = 4*(1:9)
obs_ratio_list = c(0.1, 0.3, 0.5)


s_list_plot_error = c(8, 16, 24, 32)

data = read_csv(paste0("data_", city, "_income.csv"))
result = readRDS(paste0("result_", city, "_rep", rep_num, ".rds"))

nuisance = "grf"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")

nuisance = "xgboost"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")

nuisance = "glmnet"
plts = plot_errors_new(result, nuisance, s_list_plot_error)
grid_arrange_shared_legend(plts[[1]], plts[[2]], nrow = 2, ncol = 1)
R2 = extract_R2(result, obs_ratio_list, nuisance, data, s_length_list)$R1
R2 %>% ggplot(aes(x = type, y = Rsquare)) + 
  geom_col() + facet_wrap(~obs_ratio) + ylab("Average R Square")

