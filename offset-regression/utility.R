### Data Simulation 
simulate_XR <- function(N, p, r_N, mu1, mu0, sigma1, sigma0){
  R = rbinom(N,1,r_N)
  X = matrix(0, N, p)
  X[R == 1, ] = mvrnorm(sum(R), mu1, sigma1 * diag(p))
  X[R == 0, ] = mvrnorm(N - sum(R), mu0, sigma0 * diag(p))
  data.frame(X, R, r_N = mean(R))
}
truncate <- function(e, thres){
  e[e<=thres] = thres
  e[e>=1-thres] = 1-thres
  e
}
simuate_T <- function(df_XR, eta, thres = 0.05){
  p = length(eta)
  Xname = paste0("X", 1:p)
  X = as.matrix(df_XR[, Xname])
  R = df_XR$R
  
  e1 = 1/(1 + exp(-c(X[R==1, ] %*% eta)))
  e1 = truncate(e1, thres)
  e0 = 1/(1 + exp(-c(X[R==0, ] %*% eta)))
  e0 = truncate(e0, thres)
  
  df_XR$e[R==1] = e1; df_XR$e[R==0] = e0
  df_XR$Treat = rbinom(nrow(X), 1, df_XR$e)
  
  df_XR
}
simulate_S <- function(df_XRT, ps, mus1, mus0){
  N = nrow(df_XRT)
  treat = df_XRT$Treat
  S0 = matrix(rnorm(N*ps, mean = mus0, sd = 1), N, ps)
  S1 = matrix(rnorm(N*ps, mean = mus1, sd = 1), N, ps)
  S = matrix(0, N, ps)
  S[treat==1, ] = S1[treat==1, ]; S[treat==0, ] = S0[treat==0, ]
  colnames(S) = paste0("S", 1:ps); colnames(S0) = paste0("S0_", 1:ps); colnames(S1) = paste0("S1_", 1:ps)
  
  list(S = S, S0 = S0, S1 = S1)
}
simulate_SY <- function(df_XRT, ps, mus1, mus0, omegas, beta, alpha){
  px = length(beta)
  N = nrow(df_XRT)
  
  Xname = paste0("X", 1:px)
  X = df_XRT[, Xname]
  treat = df_XRT$Treat
  Y = rep(0, N)
  
  
  Slist = simulate_S(df_XRT, ps, mus1, mus0)
  
  Y1 = 1 + omegas * rowMeans(Slist$S1) + as.matrix(X) %*% beta + 
    as.matrix(X)^2 %*% alpha + rnorm(nrow(df_XRT))
  Y0 = -1 - omegas * rowMeans(Slist$S0) + as.matrix(X) %*% beta + 
    as.matrix(X)^2 %*% alpha + rnorm(nrow(df_XRT))
  Y[treat == 1] = Y1[treat==1]; Y[treat == 0] = Y0[treat==0]
  df_XRT$Y1 = Y1; df_XRT$Y0 = Y0; df_XRT$Y = Y
  
  cbind(df_XRT, Slist$S, Slist$S1, Slist$S0)
}



### Calculate true values of ATE and the nuisance functions
calculate_true_mean <- function(omegas, mus1, mus0){
  2 + omegas*mus1 + omegas*mus0
}
generate_rN <- function(N, rate = 1/3, const = 1/2){
  # generate the marginal labeling proportion r_N (refered to as pi_N in the paper)
  const * N^{-rate}
}
calculate_true_coef <- function(mu1, mu0, sigma1, sigma0){
  # calculate the true coefficients of the offset logistic regression 
  # corresponding to our data generating process
  p = length(mu0)
  gamma2 = 1/2 * rep(1/sigma0 - 1/sigma1, p)
  gamma1 = -mu1*(1/sigma0 - 1/sigma1) + (mu1 - mu0) * (1/sigma0)
  gamma0 = 1/2 * (sum(mu0^2/sigma0) - sum(mu1^2/sigma1)) + 1/2 * p * log(sigma0/sigma1)
  list(gamma2 = gamma2, gamma1 = gamma1, gamma0 = gamma0)
}
calculate_true_r_lambda <- function(df, gamma, r_N){
  # calculate the true lambda*_N and r*_N induced by the observation model 
  # given the gamma coefficients calculated in calculate_true_coef() and 
  # the marginal proportion r_N 
  
  X = df[, paste0("X", 1:p)]
  
  log_lambda = as.matrix(X^2)%*%gamma$gamma2 + as.matrix(X)%*%gamma$gamma1 + gamma$gamma0
  
  r = (r_N/(1-r_N))*exp(log_lambda)/(1+(r_N/(1-r_N))*exp(log_lambda))
  list(r = r, log_lambda = log_lambda)
}
calculate_true_muxs <- function(df, p, ps, omegas, mus1, mus0, beta, alpha){
  X = df[, paste0("X", 1:p)]; S = df[, paste0("S", 1:ps)] 
  
  muxs1 = 1 + omegas * rowMeans(S) + as.matrix(X) %*% beta + as.matrix(X)^2 %*% alpha 
  muxs0 = -1 - omegas * rowMeans(S) + as.matrix(X) %*% beta + as.matrix(X)^2 %*% alpha 
  
  list(muxs1 = muxs1, muxs0 = muxs0)
}
calculate_true_mux <- function(df, p, ps, omegas, mus1, mus0, beta, alpha){
  X = df[, paste0("X", 1:p)]
  
  mux1 = 1 + omegas * mus1 + as.matrix(X) %*% beta + as.matrix(X)^2 %*% alpha 
  mux0 = -1 - omegas * mus0 + as.matrix(X) %*% beta + as.matrix(X)^2 %*% alpha  
  list(mux1 = mux1, mux0 = mux0)
}


### Parametric Model Nuisance Estimation
est_r_logistic_square_offset <- function(df, p){
  Xname = paste0("X", 1:p)
  form1 = paste(Xname, collapse = " + ")
  form2 = paste0("I(", Xname, "^2)", collapse = " + ")
  form = formula(paste("R ~", paste(c(form1, form2, "offset(log(r_N/(1-r_N)))"), collapse = " + ")))
  fit = glm(form, data = df, family=binomial)
  fit
}
est_e_logistic <- function(df, p){
  Xname = paste0("X", 1:p)
  form = formula(paste("Treat ~ 0 + ", paste(Xname, collapse = " + "), collapse = " + "))
  fit = glm(form, data = df, family=binomial)
  fit
}
est_muxs_square <- function(df, p, ps){
  Xname = paste0("X", 1:p); Sname = paste0("S", 1:ps)
  form1 = paste(Xname, collapse = " + ")
  form2 = paste0("I(", Xname, "^2)", collapse = " + ")
  form3 = paste(Sname, collapse = " + ")
  form = formula(paste("Y ~", paste(c(form1, form2, form3), collapse = "+")))
  fit = lm(form, data = df, subset = (R == 1))
  fit
}
est_mux_square <- function(df, p, ps){
  Xname = paste0("X", 1:p); 
  form1 = paste(Xname, collapse = " + ")
  form2 = paste0("I(", Xname, "^2)", collapse = " + ")
  form = formula(paste("Y ~", paste(c(form1, form2), collapse = "+")))
  fit = lm(form, data = df, subset = (R == 1))
  fit
}


### GBM Nuisance Estimation
est_r_gbm_offset <- function(df, p){
  Xname = paste0("X", 1:p)
  form1 = paste(Xname, collapse = " + ")
  form = formula(paste("R ~", paste(c(form1, "offset(log(r_N/(1-r_N)))"), collapse = " + ")))
  fit = gbm(formula = form, data = df, distribution = "bernoulli", n.trees = 1000, shrinkage = 0.05, bag.fraction = 0.9, n.minobsinnode = 5)
  fit 
}
est_e_gbm <- function(df, p){
  Xname = paste0("X", 1:p)
  form = formula(paste("Treat ~ ", paste(Xname, collapse = " + "), collapse = " + "))
  fit = gbm(formula = form, data = df, distribution = "bernoulli", n.trees = 1000, shrinkage = 0.05, bag.fraction = 0.9, n.minobsinnode = 5)
  fit
}
est_muxs_gbm <- function(df, p, ps){
  Xname = paste0("X", 1:p); Sname = paste0("S", 1:ps)
  form = formula(paste("Y ~ ", paste(c(Xname, Sname), collapse = " + "), collapse = " + "))
  fit = gbm(form, data = df[df$R == 1, ],
            distribution = "gaussian", n.trees = 1000, shrinkage = 0.05, bag.fraction = 0.9, n.minobsinnode = 5)
  fit
}
est_mux_gbm <- function(df, p, ps){
  Xname = paste0("X", 1:p)
  form = formula(paste("Y ~ ", paste(c(Xname), collapse = " + "), collapse = " + "))
  fit = gbm(form, data = df[df$R == 1, ],
            distribution = "gaussian", n.trees = 1000, shrinkage = 0.05, bag.fraction = 0.9, n.minobsinnode = 5)
  fit
}


### Extract nuisance estimation results 
extract_r_logistic_offset_fit <- function(fit, valid){
  pi_hat = predict(fit,newdata=valid,type="response") 
  log_lambda_hat = predict(fit,newdata=valid,type="link") - log(valid$r_N[1]/(1-valid$r_N[1])) 
  list(pi_hat = pi_hat, log_lambda_hat = log_lambda_hat)
}
extract_r_gbm_offset_fit <- function(fit, valid){
  log_lambda_hat = predict(fit,newdata=valid,type="link") 
  pi_hat = 1/(1+exp(-(log_lambda_hat + log(valid$r_N[1]/(1-valid$r_N[1])))))
  list(pi_hat = pi_hat, log_lambda_hat = log_lambda_hat)
}
mse <- function(t){
  sqrt(mean(t^2))
}

### ATE estimation given nuisance estimation results 
# estimate_dr <- function(r_hat0, e_hat0, muxs1_hat0, muxs0_hat0, mux1_hat0, mux0_hat0, df){
#   R = df$R; Y = df$Y; Treat = df$Treat
#   temp = (mux1_hat0 - mux0_hat0) +
#     ((Treat * R)/(e_hat0 * r_hat0)) * (Y - muxs1_hat0) -
#     (((1 - Treat) * R)/((1 - e_hat0) * r_hat0)) * (Y - muxs0_hat0) +
#     (Treat/e_hat0) * (muxs1_hat0 - mux1_hat0) - ((1-Treat)/(1-e_hat0)) * (muxs0_hat0 - mux0_hat0)
#   mean(temp)
# }
estimate_dr_sn <- function(r_hat0, e_hat0, muxs1_hat0, muxs0_hat0, mux1_hat0, mux0_hat0, df){
  # doubly robust estimator with self-normalization weights: this tends to be more stable than the vanilla
  # estimate_dr function without self-normalization si we use this one in the simulations 
  R = df$R; Y = df$Y; Treat = df$Treat
  IPW1 = mean((Treat * R)/(e_hat0 * r_hat0)); IPW2 = mean(((1 - Treat) * R)/((1 - e_hat0) * r_hat0))
  IPW3 = mean(Treat/e_hat0); IPW4 = mean((1-Treat)/(1-e_hat0))
  
  temp = (mux1_hat0 - mux0_hat0) + 
    ((Treat * R)/(e_hat0 * r_hat0)) * (Y - muxs1_hat0)/IPW1 -  
    (((1 - Treat) * R)/((1 - e_hat0) * r_hat0)) * (Y - muxs0_hat0)/IPW2 + 
    (Treat/e_hat0) * (muxs1_hat0 - mux1_hat0)/IPW3 - ((1-Treat)/(1-e_hat0)) * (muxs0_hat0 - mux0_hat0)/IPW4
  mean(temp)
}
estimate_dr_std_limit_sn <- function(r_hat0, e_hat0, muxs1_hat0, muxs0_hat0, mux1_hat0, mux0_hat0, df, est){
  # estimate the standard error 
  
  R = df$R; Y = df$Y; Treat = df$Treat
  IPW1 = mean((Treat * R)/(e_hat0 * r_hat0)); IPW2 = mean(((1 - Treat) * R)/((1 - e_hat0) * r_hat0))
  
  temp = ((Treat * R)/(e_hat0 * r_hat0)) * (Y - muxs1_hat0)/IPW1 -  
    (((1 - Treat) * R)/((1 - e_hat0) * r_hat0)) * (Y - muxs0_hat0)/IPW2
  sqrt((sum(R)/length(Y))*var(temp))/sqrt(sum(R))
}


