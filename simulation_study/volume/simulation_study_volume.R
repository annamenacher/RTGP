##################################################
### Script: Simulation study (volumetric data) ###
##################################################

# Use this file to get parameter estimates, prediction results, and inference results for GPR + Normal, GPR + Horseshoe, STGP, RTGP-Gibbs, and RTGP-VI.

# Either use array ID from cluster batch submission command OR set to dataset that should
# be evaluated.
dataset = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n_train = 500 
poly_degree = 10
signal = 1
noise = 0.2
effect = 124

#################
### Libraries ###
#################

library(MASS)
library(BayesGPfit)
library(parallel)
library(fields)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(truncnorm)
library(fastBayesReg)
library(STGP)
library(ggplot2)
library(ggpointdensity)
library(xtable)

#################
### Functions ###
#################


plot_2D_funcs = function(f,grids,names=NULL,nrow=NULL,ncol=NULL){
  fdat = data.frame(f,v1=grids[,1],v2=grids[,2])
  if(!is.null(names)){
    names(fdat) = c(names,"v1","v2")
  }
  fdat1 = melt(fdat,id.vars=c("v1","v2"))
  return(ggplot(fdat1,aes(x=v1,y=v2,fill=value))+facet_wrap(~variable,nrow=nrow,ncol=ncol)+coord_equal()+geom_tile()+scale_fill_gradientn(colours =GP.create.cols()))
}

#############
### Paths ###
#############

# Path to folder where data is saved in.
path_data = ""

# Path to folder to save results in.
path_core = ""

path_results = path =  sprintf("%s/N%s_poly%s_signal%s_noise%s_effect%s/", path_core, n_train, poly_degree, signal, noise, effect)
dir.create(file.path(path_results), showWarnings = FALSE)

path_results = path =  sprintf("%s/N%s_poly%s_signal%s_noise%s_effect%s/dataset%s/", path_core, n_train, poly_degree, signal, noise, effect, dataset)
dir.create(file.path(path_results), showWarnings = FALSE)

# Save GPR + Normal & GPR + Horseshoe.
path_results_NM = path =  sprintf("%s/N%s_poly%s_signal%s_noise%s_effect%s/dataset%s/NM/", path_core, n_train, poly_degree, signal, noise, effect, dataset)
dir.create(file.path(path_results_NM), showWarnings = FALSE)

# Save STGP.
path_results_STGP = path =  sprintf("%s/N%s_poly%s_signal%s_noise%s_effect%s/dataset%s/STGP/", path_core, n_train, poly_degree, signal, noise, effect, dataset)
dir.create(file.path(path_results_STGP), showWarnings = FALSE)

# Save RTGP-Gibbs.
path_results_RTGP_Gibbs = path =  sprintf("%s/N%s_poly%s_signal%s_noise%s_effect%s/dataset%s/RTGP_Gibbs/", path_core, n_train, poly_degree, signal, noise, effect, dataset)
dir.create(file.path(path_results_RTGP_Gibbs), showWarnings = FALSE)

# Save RTGP-VI.
path_results_RTGP_VI = path =  sprintf("%s/N%s_poly%s_signal%s_noise%s_effect%s/dataset%s/RTGP_VI/", path_core, n_train, poly_degree, signal, noise, effect, dataset)
dir.create(file.path(path_results_RTGP_VI), showWarnings = FALSE)

# Save tables.
path_results_all =  sprintf("%s/N%s_poly%s_signal%s_noise%s_effect%s/dataset%s/all/", path_core, n_train, poly_degree, signal, noise, effect, dataset)
dir.create(file.path(path_results_all), showWarnings = FALSE)

############
### Data ###
############

# Load dataset.
load(paste0(path_data, "N", toString(n_train), "_poly", toString(poly_degree), "_signal", toString(signal), "_noise", toString(noise), "_effect", toString(effect), "/dataset", toString(dataset), "/data.RData"))

X = data$X 
y = data$y 
v_list = data$v_list
Psi = data$Psi 
lambda_l = data$lambda 
beta_true = data$beta_true
alpha_true = data$alpha_true 
true_sigma2 = data$true_sigma2 
y_test = data$y_test 
img_test = data$img_test 
Z_test = data$Z_test 
p = data$p 
n_train = data$n_train 
n_test = data$n_test  
n = data$n
poly_degree = data$poly_degree 
a = data$a 
b = data$b 
signal = data$signal  
noise = data$noise  
r = data$r 
effect = data$effect 
beta0_true = data$beta0_true 

sqrt_lambda = sqrt(lambda_l)
Bases = Psi*sqrt_lambda
L = nrow(Bases)
Z = cbind(1,t(Bases%*%t(X)))

###########
### GPR ###
###########

# Create list for time results.
result_time = list()

# Set current results path.
path = path_results_NM
# Create list for GPR results.
result_nm = list()

# Start time for baseline GPRs.
time.start = Sys.time()

# Fit GPR + Horseshoe, GPR + Normal, GPR + MFVB.
hs_fit_SOI = fast_horseshoe_lm(y,Z)
nm_fit_SOI = fast_normal_lm(y,Z)
mfvb_fit_SOI = fast_mfvb_normal_lm(y,Z)

# End time for baseline GPRs.
time.end = Sys.time()
# Duration of estimating model parameters of GPRs.
time = time.end - time.start

# Save time result of GPRs.
result_time$time_nm = time

# Concatenate basis coefficients of GPRs.
theta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-1],
                       NM = nm_fit_SOI$post_mean$betacoef[-1],
                       MFVB = mfvb_fit_SOI$post_mean$betacoef[-1])

# Concatenate image coefficients of GPRs.
beta_fit = data.frame(HS = crossprod(Bases,hs_fit_SOI$post_mean$betacoef[-1]),
                      NM = crossprod(Bases,nm_fit_SOI$post_mean$betacoef[-1]),
                      MFVB = crossprod(Bases,mfvb_fit_SOI$post_mean$betacoef[-1]))

# Save results in list.
result_nm$beta_fit = beta_fit
result_nm$theta_fit = theta_fit
result_nm$time = time
save(result_nm, file = sprintf("%sresult_nm.RData", path_results_NM))

# Plot estimated and true image coefficients.
p_plot = plot_2D_funcs(cbind(beta_fit,beta_true),v_list)
ggsave(p_plot, filename=paste0(path_results_NM, "estimated_beta_NM_HS_MFVB.png"), type="cairo-png")

# Plot scatterplot of true beta vs fitted beta of GPR + Normal.
df = data.frame(beta_true = beta_true, beta_method = beta_fit$NM)
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta True') +
  ylab('Beta NM') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
p_plot
ggsave(p_plot, filename=paste0(path_results_NM, "scatterplot_beta_true_vs_beta_NM.png"), type="cairo-png")

# True threshold.
threshold_true = signal/2

# True significance mask.
mask = ifelse(abs(alpha_true)>threshold_true, 1, 0)

# Plot true significance map.
p_plot = plot_2D_funcs(mask, v_list)
ggsave(p_plot, filename=paste0(path, "truth_binary_significance.png"), type="cairo-png")

### GPR + Normal ###

# Evaluate bias of parameters
print('Bias: NM')
(bias_nm_no_effect = mean((beta_fit$NM - beta_true)[mask == 0]))
(bias_nm_effect = mean((beta_fit$NM - beta_true)[mask == 1]))
(bias_nm_total = mean((beta_fit$NM - beta_true)))

# Evaluate MSE of parameters
print('MSE: NM')
(mse_nm_no_effect = mean((beta_fit$NM - beta_true)[mask == 0]^2))
(mse_nm_effect = mean((beta_fit$NM - beta_true)[mask == 1]^2))
(mse_nm_total = mean((beta_fit$NM - beta_true)^2))

# Evaluate inference results
print('Inference results: NM')
t_beta_mean = apply(crossprod(Bases,nm_fit_SOI$mcmc$betacoef[-1,]), 1, mean)
t_beta_sd = apply(crossprod(Bases,nm_fit_SOI$mcmc$betacoef[-1,]), 1, sd)
t = t_beta_mean / t_beta_sd
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr = ifelse(p_fdr < 0.05, 1, 0)
TP = t_bin_fdr + mask
TP = length(TP[TP == 2])
TN = t_bin_fdr + mask
TN = length(TN[TN == 0])
FP = t_bin_fdr - mask
FP = length(FP[FP == 1])
FN = t_bin_fdr - mask
FN = length(FN[FN == -1])
(TPR_nm = TP/(TP + FN))
(TDR_nm = TP/(TP + FP))
(FPR_nm = FP/(FP + TN))
(FDR_nm = FP/(FP + TP))

# Plot estimated beta coefficients of GPR + Normal.
p_plot = plot_2D_funcs(beta_fit$NM, v_list)
ggsave(p_plot, filename=paste0(path, "NM_beta_estimated.png"), type="cairo-png")
# Plot binary significance map of GPR + Normal.
p_plot = plot_2D_funcs(t_bin_fdr, v_list)
ggsave(p_plot, filename=paste0(path, "NM_binary_significance.png"), type="cairo-png")
# Plot test statistics map of GPR + Normal.
p_plot = plot_2D_funcs(t, v_list)
ggsave(p_plot, filename=paste0(path, "NM_test_statistic.png"), type="cairo-png")

# Evaluate R2 (train) of prediction results
print('R2 Train: NM')
pred_nm = nm_fit_SOI$post_mean$betacoef[1] + crossprod(t(X),beta_fit$NM)
(R2_nm = cor(pred_nm,y)^2)

# Evaluate bias (train) of prediction results
print('Bias Prediction Train: NM')
(bias_y_pred_nm = mean(pred_nm - y))

# Evaluate MSE (train) of prediction results
print('MSE Prediction Train: NM')
(mse_y_pred_nm = mean((pred_nm - y)^2))

# Evaluate R2 (test) of prediction results
print('R2 Test: NM')
pred_nm_test = nm_fit_SOI$post_mean$betacoef[1] + crossprod(data$img_test,beta_fit$NM)
(R2_nm_test = cor(pred_nm_test,data$y_test)^2)

# Evaluate bias (test) of prediction results
print('Bias Prediction Test: NM')
(bias_y_pred_nm_test = mean(pred_nm_test - data$y_test))

# Evaluate MSE (test) of prediction results
print('MSE Prediction Test: NM')
(mse_y_pred_nm_test = mean((pred_nm_test - data$y_test)^2))
 

############
### STGP ###
############

# Set current results path.
path = path_results_STGP

sim.STGP.new = function(m = 10){
  
  # Function generates knots on a spatial 2D grid for the kernel convolution approximation of a GP
  # from the STGP package from Jian Kang.
  
  ## Settings
  
  s  <- (1:sqrt(p)) / (sqrt(p)+1)            
  S  <- expand.grid(x=s, y=s)
  m2 <- m*m
  
  ## Set up the basis functions:
  
  knots <- seq(0,1,length=m)
  Grid  <- expand.grid(x = knots, y = knots)
  bw    <- min(dist(Grid))
  
  AA    <- find_neighbors_2D(Grid, bw+0.05)
  A     <- ADJ2D(m2, AA)
  
  D     <- rdist(S,Grid)
  B     <- exp(-0.5*(D/bw)^2)
  B     <- ifelse(D>2.5*bw,0,B)
  
  output <- list(B=B,
                 ADJ=A,
                 voxel=S)
  
  return (output)
}

# Generate knots.
sim_stgp = sim.STGP.new(m=10)

# Save the generated knots for STGP model estimation.
save(sim_stgp, file = sprintf("%sdata_STGP.RData", path))

# Start time for the estimation of STGP parameters.
time.start = Sys.time()

# STGP model estimation via STGP package.
model_stgp = STGP(y,X,sim_stgp$B,sim_stgp$ADJ,t_min=0,t_max=2,show_iters=50, iters = 5000, burn = 2500, Xp = rbind(X, t(data$img_test)))

# End time for the estimation of STGP parameters.
time.end = Sys.time()
# Time needed to estimte parameters of the STGP model.
time = time.end - time.start

# Save the time needed to estimate STGP model.
result_time$time_stgp = time
model_stgp$time = time

# Save the results from the STGP model estimation.
save(model_stgp, file = sprintf("%sresult_STGP.RData", path))

# Beta coefficients (posterior mean of MCMC chain)
beta_stgp = as.vector(model_stgp$beta.mn)
# Beta coefficients (posterior variance of MCMC chain)
beta_var_stgp = as.vector(model_stgp$beta.var)
# Probability of nonzero beta coefficients (of MCMC chain)
prob_beta_nonzero_stgp = as.vector(model_stgp$beta.nz)
# Binary significance of STGP model. 
t_binary = ifelse(prob_beta_nonzero_stgp > 0.5, 1, 0)

# Plot estimated beta map.
p_plot = plot_2D_funcs(beta_stgp, v_list)
ggsave(p_plot, filename=paste0(path, "STGP_beta_estimated.png"), type="cairo-png")
# Plot estimated beta variance map.
p_plot = plot_2D_funcs(beta_var_stgp, v_list)
ggsave(p_plot, filename=paste0(path, "STGP_beta_estimated_variance.png"), type="cairo-png")
# Plot estimated probability of nonzero beta map.
p_plot = plot_2D_funcs(prob_beta_nonzero_stgp, v_list)
ggsave(p_plot, filename=paste0(path, "STGP_prob_beta_nonzero.png"), type="cairo-png")
# Plot estimated binary significance map.
p_plot = plot_2D_funcs(t_binary, v_list)
ggsave(p_plot, filename=paste0(path, "STGP_binary_significance.png"), type="cairo-png")

# Plot scatterplot of true beta vs estimated beta.
df = data.frame(beta_true = beta_true, beta_method = beta_stgp)
      p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
        geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
        geom_abline(intercept = 0) +
        xlim(min(df), max(df)) +
        ylim(min(df), max(df)) +
        xlab('Beta True') +
        ylab('Beta STGP') +
        theme(text = element_text(size=35)) +
        guides(fill = guide_legend(title=''))
      ggsave(p_plot, filename=paste0(path, "scatterplot_beta_true_vs_beta_stgp.png"), type="cairo-png")

# Evaluate bias of parameters
print('Bias: STGP')
(bias_stgp_no_effect = mean((beta_stgp - beta_true)[mask == 0]))
(bias_stgp_effect = mean((beta_stgp - beta_true)[mask == 1]))
(bias_stgp_total = mean((beta_stgp - beta_true)))

# Evaluate MSE of parameters
print('MSE: STGP')
(mse_stgp_no_effect = mean((beta_stgp - beta_true)[mask == 0]^2))
(mse_stgp_effect = mean((beta_stgp - beta_true)[mask == 1]^2))
(mse_stgp_total = mean((beta_stgp - beta_true)^2))

# Evaluate inference results
print('Inference results: STGP')
t_bin_fdr = t_binary
TP = t_bin_fdr + mask
TP = length(TP[TP == 2])
TN = t_bin_fdr + mask
TN = length(TN[TN == 0])
FP = t_bin_fdr - mask
FP = length(FP[FP == 1])
FN = t_bin_fdr - mask
FN = length(FN[FN == -1])
(TPR_stgp = TP/(TP + FN))
(TDR_stgp = TP/(TP + FP))
(FPR_stgp = FP/(FP + TN))
(FDR_stgp = FP/(FP + TP))

# Evaluate R2 (train) of prediction results
print('R2 Train: STGP')
pred_stgp = model_stgp$fit.mn[1:n_train,1]
(R2_stgp = cor(pred_stgp,y)^2)

# Evaluate bias (train) of prediction results
print('Bias Prediction Train: STGP')
(bias_y_pred_stgp = mean(pred_stgp - y))

# Evaluate MSE (train) of prediction results
print('MSE Prediction Train: STGP')
(mse_y_pred_stgp = mean((pred_stgp - y)^2))

# Evaluate R2 (test) of prediction results
print('R2 Test: STGP')
pred_stgp_test = model_stgp$fit.mn[(n_train+1):n,1]
(R2_stgp_test = cor(pred_stgp_test,data$y_test)^2)

# Evaluate bias (test) of prediction results
print('Bias Prediction Test: STGP')
(bias_y_pred_stgp_test = mean(pred_stgp_test - data$y_test))

# Evaluate MSE (test) of prediction results
print('MSE Prediction Test: STGP')
(mse_y_pred_stgp_test = mean((pred_stgp_test - data$y_test)^2))


###############
### RTGP-VI ###
###############

# Maximum number of iterations
n_iter = 10000
# Threshold for ELBO difference
eps = 0.1

# Functions for RTGP-VI model estimation.
phi_function = function(x){
  
  # Standard normal density function
  
  x = (1/sqrt(2*pi)) * exp((-0.5)*x^2)
  return(x)
}

log_likelihood = function(y, Psi, lambda_l, theta, beta0, X, sigma_epsilon, alpha, threshold){
  
  # Log likelihood function of RTGP model
  
  beta = t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta
  beta_thresh = matrix(beta * ifelse(abs(alpha) > threshold, 1, 0), p, 1)
  mu = beta0 + X%*%beta_thresh
  ll = sum(dnorm(y, mu, sigma_epsilon, log=TRUE))
  return(ll)
}

joint_density = function(){
  
  # Joint density function of prior and likelihood of RTGP model
  
  beta = t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta
  beta_thresh = matrix(beta * ifelse(abs(alpha) > threshold, 1, 0), p, 1)
  mu = beta0 + X%*%beta_thresh
  joint = sum(dnorm(y, mu, sigma_epsilon, log=TRUE)) + dnorm(beta0, 0, sigma_beta0, log = T) + 
    sum(dnorm(theta, 0, sigma_beta, log = T)) + 
    sum(dnorm(alpha, beta, sigma_alpha, log = T)) +
    (-1)*(a_beta + 1)*log(sigma_beta^2) - ((b_beta) / sigma_beta^2)  +
    (-1)*(a_alpha + 1)*log(sigma_alpha^2) - ((b_alpha) / sigma_alpha^2)  +
    (-1)*(a_epsilon + 1)*log(sigma_epsilon^2) - ((b_epsilon) / sigma_epsilon^2) 
  return(joint)
}

# Constants and hyperparameters
# Hyperparameters of Hyperprior on variances sigma_alpha^2, sigma_beta^2, sigma_epsilon^2
a_beta = b_beta = 0.001
a_alpha = b_alpha = 0.001
a_epsilon = b_epsilon = 0.001
# Boundaries of Uniform prior on threshold lambda (lower bound: 25% quantile of estimated beta coefficients of GPR + Normal,
# upper bound: 90% quantile of estimated beta coefficients of GPR + Normal)
t_min = max(beta_fit$NM) * 0.25
t_max =  max(beta_fit$NM) * 0.9
# Set equidistant threshold options for Uniform prior on lambda.
threshold_options = seq(t_min + 0.01, t_max, length.out = 10)
# Hyperparameter of variance of prior on intercept beta0.
sigma_beta0 = 10

# Initialisation of parameters
threshold = (t_min + t_max) / 2
theta = matrix(theta_fit$NM,L,1)
beta0 =  nm_fit_SOI$post_mean$betacoef[1]
beta = beta_fit$NM
alpha = beta
sigma_epsilon = sigma_epsilon_prior = sqrt(nm_fit_SOI$post_mean$sigma2_eps)
sigma_alpha = sigma_alpha_prior = 0.1
sigma_beta = sigma_beta_prior = 1

# Start time of RTGP model estimation
time.start = Sys.time()
# Calculate true log likelihood of simulation study.
log_lik_true = sum(dnorm(y, beta0_true + X%*%matrix(beta_true,p,1), sqrt(true_sigma2), log=TRUE))

# Inverse of variance of sigma_alpha^2, sigma_beta^2, sigma_epsilon^2
exp_1_sigma_epsilon2 = 1/sigma_epsilon^2
exp_1_sigma_alpha2 = 1/sigma_alpha^2
exp_1_sigma_beta2 = 1/sigma_beta^2

# Variational posterior variance of beta0 
Sigma_beta0 = 1/(n_train*exp_1_sigma_epsilon2 + (1/sigma_beta0^2))
# Variational posterior variance of theta 
Sigma_theta = solve(exp_1_sigma_epsilon2 * (t(X%*%diag(ifelse(abs(as.vector(alpha))>threshold,1,0),p)%*%t(Psi)%*%diag(sqrt(lambda_l),L))%*%X%*%diag(ifelse(abs(as.vector(alpha))>threshold,1,0),p)%*%t(Psi)%*%diag(sqrt(lambda_l),L) + diag(exp_1_sigma_beta2,L) + exp_1_sigma_alpha2*t(t(Psi)%*%diag(sqrt(lambda_l),L))%*%t(Psi)%*%diag(sqrt(lambda_l),L)))

# Initialisation of weights of mixture of truncated normal distributions on alpha
w = matrix(NA, 3, p)
w[1,] = as.vector(pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2)))
w[2,] = as.vector(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2)))
w[3,] = as.vector(1 - pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)))
w_sum = apply(w, 2, sum)
w = sweep(w, 2, w_sum, FUN="/")

# Initialisation of mixture of truncated normal distributions on alpha 
exp_alpha = matrix(NA, 3, p)
exp_alpha[1,] = beta - sigma_alpha * (phi_function((-threshold - beta) / sigma_alpha) / pnorm((-threshold - beta) / sigma_alpha))
exp_alpha[2,] = beta - sigma_alpha * (phi_function((threshold - beta) / sigma_alpha) - phi_function((-threshold - beta) / sigma_alpha)) / (pnorm((threshold - beta) / sigma_alpha) - pnorm((-threshold - beta) / sigma_alpha))
exp_alpha[3,] = beta + sigma_alpha * phi_function((threshold - beta) / sigma_alpha) / (1 - pnorm((threshold - beta) / sigma_alpha))
  
# Initialisation of difference in ELBO
diff = 100
# Initialisation of ELBO
ELBO = -100000000000000000
# Start iteration counter
iter = 0

# Save VI progression to check continuous increase in ELBO.
log_likelihood_chain = numeric(n_iter)
joint_chain = numeric(n_iter)
diff_chain = numeric(n_iter)
threshold_chain = numeric(n_iter)
ELBO_chain = numeric(n_iter)

# Run until convergence is reached. 
while((diff > eps) && (iter != n_iter)){

  tryCatch({
    
    # Increase iteration count by 1.
    iter = iter + 1
    
    # Save values of beta coefficients and ELBO of previous iteration. 
    beta_old = beta
    ELBO_old = ELBO
    
    # Calculate log likelihood.
    ll = log_likelihood(y, Psi, lambda_l, theta, beta0, X, sigma_epsilon, alpha, threshold)
    # Calculate joint.
    joint = joint_density()
    
    # Update beta0 
    exp_I_abs_alpha_greater_lamda = w[1,] + w[3,]
    Sigma_beta0 = 1/(n_train*exp_1_sigma_epsilon2 + (1/(sigma_beta0^2)))
    beta0 = as.vector( Sigma_beta0 * (exp_1_sigma_epsilon2 * sum(y - X%*%(exp_I_abs_alpha_greater_lamda * beta))))
    
    # Update theta
    Sigma_theta = solve(exp_1_sigma_epsilon2 * t(X %*% diag(as.vector(exp_I_abs_alpha_greater_lamda), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) %*% (X %*% diag(as.vector(exp_I_abs_alpha_greater_lamda), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) + diag(exp_1_sigma_beta2, L) + exp_1_sigma_alpha2*t(t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)) %*% (t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)))
    theta = Sigma_theta %*% t(exp_1_sigma_epsilon2 * t(y - beta0) %*% (X%*%diag(as.vector(exp_I_abs_alpha_greater_lamda),p)%*%t(Psi)%*%diag(as.vector(sqrt(lambda_l)),L)) + exp_1_sigma_alpha2 * t(alpha)%*%t(Psi)%*%diag(as.vector(sqrt(lambda_l)),L))  
    beta = t(Psi)%*%(sqrt(lambda_l) * theta)

    # Update alpha
    exp_beta2 = diag((t(Psi) %*% diag(sqrt(lambda_l),L)) %*% (Sigma_theta + theta %*% t(theta)) %*% t(t(Psi) %*% diag(sqrt(lambda_l),L)))
    c = exp((-0.5) * exp_1_sigma_epsilon2 * 2 * t(t(y-beta0) %*% (-X)) * beta  - 0.5 * exp_1_sigma_epsilon2 * colSums(X^2) * exp_beta2-0.5 * exp_1_sigma_epsilon2 * colSums(sweep(X, 2, beta, "*")*(sweep(-sweep(X, 2, (beta * exp_I_abs_alpha_greater_lamda), "*") , 1, X %*% (beta * exp_I_abs_alpha_greater_lamda), "+"))))
    
    w = matrix(NA, 3, p)
    w[1,] = as.vector(pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) * c
    w[2,] = as.vector(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) 
    w[3,] = as.vector(1 - pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2))) * c 
    w_sum = apply(w, 2, sum)
    w = sweep(w, 2, w_sum, FUN="/")
    
    exp_alpha = matrix(NA, 3, p)
    exp_alpha[1,] = beta - sigma_alpha * (phi_function((-threshold - beta) / sigma_alpha) / pnorm((-threshold - beta) / sigma_alpha))
    exp_alpha[2,] = beta - sigma_alpha * (phi_function((threshold - beta) / sigma_alpha) - phi_function((-threshold - beta) / sigma_alpha)) / (pnorm((threshold - beta) / sigma_alpha) - pnorm((-threshold - beta) / sigma_alpha))
    exp_alpha[3,] = beta + sigma_alpha * phi_function((threshold - beta) / sigma_alpha) / (1 - pnorm((threshold - beta) / sigma_alpha))
    exp_alpha[is.na(exp_alpha)] = 0
    exp_alpha[exp_alpha == Inf] = 0
    exp_alpha[exp_alpha == -Inf] = 0
    
    alpha = matrix(w[1,] * exp_alpha[1,] + w[2,] * exp_alpha[2,] + w[3,] * exp_alpha[3,], p, 1)    #alpha = exp_alpha
    
    # Update threshold
    threshold_t = numeric(length(threshold_options))
    for(k in 1:length(threshold_options)){
      threshold_t[k] = (-0.5) * exp_1_sigma_epsilon2 * sum((y - beta0 - X%*%(beta * ifelse(abs(alpha) > threshold_options[k], 1, 0)))^2)
    }

    logProb = threshold_t - max(threshold_t)
    Prob = exp(logProb)/sum(exp(logProb))
    threshold = sample(threshold_options, size=1, prob=Prob)   

    # Update sigma_epsilon
    w = matrix(NA, 3, p)
    w[1,] = as.vector(pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) * c
    w[2,] = as.vector(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) 
    w[3,] = as.vector(1 - pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2))) * c 
    w_sum = apply(w, 2, sum)
    w = sweep(w, 2, w_sum, FUN="/")
    exp_I_abs_alpha_greater_lamda = w[1,] + w[3,]
    var_beta_alpha = exp_beta2 * exp_I_abs_alpha_greater_lamda - beta^2 * exp_I_abs_alpha_greater_lamda^2
    exp_beta_alpha_X = sum(sweep(X^2, 2, var_beta_alpha, "*")) + sum(diag((X%*%(beta * exp_I_abs_alpha_greater_lamda)) %*% t(X%*%(beta * exp_I_abs_alpha_greater_lamda))))
    
    # Update sigma_epsilon
    exp_1_sigma_epsilon2 = as.vector((n_train + a_epsilon) / (b_epsilon + 0.5 * (t(y)%*%y + n_train * (Sigma_beta0 + beta0^2) + exp_beta_alpha_X - 2 * t(y-beta0)%*%(X%*%(beta * exp_I_abs_alpha_greater_lamda)) - sum(2* beta0 * y))))
    sigma_epsilon = as.vector(sqrt(1/exp_1_sigma_epsilon2))
    
    # Update sigma_beta
    exp_1_sigma_beta2 = as.vector((L/2 + a_beta) / (0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta))
    sigma_beta = as.vector(sqrt(1/exp_1_sigma_beta2))
    
    # Update sigma_alpha
    var_alpha = matrix(NA, 3, p)
    var_alpha[1,] = sigma_alpha^2 * (1 - ((-threshold - beta)/sigma_alpha) * phi_function(((-threshold - beta)/sigma_alpha)) / (pnorm(((-threshold - beta)/sigma_alpha)) + 10^(-10)) - (phi_function(((-threshold - beta)/sigma_alpha)) / (pnorm(((-threshold - beta)/sigma_alpha)) + 10^(-10)))^2)
    var_alpha[2,] = sigma_alpha^2 * (1 - (((threshold - beta)/sigma_alpha) * phi_function(((threshold - beta)/sigma_alpha)) - ((-threshold - beta)/sigma_alpha) * phi_function(((-threshold - beta)/sigma_alpha))) / (pnorm(((threshold - beta)/sigma_alpha)) - pnorm(((-threshold - beta)/sigma_alpha)) + 10^(-10)) - ((phi_function(((threshold - beta)/sigma_alpha)) - phi_function(((-threshold - beta)/sigma_alpha))) / (pnorm(((threshold - beta)/sigma_alpha)) - pnorm(((-threshold - beta)/sigma_alpha)) + 10^(-10)))^2)
    var_alpha[3,] = sigma_alpha^2 * (1 + (((threshold - beta)/sigma_alpha)) * phi_function(((threshold - beta)/sigma_alpha)) / (1 - pnorm(((threshold - beta)/sigma_alpha)) + 10^(-10)) - (phi_function(((threshold - beta)/sigma_alpha)) / (1 - pnorm(((threshold - beta)/sigma_alpha)) + 10^(-10)))^2)
    
    exp_alpha2 = w[1,] * (var_alpha[1,] + exp_alpha[1,]^2) + w[2,] * (var_alpha[2,] + exp_alpha[2,]^2) + w[3,] * (var_alpha[3,] + exp_alpha[3,]^2) 
    exp_1_sigma_alpha2 = (p/2 + a_alpha) / (b_alpha + 0.5 * sum((alpha - beta)^2))
    sigma_alpha = as.vector(sqrt(1/exp_1_sigma_alpha2))
    
    beta_thresholded =  as.vector( beta * ifelse(abs(as.vector(alpha))>threshold, 1, 0) )
    alpha_thresholded = as.vector( ifelse(abs(alpha)>threshold, 1, 0) )

    log_likelihood_chain[(iter)] = ll
    joint_chain[iter] = joint
    threshold_chain[iter] = threshold
    
    # Calculate ELBO
    exp_ln_sigma_alpha2 = log(0.5 * sum((alpha-beta)^2) + b_alpha) - digamma(p/2 + a_alpha)
    exp_ln_sigma_epsilon2 = log(0.5 * (t(y)%*%y + n_train * (Sigma_beta0 + beta0^2) + exp_beta_alpha_X - 2 * t(y-beta0)%*%(X%*%(beta * exp_I_abs_alpha_greater_lamda)) - sum(2* beta0 * y)) + b_epsilon) - digamma(n_train/2 + a_epsilon)
    exp_ln_sigma_beta2 = log(0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta) - digamma(L/2 + a_beta)
    part = log(pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2)))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part1 = sum(w[1,] * part) 
    part = -(-threshold - beta)*sqrt(exp_1_sigma_alpha2) * phi_function((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) / (2*pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part1 = alpha_part1 + sum(w[1,] * part) 
    part = log(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2))) - pnorm((-threshold-beta)*sqrt(exp_1_sigma_alpha2))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part2 = sum(w[2,] * part)
    part = ((-threshold - beta)*sqrt(exp_1_sigma_alpha2) * phi_function((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) - (threshold - beta)*sqrt(exp_1_sigma_alpha2)*phi_function((threshold - beta)*sqrt(exp_1_sigma_alpha2))) / (2 * (pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2))) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part2 = alpha_part2 + sum(w[2,] * part)
    part = log(1 - pnorm((threshold-beta)*sqrt(exp_1_sigma_alpha2)))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part3 = sum(w[3,] * part)
    part = ((threshold-beta)*sqrt(exp_1_sigma_alpha2)*phi_function((threshold-beta)*sqrt(exp_1_sigma_alpha2))) / (2*(1-pnorm((threshold-beta)*sqrt(exp_1_sigma_alpha2))) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part3 = alpha_part3 + sum(w[3,] * part)
    
    ELBO = - (n_train/2) * log(2*pi) - (n_train/2) * exp_ln_sigma_epsilon2 - 0.5 * exp_1_sigma_epsilon2 * (t(y)%*%y + n_train * (Sigma_beta0 + beta0^2) + exp_beta_alpha_X - 2 * t(y-beta0)%*%(X%*%(beta * exp_I_abs_alpha_greater_lamda)) - sum(2* beta0 * y)) - 
      0.5 * log(2*pi) - 0.5 * log(sigma_beta0^2) - (1/sigma_beta0^2) * (Sigma_beta0 + beta0^2) - (L/2) * log(2*pi) - (L/2) * exp_ln_sigma_beta2 -
      0.5 * exp_1_sigma_beta2 * sum(diag(Sigma_theta) + theta^2) - (p/2) * log(2*pi) - (p/2) * exp_ln_sigma_alpha2 - 0.5 * exp_1_sigma_alpha2 * sum((alpha-beta)^2) +
      a_epsilon * log(b_epsilon) + (a_epsilon + 1) * exp_ln_sigma_epsilon2 - b_epsilon * exp_1_sigma_epsilon2 + a_beta * log(b_beta) +
      (a_beta + 1) * exp_ln_sigma_beta2 - b_beta * exp_1_sigma_beta2 + a_alpha * log(b_alpha) + (a_alpha + 1) * exp_ln_sigma_alpha2 -
      b_alpha * exp_1_sigma_alpha2 + 0.5 * log(2*pi) + 0.5 * log(Sigma_beta0) + 0.5 + (L/2) * log(2*pi) + (L/2) * log(det(Sigma_theta) + 10^(-10)) + 
      (L/2) + (p/2) * log(2*pi) + (p/2) * exp_ln_sigma_alpha2  + alpha_part1 + 
      alpha_part2 + alpha_part3 - (n_train/2 + a_epsilon) * log(b_epsilon + 0.5 * sum((y - beta0 - X%*%beta_thresholded)^2)) + 
      ((n_train/2 + a_epsilon) + 1) * exp_ln_sigma_epsilon2 +
      (b_epsilon + 0.5 * sum((y - beta0 - X%*%beta_thresholded)^2)) * exp_1_sigma_epsilon2 - (L/2 + a_beta) * log(0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta) +
      (L/2 + a_beta + 1) * exp_ln_sigma_beta2 + (0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta) * exp_1_sigma_beta2 -
      (p/2 + a_alpha) * log(0.5 * sum((alpha - beta)^2) + b_alpha)  + (p/2 + a_alpha + 1) * exp_ln_sigma_alpha2 +
      (0.5 * sum((alpha - beta)^2) + b_alpha) * exp_1_sigma_alpha2
    
    ELBO_chain[iter] = ELBO
    diff = as.vector(ELBO) - as.vector(ELBO_old)
    diff_chain[iter] = diff
  
    # Save results and plot progress every 10 iterations.
    if(iter %% 10 == 0){
      print(iter)
      print(diff)

      result = list()
      result$beta0 = beta0
      result$theta = theta
      result$Sigma_theta = Sigma_theta
      result$beta = beta
      result$alpha = alpha
      result$sigma_epsilon = sigma_epsilon
      result$sigma_beta = sigma_beta
      result$sigma_alpha = sigma_alpha
      result$beta_thresholded = beta_thresholded
      result$alpha_thresholded = alpha_thresholded
      result$log_likelihood_chain = log_likelihood_chain
      result$joint_chain = joint_chain
      result$diff_chain = diff_chain
      result$threshold_chain = threshold_chain
      result$threshold = threshold 
      result$threshold_options = threshold_options
      result$threshold_probabilities = Prob
      result$ELBO_chain = ELBO_chain
      
      save(result, file = sprintf("%sresult_rtgp_vi.RData", path_results_RTGP_VI))

      path = path_results_RTGP_VI

      p_plot = plot_2D_funcs(beta , v_list)
      ggsave(p_plot, filename=paste0(path, "RTGP_VI_mean_beta_unthresholded.png"), type="cairo-png")
      
      p_plot = plot_2D_funcs(beta_thresholded , v_list)
      ggsave(p_plot, filename=paste0(path, "RTGP_VI_mean_beta_thresholded.png"), type="cairo-png")
      
      p_plot = plot_2D_funcs(alpha_thresholded , v_list)
      ggsave(p_plot, filename=paste0(path, "RTGP_VI_prob_alpha_greater_threshold.png"), type="cairo-png")
      
      df = data.frame('joint' = joint_chain[1:iter],
                      'Iteration' = 1:iter)
      p_plot = ggplot(df, aes(x = Iteration, y = joint)) + 
        geom_line(aes()) +
        theme_minimal() + 
        theme(panel.background = element_rect(fill = "white",
                                              colour = "white",
                                              size = 0.5, linetype = "solid"),
              plot.background = element_rect(fill = "white")) +
        ggtitle('Log-Joint')
      ggsave(p_plot, filename=paste0(path, "RTGP_VI_log_joint_with_true_joint.png"), type="cairo-png")
      
      df = data.frame('log_likelihood' = log_likelihood_chain[1:iter],
                      'Iteration' = 1:iter)
      p_plot = ggplot(df, aes(x = Iteration, y = log_likelihood)) + 
        geom_hline(yintercept = log_lik_true, color = 'red') +
        geom_line(aes()) +
        theme_minimal() + 
        theme(panel.background = element_rect(fill = "white",
                                              colour = "white",
                                              size = 0.5, linetype = "solid"),
              plot.background = element_rect(fill = "white")) +
        ggtitle('Log-Likelihood')
      ggsave(p_plot, filename=paste0(path, "RTGP_VI_log_likelihood_with_true_log_lik.png"), type="cairo-png")
      
      df = data.frame('diff' = diff_chain[2:iter],
                      'Iteration' = 2:iter)
      p_plot = ggplot(df, aes(x = Iteration, y = diff)) + 
        geom_line(aes()) +
        theme_minimal() + 
        theme(panel.background = element_rect(fill = "white",
                                              colour = "white",
                                              size = 0.5, linetype = "solid"),
              plot.background = element_rect(fill = "white")) +
        ggtitle('Difference in ELBO')
      ggsave(p_plot, filename=paste0(path, "RTGP_VI_diff_ELBO.png"), type="cairo-png")
      
      df = data.frame('threshold' = threshold_chain[1:iter],
                      'Iteration' = 1:iter)
      p_plot = ggplot(df, aes(x = Iteration, y = threshold)) + 
        geom_line(aes()) +
        theme_minimal() + 
        theme(panel.background = element_rect(fill = "white",
                                              colour = "white",
                                              size = 0.5, linetype = "solid"),
              plot.background = element_rect(fill = "white")) +
        ggtitle('Threshold')
      ggsave(p_plot, filename=paste0(path, "RTGP_VI_threshold.png"), type="cairo-png")
    
    
    df = data.frame('ELBO' = ELBO_chain[2:iter],
                    'Iteration' = 2:iter)
    p_plot = ggplot(df, aes(x = Iteration, y = ELBO)) + 
      geom_line(aes()) +
      theme_minimal() + 
      theme(panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            plot.background = element_rect(fill = "white")) +
      ggtitle('ELBO')
    ggsave(p_plot, filename=paste0(path, "RTGP_VI_ELBO.png"), type="cairo-png")
  }
  
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  time.end = Sys.time()
  time = time.end - time.start
  result_time$time_rtgp_vi = time

}

# Save final results once algorithm has converged.
result = list()
result$beta0 = beta0
result$theta = theta
result$Sigma_theta = Sigma_theta
result$beta = beta
result$alpha = alpha
result$sigma_epsilon = sigma_epsilon
result$sigma_beta = sigma_beta
result$sigma_alpha = sigma_alpha
result$beta_thresholded = beta_thresholded
result$alpha_thresholded = alpha_thresholded
result$log_likelihood_chain = log_likelihood_chain
result$joint_chain = joint_chain
result$diff_chain = diff_chain
result$time = time
result$threshold_chain = threshold_chain
result$threshold = threshold 
result$threshold_options = threshold_options
result$threshold_probabilities = Prob
result$ELBO_chain = ELBO_chain

save(result, file = sprintf("%sresult_rtgp_vi.RData", path_results_RTGP_VI))


p_plot = plot_2D_funcs(beta , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_VI_mean_beta_unthresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(beta_thresholded , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_VI_mean_beta_thresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(alpha_thresholded , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_VI_prob_alpha_greater_threshold.png"), type="cairo-png")

df = data.frame('joint' = joint_chain[1:iter],
                'Iteration' = 1:iter)
p_plot = ggplot(df, aes(x = Iteration, y = joint)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white")) +
  ggtitle('Log-Joint')
ggsave(p_plot, filename=paste0(path, "RTGP_VI_log_joint_with_true_joint.png"), type="cairo-png")

df = data.frame('log_likelihood' = log_likelihood_chain[1:iter],
                'Iteration' = 1:iter)
p_plot = ggplot(df, aes(x = Iteration, y = log_likelihood)) + 
  geom_hline(yintercept = log_lik_true, color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white")) +
  ggtitle('Log-Likelihood')
ggsave(p_plot, filename=paste0(path, "RTGP_VI_log_likelihood_with_true_log_lik.png"), type="cairo-png")

df = data.frame('diff' = diff_chain[2:iter],
                'Iteration' = 2:iter)
p_plot = ggplot(df, aes(x = Iteration, y = diff)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white")) +
  ggtitle('Difference in ELBO')
ggsave(p_plot, filename=paste0(path, "RTGP_VI_diff_ELBO.png"), type="cairo-png")

df = data.frame('threshold' = threshold_chain[1:iter],
                'Iteration' = 1:iter)
p_plot = ggplot(df, aes(x = Iteration, y = threshold)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white")) +
  ggtitle('Threshold')
ggsave(p_plot, filename=paste0(path, "RTGP_VI_threshold.png"), type="cairo-png")


df = data.frame('ELBO' = ELBO_chain[2:iter],
                'Iteration' = 2:iter)
p_plot = ggplot(df, aes(x = Iteration, y = ELBO)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white")) +
  ggtitle('ELBO')
ggsave(p_plot, filename=paste0(path, "RTGP_VI_ELBO.png"), type="cairo-png")


# Evaluate bias of parameters
print('Bias: RTGP-VI')
(bias_rtgp_vi_no_effect = mean((beta_thresholded - beta_true)[mask == 0]))
(bias_rtgp_vi_effect = mean((beta_thresholded  - beta_true)[mask == 1]))
(bias_rtgp_vi_total = mean((beta_thresholded  - beta_true)))

# Evaluate MSE of parameters
print('MSE: RTGP-VI')
(mse_rtgp_vi_no_effect = mean((beta_thresholded  - beta_true)[mask == 0]^2))
(mse_rtgp_vi_effect = mean((beta_thresholded  - beta_true)[mask == 1]^2))
(mse_rtgp_vi_total = mean((beta_thresholded  - beta_true)^2))

# Evaluate inference results
print('Inference results: RTGP-VI')
t_bin_fdr = alpha_thresholded
TP = t_bin_fdr + mask
TP = length(TP[TP == 2])
TN = t_bin_fdr + mask
TN = length(TN[TN == 0])
FP = t_bin_fdr - mask
FP = length(FP[FP == 1])
FN = t_bin_fdr - mask
FN = length(FN[FN == -1])
(TPR_rtgp_vi = TP/(TP + FN))
(TDR_rtgp_vi = TP/(TP + FP))
(FPR_rtgp_vi = FP/(FP + TN))
(FDR_rtgp_vi = FP/(FP + TP))

p_plot = plot_2D_funcs(t_bin_fdr, v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_VI_binary_significance.png"), type="cairo-png")

# Evaluate R2 (train) of prediction results
print('R2 Train: RTGP-VI')
pred_rtgp_vi = beta0 + X%*%matrix(beta_thresholded , length(beta_true), 1)
(R2_rtgp_vi = cor(pred_rtgp_vi,y)^2)

# Evaluate bias (train) of prediction results
print('Bias Prediction Train: RTGP-VI')
(bias_y_pred_rtgp_vi = mean(pred_rtgp_vi - y))

# Evaluate MSE (train) of prediction results
print('MSE Prediction Train: RTGP-VI')
(mse_y_pred_rtgp_vi = mean((pred_rtgp_vi - y)^2))

# Evaluate R2 (test) of prediction results
print('R2 Test: RTGP-VI')
pred_rtgp_vi_test = beta0 + t(data$img_test)%*%matrix(beta_thresholded , length(beta_true), 1)
(R2_rtgp_vi_test = cor(pred_rtgp_vi_test,data$y_test)^2)

# Evaluate bias (test) of prediction results
print('Bias Prediction Test: RTGP-VI')
(bias_y_pred_rtgp_vi_test = mean(pred_rtgp_vi_test - data$y_test))

# Evaluate MSE (test) of prediction results
print('MSE Prediction Test: RTGP-VI')
(mse_y_pred_rtgp_vi_test = mean((pred_rtgp_vi_test - data$y_test)^2))



##################
### RTGP-Gibbs ###
##################

# Number of MCMC iterations.
n_iter = 5000
# Number of burn-in iterations.
burn_in = 2500

# Functions 
y_tilde_function = function(j){
  
  # Intermediary quantity to calculate thresholded beta * input image.
  
  x = beta[j] * ifelse(abs(alpha[j]) > threshold, 1, 0) * X[,j] 
  return(x)
}

M1_function = function(j){
  
  # Intermediary quantity to calculate thresholded beta * input image * beta.
  
  x = y_tilde[,j] * X[,j] * beta[j]
  return(x)
}

M2_function = function(j){
  
  # Intermediary quantity to calculate sum across subjects of input image^2.
  
  x = sum(X[,j]^2)
  return(x)
}

log_likelihood = function(y, Psi, lambda_l, theta, beta0, X, sigma_epsilon, alpha, threshold){
  
  # Function to calculate log-likelihood.
  
  beta = t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta
  beta_thresh = matrix(beta * ifelse(abs(alpha) > threshold, 1, 0), p, 1)
  mu = beta0 + X%*%beta_thresh
  ll = sum(dnorm(y, mu, sigma_epsilon, log=TRUE))
  return(ll)
}

joint_density = function(){
  
  # Function to calculate joint consisting of prior and likelihood.
  
  beta = t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta
  beta_thresh = matrix(beta * ifelse(abs(alpha) > threshold, 1, 0), p, 1)
  mu = beta0 + X%*%beta_thresh
  joint = sum(dnorm(y, mu, sigma_epsilon, log=TRUE)) + dnorm(beta0, 0, sigma_beta0, log = T) + 
    sum(dnorm(theta, 0, sigma_beta, log = T)) + 
    sum(dnorm(alpha, beta, sigma_alpha, log = T)) +
  (-1)*(a_beta + 1)*log(sigma_beta^2) - ((b_beta) / sigma_beta^2)  +
  (-1)*(a_alpha + 1)*log(sigma_alpha^2) - ((b_alpha) / sigma_alpha^2)  +
  (-1)*(a_epsilon + 1)*log(sigma_epsilon) - ((b_epsilon) / sigma_epsilon)  
  return(joint)
}
 
sample_alpha = function(j){
  
  # Function to sample alpha from mixture distribution with the estimated probability w.
  
  alpha_j = sample(x = a_all[,j], size = 1, replace = F, prob = w[,j])
  return(alpha_j)
}

# Constants and hyperparameters.
# Hyperparameters of priors on sigma_alpha, sigma_beta, sigma_epsilon.
a_beta = b_beta = 0.001
a_alpha = b_alpha = 0.001
a_epsilon = b_epsilon = 0.001
# Hyperparameter of prior on intercept beta0.
sigma_beta0 = 10
# Boundaries of Uniform prior on threshold lambda (lower bound: 25% quantile of estimated beta coefficients of GPR + Normal,
# upper bound: 90% quantile of estimated beta coefficients of GPR + Normal)
t_min = max(beta_fit$NM) * 0.25
t_max =  max(beta_fit$NM) * 0.9
# Set equidistant threshold options for Uniform prior on lambda.
threshold_options = seq(t_min + 0.01, t_max, length.out = 10)

# Initialisation of parameters drawn from prior distributions or set to the GPR + Normal result.
threshold = runif(1, t_min, t_max)
sigma_epsilon = sigma_epsilon_prior = rgamma(1, shape = 1, rate = 1)
sigma_epsilon = sigma_epsilon_prior = 1/sqrt(sigma_epsilon)
sigma_alpha = sigma_alpha_prior = rgamma(1, shape = 1, rate = 1)
sigma_alpha = sigma_alpha_prior = 1/sqrt(sigma_alpha)
sigma_beta = sigma_beta_prior = rgamma(1, shape = 1, rate = 1)
sigma_beta = sigma_beta_prior = 1/sqrt(sigma_beta)
theta = matrix(theta_fit$NM + rnorm(L, 0, 0.01),L,1)
beta0 =  rnorm(1, 0, sigma_beta0)
beta = beta_fit$NM + rnorm(p, 0, 0.01)
alpha = rnorm(p, beta, sigma_alpha)
beta_thresholded = beta * ifelse(abs(alpha)>threshold,1,0)
alpha_thresholded = ifelse(abs(alpha)>threshold,1,0)

# Create and save a list of initialised parameter values.
init_list = list()
init_list$a_beta = a_beta
init_list$a_epsilon = a_epsilon
init_list$a_alpha = a_alpha
init_list$threshold = threshold
init_list$sigma_epsilon = sigma_epsilon
init_list$sigma_beta = sigma_beta
init_list$sigma_alpha = sigma_alpha
init_list$theta = theta
init_list$theta = beta
init_list$theta = alpha
init_list$theta = beta0
init_list$beta_thresholded = beta_thresholded
init_list$alpha_thresholded = alpha_thresholded
save(init_list, file = sprintf("%sinit_list_rtgp_gibbs.RData", path_results ))

# Create quantities to save MCMC chains.
sigma_beta_chain = numeric(n_iter)
sigma_alpha_chain = numeric(n_iter)
sigma_epsilon_chain = numeric(n_iter)
alpha_chain = matrix(NA, nrow = (n_iter), ncol = p)
theta_chain = matrix(NA, nrow = (n_iter), ncol = L)
beta_chain = matrix(NA, nrow = (n_iter), ncol = p)
beta_thresholded_chain = matrix(NA, nrow = (n_iter), ncol = p)
alpha_thresholded_chain = matrix(NA, nrow = (n_iter), ncol = p)
threshold_chain  = numeric(n_iter)
beta0_chain = numeric(n_iter)
log_likelihood_chain = numeric(n_iter )
joint_chain = numeric(n_iter)

# Start time for RTGP model estimation.
time.start = Sys.time()
# Calculate true log likelihood.
log_lik_true = sum(dnorm(y, beta0_true + X%*%matrix(beta_true,p,1), sqrt(true_sigma2), log=TRUE))

# Set path to RTGP-Gibbs path. 
path = path_results_RTGP_Gibbs

# Run for set amount of iterations.
for(iter in 1:n_iter){
  
  tryCatch({
    
    if(iter %% 50 == 0){
      print(paste0('Iteration: ', toString(iter)))
    }
    
    # Calculate log-likelihood.
    ll = log_likelihood(y, Psi, lambda_l, theta, beta0, X, sigma_epsilon, alpha, threshold)
    # Calculate joint density.
    joint = joint_density()
    
    # Draw beta0.
    beta0_sd = sqrt(1/((n_train/sigma_epsilon^2) + (1/sigma_beta0^2)))
    beta0_mean = beta0_sd^2 * (1/sigma_epsilon^2) * sum(y - X %*% diag(as.vector(ifelse(abs(alpha) > threshold, 1, 0)), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta)
    beta0 = rnorm(1, beta0_mean, beta0_sd)
    
    # Draw theta.
    theta_sigma = solve((1/sigma_epsilon^2) * t(X %*% diag(as.vector(ifelse(abs(alpha) > threshold, 1, 0)), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) %*% (X %*% diag(as.vector(ifelse(abs(alpha) > threshold, 1, 0)), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) + diag((1/sigma_beta^2), L) + (1/sigma_alpha^2)*t(t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)) %*% (t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)))
    theta_mean = theta_sigma %*% t((1/sigma_epsilon^2) * t(y - beta0) %*% (X %*% diag(as.vector(ifelse(abs(alpha) > threshold, 1, 0)), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) + (1/sigma_alpha^2) * t(alpha)%*%(t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)))
    theta = matrix(mvrnorm(1, theta_mean, theta_sigma), L, 1)
    
    # Calculate beta through theta.
    beta = t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta
    
    # Draw alpha.
    sum_y_tilde = apply(matrix(1:p, p, 1), 1, y_tilde_function)
    y_tilde = apply(matrix(1:p, p, 1), 1, function(j) y - rowSums(sum_y_tilde) + sum_y_tilde[,j] - beta0)
    M1 = colSums(apply(matrix(1:p, p, 1), 1, M1_function))
    M2 = apply(matrix(1:p, p, 1), 1, M2_function)
    
    c = exp((-0.5)*(1/sigma_epsilon^2) * (beta^2 * M2 - 2*M1))
    w = matrix(NA, 3, p)
    w[1,] = as.vector(pnorm(-(threshold + beta)/sigma_alpha))*c
    w[2,] = as.vector(pnorm((threshold - beta)/sigma_alpha) - pnorm(-(threshold + beta)/sigma_alpha))
    w[3,] = as.vector(1 - pnorm((threshold - beta)/sigma_alpha))*c
    w_sum = apply(w, 2, sum)
    w = sweep(w, 2, w_sum, FUN="/")
    
    a_all = matrix(NA, 3, p)
    a_all[1,] = rtruncnorm(1, a=-Inf, b=-threshold, mean = as.vector(beta), sd = sigma_alpha)
    a_all[2,] = rtruncnorm(1, a=-threshold, b=threshold, mean = as.vector(beta), sd = sigma_alpha)
    a_all[3,] =  rtruncnorm(1, a=-threshold, b=Inf, mean = as.vector(beta), sd = sigma_alpha)

    alpha = apply(matrix(1:p,p,1), 1, sample_alpha)
    alpha = matrix(alpha, p, 1)
    
    # Draw sigma_beta.
    sigma_beta = rgamma(1, shape = (a_beta+L/2), rate = (b_beta + 0.5*(sum(theta^2))))
    sigma_beta = 1/sqrt(sigma_beta)
    
    # Draw sigma_epsilon.
    sigma_epsilon = rgamma(1, shape = (a_epsilon+n_train/2), rate = (b_epsilon + 0.5*(sum((y - X %*% diag(as.vector(ifelse(abs(alpha) > threshold, 1, 0)), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta - beta0)^2))))
    sigma_epsilon = 1/sqrt(sigma_epsilon)
    
    # Draw sigma_alpha.
    sigma_alpha = rgamma(1, shape = (a_alpha+p/2), rate = (b_alpha + 0.5*sum((alpha - t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L) %*% theta)^2)))
    sigma_alpha = 1/sqrt(sigma_alpha)
   
    # Draw threshold.
    Q = numeric(length(threshold_options))
    for(k in 1:length(threshold_options)){
      Q[k] = (1/(2*(sigma_epsilon^2)))*sum((y - beta0 - X%*%matrix(beta * ifelse(abs(alpha)>threshold_options[k], 1, 0),p,1))^2)
    }
    
    logProb <- -Q
    logProb <- logProb - max(logProb)
    Prob <- exp(logProb)/sum(exp(logProb))
    threshold <- sample(threshold_options, size=1, prob=Prob)

    # Calculate thresholded quantities.
    beta_thresholded = beta * ifelse(abs(alpha)>threshold, 1, 0)
    alpha_thresholded = ifelse(abs(alpha)>threshold, 1, 0)
    
    # Save estimated paramters of MCMC iteration in chain.
    sigma_alpha_chain[(iter)] = sigma_alpha
    sigma_beta_chain[(iter)] = sigma_beta
    sigma_epsilon_chain[(iter)] = sigma_epsilon
    beta0_chain[(iter)] = as.vector(beta0)
    alpha_chain[(iter),] = as.vector(alpha)
    theta_chain[(iter),] = as.vector(theta)
    beta_chain[(iter),] = as.vector(beta)
    beta_thresholded_chain[(iter),] = as.vector(beta_thresholded)
    alpha_thresholded_chain[(iter),] = as.vector(alpha_thresholded)
    threshold_chain[(iter)] = threshold
    log_likelihood_chain[(iter)] = ll
    joint_chain[(iter)] = joint 

# Save results and plot progress every 100 iterations. 
if(iter %% 100 == 0){
   
  time.end = Sys.time()
  time = time.end - time.start
  result_time$time_rtgp = time
  result = list()
  result$sigma_alpha_chain = sigma_alpha_chain
  result$sigma_beta_chain = sigma_beta_chain
  result$sigma_epsilon_chain = sigma_epsilon_chain
  result$beta0_chain = beta0_chain
  result$alpha_chain = alpha_chain
  result$theta_chain = theta_chain
  result$beta_chain = beta_chain
  result$beta_thresholded_chain = beta_thresholded_chain
  result$alpha_thresholded_chain = alpha_thresholded_chain
  result$log_likelihood = log_likelihood_chain
  result$log_lik_true = log_lik_true
  result$joint = joint_chain
  result$threshold_chain = threshold_chain
  result$time = time
save(result, file = sprintf("%sresult_rtgp_gibbs.RData", path_results_RTGP_Gibbs ))

p_plot = plot_2D_funcs(apply(result$beta_chain[1:(iter),],2,mean) , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_mean_beta_unthresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(apply(result$beta_thresholded_chain[1:(iter),], 2, mean) , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_mean_beta_thresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(apply(result$beta_chain[1:(iter),],2,sd)^2 , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_var_beta_unthresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(apply(result$beta_thresholded_chain[1:(iter),],2,sd)^2 , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_var_beta_thresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(apply(alpha_thresholded_chain[1:(iter),],2,mean) , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_prob_alpha_greater_threshold.png"), type="cairo-png")

idx_beta = which(beta_true > 0, arr.ind = T)

df = data.frame('joint' = result$joint[1:(iter)],
                'Iteration' = 1:(iter))
p_plot = ggplot(df, aes(x = Iteration, y = joint)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('Log-Joint')
ggsave(p_plot, filename=paste0(path, "RTGP_log_joint_with_true_joint.png"), type="cairo-png")

df = data.frame('log_likelihood' = result$log_likelihood[1:(iter)],
                'Iteration' = 1:(iter))
p_plot = ggplot(df, aes(x = Iteration, y = log_likelihood)) + 
  geom_hline(yintercept = result$log_lik_true, color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('Log-Likelihood')
ggsave(p_plot, filename=paste0(path, "RTGP_log_likelihood_with_true_log_lik.png"), type="cairo-png")

## Trace plot (threshold)
df = data.frame('threshold' = result$threshold_chain[1:(iter)],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = threshold)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('threshold')
ggsave(p_plot, filename=paste0(path, "RTGP_threshold.png"), type="cairo-png")

## Trace plots (unthresholded betas)
df = data.frame('Beta' = result$beta_chain[1:(iter),idx_beta[1]],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = Beta)) + 
  geom_hline(yintercept = beta_true[idx_beta[1]], color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle(paste0('Beta (Voxel with Effect:', toString(idx_beta[1]) ,')'))
ggsave(p_plot, filename=paste0(path, "RTGP_beta_unthresholded_effect_voxel", toString(idx_beta[1]), ".png"), type="cairo-png")

df = data.frame('Beta' = result$beta_chain[1:(iter),idx_beta[100]],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = Beta)) + 
  geom_hline(yintercept = beta_true[idx_beta[100]], color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle(paste0('Beta (Voxel with Effect:', toString(idx_beta[100]) ,')'))
ggsave(p_plot, filename=paste0(path, "RTGP_beta_unthresholded_effect_voxel", toString(idx_beta[100]), ".png"), type="cairo-png")

## Trace plots (alpha)
df = data.frame('Alpha' = result$alpha_chain[1:(iter),idx_beta[1]],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = Alpha)) + 
  geom_hline(yintercept = alpha_true[idx_beta[1]], color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle(paste0('Alpha (Voxel with Effect:', toString(idx_beta[1]) ,')'))
ggsave(p_plot, filename=paste0(path, "RTGP_alpha_effect_voxel", toString(idx_beta[1]), ".png"), type="cairo-png")

df = data.frame('Alpha' = result$alpha_chain[1:(iter),idx_beta[100]],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = Alpha)) + 
  geom_hline(yintercept = alpha_true[idx_beta[1]], color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle(paste0('Alpha (Voxel with Effect:', toString(idx_beta[100]) ,')'))
ggsave(p_plot, filename=paste0(path, "RTGP_alpha_effect_voxel", toString(idx_beta[100]), ".png"), type="cairo-png")

## Trace plot (sigma_epsilon)
df = data.frame('sigma_epsilon' = result$sigma_epsilon_chain[1:(iter)],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = sigma_epsilon)) + 
  geom_hline(yintercept = sqrt(true_sigma2), color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('sigma_epsilon')
ggsave(p_plot, filename=paste0(path, "RTGP_simga_epsilon.png"), type="cairo-png")

## Trace plot (sigma_beta)
df = data.frame('sigma_beta' = result$sigma_beta_chain[1:(iter)],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = sigma_beta)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('sigma_beta')
ggsave(p_plot, filename=paste0(path, "RTGP_simga_beta.png"), type="cairo-png")

## Trace plot (sigma_alpha)
df = data.frame('sigma_alpha' = result$sigma_alpha_chain[1:(iter)],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = sigma_alpha)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('sigma_alpha')
ggsave(p_plot, filename=paste0(path, "RTGP_simga_alpha.png"), type="cairo-png")

## Trace plot (beta0)
df = data.frame('beta0' = result$beta0_chain[1:(iter)],
                'Iteration' = 1:(iter))

p_plot = ggplot(df, aes(x = Iteration, y = beta0)) + 
  geom_hline(yintercept = beta0_true, color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('beta0')
ggsave(p_plot, filename=paste0(path, "RTGP_beta0.png"), type="cairo-png")


df = data.frame(beta_true = beta_true, beta_method = apply(result$beta_thresholded_chain[1:(iter),], 2, mean))
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta True') +
  ylab('Beta RTGP') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
ggsave(p_plot, filename=paste0(path, "scatterplot_beta_true_vs_beta_rtgp.png"), type="cairo-png")

}
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
 
}

# Save results and plot results with burn-in taken into account. 

time.end = Sys.time()
time = time.end - time.start
result_time$time_rtgp = time
result = list()
result$sigma_alpha_chain = sigma_alpha_chain
result$sigma_beta_chain = sigma_beta_chain
result$sigma_epsilon_chain = sigma_epsilon_chain
result$beta0_chain = beta0_chain
result$alpha_chain = alpha_chain
result$theta_chain = theta_chain
result$beta_chain = beta_chain
result$beta_thresholded_chain = beta_thresholded_chain
result$alpha_thresholded_chain = alpha_thresholded_chain
result$log_likelihood = log_likelihood_chain
result$log_lik_true = log_lik_true
result$joint = joint_chain
result$threshold_chain = threshold_chain
result$time = time
save(result, file = sprintf("%sresult_rtgp_gibbs.RData", path_results_RTGP_Gibbs ))

p_plot = plot_2D_funcs(apply(result$beta_chain[(burn_in + 1):(iter),],2,mean) , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_mean_beta_unthresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(apply(result$beta_thresholded_chain[(burn_in + 1):(iter),], 2, mean) , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_mean_beta_thresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(apply(result$beta_chain[(burn_in + 1):(iter),],2,sd)^2 , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_var_beta_unthresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(apply(result$beta_thresholded_chain[(burn_in + 1):(iter),],2,sd)^2 , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_var_beta_thresholded.png"), type="cairo-png")

p_plot = plot_2D_funcs(apply(alpha_thresholded_chain[(burn_in + 1):(iter),],2,mean) , v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_prob_alpha_greater_threshold.png"), type="cairo-png")

idx_beta = which(beta_true > 0, arr.ind = T)

df = data.frame('joint' = result$joint[(burn_in + 1):(iter)],
                'Iteration' = (burn_in + 1):(iter))
p_plot = ggplot(df, aes(x = Iteration, y = joint)) + 
  geom_line(aes()) +
  #ylim(-2.2, 2.2) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('Log-Joint')
ggsave(p_plot, filename=paste0(path, "RTGP_log_joint_with_true_joint.png"), type="cairo-png")

df = data.frame('log_likelihood' = result$log_likelihood[(burn_in + 1):(iter)],
                'Iteration' = (burn_in + 1):(iter))
p_plot = ggplot(df, aes(x = Iteration, y = log_likelihood)) + 
  geom_hline(yintercept = result$log_lik_true, color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('Log-Likelihood')
ggsave(p_plot, filename=paste0(path, "RTGP_log_likelihood_with_true_log_lik.png"), type="cairo-png")

## Trace plot (threshold)
df = data.frame('threshold' = result$threshold_chain[(burn_in + 1):(iter)],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = threshold)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('threshold')
ggsave(p_plot, filename=paste0(path, "RTGP_threshold.png"), type="cairo-png")

## Trace plots (unthresholded betas)
df = data.frame('Beta' = result$beta_chain[(burn_in + 1):(iter),idx_beta[1]],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = Beta)) + 
  geom_hline(yintercept = beta_true[idx_beta[1]], color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle(paste0('Beta (Voxel with Effect:', toString(idx_beta[1]) ,')'))
ggsave(p_plot, filename=paste0(path, "RTGP_beta_unthresholded_effect_voxel", toString(idx_beta[1]), ".png"), type="cairo-png")

df = data.frame('Beta' = result$beta_chain[(burn_in + 1):(iter),idx_beta[100]],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = Beta)) + 
  geom_hline(yintercept = beta_true[idx_beta[100]], color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle(paste0('Beta (Voxel with Effect:', toString(idx_beta[100]) ,')'))
ggsave(p_plot, filename=paste0(path, "RTGP_beta_unthresholded_effect_voxel", toString(idx_beta[100]), ".png"), type="cairo-png")

## Trace plots (alpha)
df = data.frame('Alpha' = result$alpha_chain[(burn_in + 1):(iter),idx_beta[1]],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = Alpha)) + 
  geom_hline(yintercept = alpha_true[idx_beta[1]], color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle(paste0('Alpha (Voxel with Effect:', toString(idx_beta[1]) ,')'))
ggsave(p_plot, filename=paste0(path, "RTGP_alpha_effect_voxel", toString(idx_beta[1]), ".png"), type="cairo-png")

df = data.frame('Alpha' = result$alpha_chain[(burn_in + 1):(iter),idx_beta[100]],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = Alpha)) + 
  geom_hline(yintercept = alpha_true[idx_beta[1]], color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle(paste0('Alpha (Voxel with Effect:', toString(idx_beta[100]) ,')'))
ggsave(p_plot, filename=paste0(path, "RTGP_alpha_effect_voxel", toString(idx_beta[100]), ".png"), type="cairo-png")

## Trace plot (sigma_epsilon)
df = data.frame('sigma_epsilon' = result$sigma_epsilon_chain[(burn_in + 1):(iter)],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = sigma_epsilon)) + 
  geom_hline(yintercept = sqrt(true_sigma2), color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('sigma_epsilon')
ggsave(p_plot, filename=paste0(path, "RTGP_simga_epsilon.png"), type="cairo-png")

## Trace plot (sigma_beta)
df = data.frame('sigma_beta' = result$sigma_beta_chain[(burn_in + 1):(iter)],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = sigma_beta)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('sigma_beta')
ggsave(p_plot, filename=paste0(path, "RTGP_simga_beta.png"), type="cairo-png")

## Trace plot (sigma_alpha)
df = data.frame('sigma_alpha' = result$sigma_alpha_chain[(burn_in + 1):(iter)],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = sigma_alpha)) + 
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('sigma_alpha')
ggsave(p_plot, filename=paste0(path, "RTGP_simga_alpha.png"), type="cairo-png")

## Trace plot (beta0)
df = data.frame('beta0' = result$beta0_chain[(burn_in + 1):(iter)],
                'Iteration' = (burn_in + 1):(iter))

p_plot = ggplot(df, aes(x = Iteration, y = beta0)) + 
  geom_hline(yintercept = beta0_true, color = 'red') +
  geom_line(aes()) +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank()) +
  ggtitle('beta0')
ggsave(p_plot, filename=paste0(path, "RTGP_beta0.png"), type="cairo-png")


df = data.frame(beta_true = beta_true, beta_method = apply(result$beta_thresholded_chain[(burn_in + 1):(iter),], 2, mean))
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta True') +
  ylab('Beta RTGP') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
ggsave(p_plot, filename=paste0(path, "scatterplot_beta_true_vs_beta_rtgp.png"), type="cairo-png")

# Evaluate bias of parameters
print('Bias: RTGP')
(bias_rtgp_no_effect = mean((apply(beta_thresholded_chain[(burn_in+1):(n_iter),], 2, mean)  - beta_true)[mask == 0]))
(bias_rtgp_effect = mean((apply(beta_thresholded_chain[(burn_in+1):(n_iter),], 2, mean)  - beta_true)[mask == 1]))
(bias_rtgp_total = mean((apply(beta_thresholded_chain[(burn_in+1):(n_iter),], 2, mean)  - beta_true)))

# Evaluate MSE of parameters
print('MSE: RTGP')
(mse_rtgp_no_effect = mean((apply(beta_thresholded_chain[(burn_in+1):(n_iter),], 2, mean)  - beta_true)[mask == 0]^2))
(mse_rtgp_effect = mean((apply(beta_thresholded_chain[(burn_in+1):(n_iter),], 2, mean)  - beta_true)[mask == 1]^2))
(mse_rtgp_total = mean((apply(beta_thresholded_chain[(burn_in+1):(n_iter),], 2, mean)  - beta_true)^2))

# Evaluate inference results
print('Inference results: RTGP')
t_bin_fdr = ifelse(apply(alpha_thresholded_chain[(burn_in + 1):(iter),],2, mean) > 0.5, 1, 0)
TP = t_bin_fdr + mask
TP = length(TP[TP == 2])
TN = t_bin_fdr + mask
TN = length(TN[TN == 0])
FP = t_bin_fdr - mask
FP = length(FP[FP == 1])
FN = t_bin_fdr - mask
FN = length(FN[FN == -1])
(TPR_rtgp = TP/(TP + FN))
(TDR_rtgp = TP/(TP + FP))
(FPR_rtgp = FP/(FP + TN))
(FDR_rtgp = FP/(FP + TP))

p_plot = plot_2D_funcs(t_bin_fdr, v_list)
ggsave(p_plot, filename=paste0(path, "RTGP_binary_significance.png"), type="cairo-png")

# Evaluate R2 (train) of prediction results
print('R2 Train: RTGP')
pred_rtgp = mean(beta0_chain[(burn_in + 1):(iter)]) + X%*%matrix(apply(beta_thresholded_chain[(burn_in + 1):(iter),], 2, mean) , length(beta_true), 1)
(R2_rtgp = cor(pred_rtgp,y)^2)

# Evaluate bias (train) of prediction results
print('Bias Prediction Train: RTGP')
(bias_y_pred_rtgp = mean(pred_rtgp - y))

# Evaluate MSE (train) of prediction results
print('MSE Prediction Train: RTGP')
(mse_y_pred_rtgp = mean((pred_rtgp - y)^2))

# Evaluate R2 (test) of prediction results
print('R2 Test: RTGP')
pred_rtgp_test = mean(beta0_chain[(burn_in + 1):(iter)]) + t(data$img_test)%*%matrix(apply(beta_thresholded_chain[(burn_in + 1):(iter),], 2, mean) , length(beta_true), 1)
(R2_rtgp_test = cor(pred_rtgp_test,data$y_test)^2)

# Evaluate bias (test) of prediction results
print('Bias Prediction Test: RTGP')
(bias_y_pred_rtgp_test = mean(pred_rtgp_test - data$y_test))

# Evaluate MSE (test) of prediction results
print('MSE Prediction Test: RTGP')
(mse_y_pred_rtgp_test = mean((pred_rtgp_test - data$y_test)^2))



##############
### Tables ###
##############

save(result_time, file = sprintf("%sresult_time.RData", path_results_all))

row_name = c('NM', 'STGP', 'RTGP-Gibbs', 'RTGP-VI')
column_name_inference = c('TPR', 'TDR', 'FPR', 'FDR')
column_name_parameters = c('Bias (no effect)', 'Bias (effect)', 'Bias (total)', 'MSE (no effect)', 'MSE (effect)', 'MSE (total)')
column_name_prediction = c('R2 (train)', 'Bias (train)', 'MSE (train)', 'R2 (test)', 'Bias (test)', 'MSE (test)')

# Inference results
matrix_inference = matrix(NA, 4, 4)
matrix_inference[1,1] = TPR_nm
matrix_inference[1,2] = TDR_nm
matrix_inference[1,3] = FPR_nm
matrix_inference[1,4] = FDR_nm
matrix_inference[2,1] = TPR_stgp
matrix_inference[2,2] = TDR_stgp
matrix_inference[2,3] = FPR_stgp
matrix_inference[2,4] = FDR_stgp
matrix_inference[3,1] = TPR_rtgp
matrix_inference[3,2] = TDR_rtgp
matrix_inference[3,3] = FPR_rtgp
matrix_inference[3,4] = FDR_rtgp
matrix_inference[4,1] = TPR_rtgp_vi
matrix_inference[4,2] = TDR_rtgp_vi
matrix_inference[4,3] = FPR_rtgp_vi
matrix_inference[4,4] = FDR_rtgp_vi
matrix_inference = data.frame(matrix_inference)
colnames(matrix_inference) = column_name_inference
rownames(matrix_inference) = row_name
xtable(matrix_inference, digits = 4)
print(xtable(matrix_inference, digits = 4), file=paste0(path_results_all, "table_inference.txt"))

# Parameter estimtes results
matrix_parameters = matrix(NA, 4, 6)
matrix_parameters[1,1] = bias_nm_no_effect
matrix_parameters[1,2] = bias_nm_effect
matrix_parameters[1,3] = bias_nm_total
matrix_parameters[1,4] = mse_nm_no_effect
matrix_parameters[1,5] = mse_nm_effect
matrix_parameters[1,6] = mse_nm_total

matrix_parameters[2,1] = bias_stgp_no_effect
matrix_parameters[2,2] = bias_stgp_effect
matrix_parameters[2,3] = bias_stgp_total
matrix_parameters[2,4] = mse_stgp_no_effect
matrix_parameters[2,5] = mse_stgp_effect
matrix_parameters[2,6] = mse_stgp_total

matrix_parameters[3,1] = bias_rtgp_no_effect
matrix_parameters[3,2] = bias_rtgp_effect
matrix_parameters[3,3] = bias_rtgp_total
matrix_parameters[3,4] = mse_rtgp_no_effect
matrix_parameters[3,5] = mse_rtgp_effect
matrix_parameters[3,6] = mse_rtgp_total

matrix_parameters[4,1] = bias_rtgp_vi_no_effect
matrix_parameters[4,2] = bias_rtgp_vi_effect
matrix_parameters[4,3] = bias_rtgp_vi_total
matrix_parameters[4,4] = mse_rtgp_vi_no_effect
matrix_parameters[4,5] = mse_rtgp_vi_effect
matrix_parameters[4,6] = mse_rtgp_vi_total
colnames(matrix_parameters) = column_name_parameters
rownames(matrix_parameters) = row_name
xtable(matrix_parameters, digits = 4)
print(xtable(matrix_parameters, digits = 4), file=paste0(path_results_all, "table_parameters.txt"))

# Prediction results
matrix_prediction = matrix(NA, 4, 6)
matrix_prediction[1,1] = R2_nm
matrix_prediction[1,2] = bias_y_pred_nm
matrix_prediction[1,3] = mse_y_pred_nm
matrix_prediction[1,4] = R2_nm_test
matrix_prediction[1,5] = bias_y_pred_nm_test
matrix_prediction[1,6] = mse_y_pred_nm_test

matrix_prediction[2,1] = R2_stgp
matrix_prediction[2,2] = bias_y_pred_stgp
matrix_prediction[2,3] = mse_y_pred_stgp
matrix_prediction[2,4] = R2_stgp_test
matrix_prediction[2,5] = bias_y_pred_stgp_test
matrix_prediction[2,6] = mse_y_pred_stgp_test

matrix_prediction[3,1] = R2_rtgp
matrix_prediction[3,2] = bias_y_pred_rtgp
matrix_prediction[3,3] = mse_y_pred_rtgp
matrix_prediction[3,4] = R2_rtgp_test
matrix_prediction[3,5] = bias_y_pred_rtgp_test
matrix_prediction[3,6] = mse_y_pred_rtgp_test

matrix_prediction[4,1] = R2_rtgp_vi
matrix_prediction[4,2] = bias_y_pred_rtgp_vi
matrix_prediction[4,3] = mse_y_pred_rtgp_vi
matrix_prediction[4,4] = R2_rtgp_vi_test
matrix_prediction[4,5] = bias_y_pred_rtgp_vi_test
matrix_prediction[4,6] = mse_y_pred_rtgp_vi_test
colnames(matrix_prediction) = column_name_prediction
rownames(matrix_prediction) = row_name
xtable(matrix_prediction, digits = 4)
print(xtable(matrix_prediction, digits = 4), file=paste0(path_results_all, "table_prediction.txt"))

result_matrix = list()
result_matrix$matrix_inference = matrix_inference
result_matrix$matrix_parameters = matrix_parameters
result_matrix$matrix_prediction = matrix_prediction
save(result_matrix, file = paste0(path_results_all, 'result_matrix.RData'))

