# Script: Perform GPR + Normal and GPR + Horseshoe model analysis for ABCD study data.

# Perform analysis for various different kernel hyperparameter settings.
sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Kernel hyperparameter settings (bandwith parameter)
phi_list = rep(c(0.01, 0.05, seq(0.1, 1, length.out = 10)), 5)
phi = phi_list[sim]
# Kernel hyperparameter settings (exponent parameter)
nu_list = c(rep(1.2,12), rep(1.4,12), rep(1.6,12), rep(1.8,12), rep(2,12))
nu = nu_list[sim]

#############
### Paths ###
#############

# Path to folder where Human Connectome workbench software is stored.
wb_path = ''

# Path to folder where GPR results will be stored in.
path_core = ''

# Path to folder where GPR results will be stored in (analysis with confounding variables)
path_results = paste0(path_core, 'L800_phi', toString(phi), '_nu', toString(nu), '/')
dir.create(file.path(path_results), showWarnings = FALSE)

# Path to folder where GPR results will be stored in (analysis without confounding variables)
path_results = paste0(path_core, 'L800_phi', toString(phi), '_nu', toString(nu), '/with_confounds/')
dir.create(file.path(path_results), showWarnings = FALSE)

# Path to folder where data is stored in.
path_data = ''

# Path to folder where eigendecompositons are stored in.
path_eigendecomp = ''
path_eigendecomp = paste0('eigendecomposition_L800_phi', toString(phi), '_nu', toString(nu), '/')

# Path to example cifti file (i.e. cope coefficients of a random subject in the ABCD data)
path_xii_file = ''

#################
### Libraries ###
#################

library(data.table)
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
library(bayestestR)

options(bitmapType='cairo-png')
library(ciftiTools)
ciftiTools.setOption('wb_path', wb_path)

################
### Analysis ###
################

######################
### with confounds ###
######################

### Left hemisphere ###

# Read in data matrix with confounding variables (training data).
X_confounds_train = fread(paste0(path_data, 'X_confounds_train.csv'))
X_confounds_train = data.matrix(X_confounds_train)[,2:ncol(X_confounds_train)]
# Read in data matrix with imaging data of left hemisphere (training data).
X_cortex_left_train = fread(paste0(path_data, 'X_left_cortex_train.csv'))
X_cortex_left_train = data.matrix(X_cortex_left_train)[,2:ncol(X_cortex_left_train)]
# Number of vertices (left hemisphere)
M = ncol(X_cortex_left_train)
# Output: Intelligence scores (training data)
y_train = as.numeric(unlist(fread(paste0(path_data, 'y_train.csv'))[,2]))

# Read in data matrix with confounding variables (validation data).
X_confounds_test = fread(paste0(path_data, "X_confounds_val.csv"))
X_confounds_test = X_confounds_test[,2:33]
X_confounds_test = data.matrix(X_confounds_test)

# Output: Intelligence scores (validation data)
y_test = read.csv("/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/not_centered/y_val.csv")[,2]

# Read in data matrix with imaging data of left hemisphere (validation data).
X_left_test = fread("/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/not_centered/X_left_cortex_val.csv")
X_left_test = X_left_test[,2:(M+1)]
X_left_test = data.matrix(X_left_test)

# Read in eigenfunctions and eigenvalues of eigendecomposition of the left hemisphere.
Psi = fread(paste0(path_eigendecomp, 'eigenvectors_left_cortex.csv'))
Psi = t(data.matrix(Psi)[,2:ncol(Psi)])
lambda_l = as.numeric(unlist(fread(paste0(path_eigendecomp, 'eigenvalues_left_cortex.csv'))[,2]))
# Square root the eigenvalues
sqrt_lambda = sqrt(lambda_l)
# Create bases through eigenfunctions and eigenvalues.
Bases = Psi*sqrt_lambda
# Number of basis functions
L = nrow(Bases)
# Transformed input data with intercept, confounding variables, and transformed input images (multiplied with bases)
Z = cbind(1, X_confounds_train, t(Bases%*%t(X_cortex_left_train)))
# Number of confounding variables 
P = ncol(X_confounds_train)

# Start time of analysis.
time.start = Sys.time()

# Perform model estimation of GPR + Normal, GPR + Horseshoe and GPR + MFVB.
hs_fit_SOI = fast_horseshoe_lm(y_train,Z, mcmc_sample = 5000, burnin = 2000)
nm_fit_SOI = fast_normal_lm(y_train,Z, mcmc_sample = 5000, burnin = 2000)
mfvb_fit_SOI = fast_mfvb_normal_lm(y_train,Z)

# End time of analysis.
time.end = Sys.time()
time = time.end - time.start

# Concatenate theta results of GPR + Normal, GPR + Horseshoe, and GPR + MFVB.
theta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-c(1:(P+1))],
                       NM = nm_fit_SOI$post_mean$betacoef[-c(1:(P+1))],
                       MFVB = mfvb_fit_SOI$post_mean$betacoef[-c(1:(P+1))])

# Concatenate beta results of GPR + Normal, GPR + Horseshoe, and GPR + MFVB.
beta_fit = data.frame(HS = crossprod(Bases,hs_fit_SOI$post_mean$betacoef[-c(1:(P+1))]),
                      NM = crossprod(Bases,nm_fit_SOI$post_mean$betacoef[-c(1:(P+1))]),
                      MFVB = crossprod(Bases,mfvb_fit_SOI$post_mean$betacoef[-c(1:(P+1))]))

# Inference results: GPR + Normal
# Determine binary significance through 95% equidistant credible interval (1: if interval does not contain 0, 0: if it does contain 0)
beta_NM = t(crossprod(Bases, nm_fit_SOI$mcmc$betacoef[(P+2):(L+P+1),]))
x=eti(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_NM = ifelse(x2==T,0,1)
# Determine binary significance through 95% highest posterior density credible interval (1: if interval does not contain 0, 0: if it does contain 0)
x=hdi(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_NM = ifelse(x2==T,0,1)
# Determine binary significance via test statistics (posterior mean / posterior standard deviation) and FDR correct p-values at 5%
t_beta_mean_NM = apply(beta_NM, 2, mean)
t_beta_sd_NM = apply(beta_NM, 2, sd)
t = t_beta_mean_NM / t_beta_sd_NM
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_NM = ifelse(p_fdr < 0.05, 1, 0)

# Inference results: GPR + Horseshoe
# Determine binary significance through 95% equidistant credible interval (1: if interval does not contain 0, 0: if it does contain 0)
beta_HS = t(crossprod(Bases, hs_fit_SOI$mcmc$betacoef[(P+2):(L+P+1),]))
x=eti(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_HS = ifelse(x2==T,0,1)
# Determine binary significance through 95% highest posterior density credible interval (1: if interval does not contain 0, 0: if it does contain 0)
x=hdi(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_HS = ifelse(x2==T,0,1)
# Determine binary significance via test statistics (posterior mean / posterior standard deviation) and FDR correct p-values at 5%
t_beta_mean_HS = apply(beta_HS, 2, mean)
t_beta_sd_HS = apply(beta_HS, 2, sd)
t = t_beta_mean_HS / t_beta_sd_HS
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_HS = ifelse(p_fdr < 0.05, 1, 0)

# Save results of left hemisphere in list.
result_left_cortex = list()
result_left_cortex$beta_HS_mcmc = beta_HS
result_left_cortex$beta_NM_mcmc = beta_NM
result_left_cortex$beta_NM_mean = t_beta_mean_NM
result_left_cortex$beta_HS_mean = t_beta_mean_HS
result_left_cortex$significance_test_statistic_NM = t_bin_fdr_NM
result_left_cortex$significance_test_statistic_HS = t_bin_fdr_HS
result_left_cortex$significance_CI_ETI_NM = CI_ETI_binary_NM
result_left_cortex$significance_CI_ETI_HS = CI_ETI_binary_HS
result_left_cortex$significance_CI_HPDI_NM = CI_HPDI_binary_NM
result_left_cortex$significance_CI_HPDI_HS = CI_HPDI_binary_HS
result_left_cortex$hs_fit_SOI = hs_fit_SOI
result_left_cortex$nm_fit_SOI = nm_fit_SOI
result_left_cortex$mfvb_fit_SOI = mfvb_fit_SOI
result_left_cortex$beta_fit = beta_fit
result_left_cortex$theta_fit = theta_fit
result_left_cortex$time = time

# Evaluate predictive results of GPR + Normal and GPR + Horseshoe models.
y_pred_train = list()
y_pred_test = list()
y_pred_train$NM = as.vector(X_confounds_train%*% matrix(result_left_cortex$nm_fit$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_cortex_left_train%*%matrix(result_left_cortex$beta_fit$NM,M,1)) + result_left_cortex$nm_fit$post_mean$betacoef[1,1]
y_pred_test$NM = as.vector(X_confounds_test%*% matrix(result_left_cortex$nm_fit$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_left_test%*%matrix(result_left_cortex$beta_fit$NM,M,1)) + result_left_cortex$nm_fit$post_mean$betacoef[1,1]

y_pred_train$HS = as.vector(X_confounds_train%*% matrix(result_left_cortex$hs_fit_SOI$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_cortex_left_train%*%matrix(result_left_cortex$beta_fit$HS,M,1)) + result_left_cortex$hs_fit_SOI$post_mean$betacoef[1,1]
y_pred_test$HS = as.vector(X_confounds_test%*% matrix(result_left_cortex$hs_fit_SOI$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_left_test%*%matrix(result_left_cortex$beta_fit$HS,M,1)) + result_left_cortex$hs_fit_SOI$post_mean$betacoef[1,1]

y_pred_train$MFVB = as.vector(X_confounds_train%*% matrix(result_left_cortex$mfvb_fit$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_cortex_left_train%*%matrix(result_left_cortex$beta_fit$MFVB,M,1)) + result_left_cortex$mfvb_fit$post_mean$betacoef[1,1]
y_pred_test$MFVB = as.vector(X_confounds_test%*% matrix(result_left_cortex$mfvb_fit$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_left_test%*%matrix(result_left_cortex$beta_fit$MFVB,M,1)) + result_left_cortex$mfvb_fit$post_mean$betacoef[1,1]

result_R2_MSE = list()
result_R2_MSE$R2_train_NM = cor(y_pred_train$NM, y_train)^2
result_R2_MSE$R2_test_NM = cor(y_pred_test$NM, y_test)^2

result_R2_MSE$mse_y_pred_train_NM = mean((y_pred_train$NM - y_train)^2)
result_R2_MSE$mse_y_pred_test_NM = mean((y_pred_test$NM - y_test)^2)

result_R2_MSE$R2_train_HS = cor(y_pred_train$HS, y_train)^2
result_R2_MSE$R2_test_HS = cor(y_pred_test$HS, y_test)^2

result_R2_MSE$mse_y_pred_train_HS = mean((y_pred_train$HS - y_train)^2)
result_R2_MSE$mse_y_pred_test_HS = mean((y_pred_test$HS - y_test)^2)

result_R2_MSE$R2_train_MFVB = cor(y_pred_train$MFVB, y_train)^2
result_R2_MSE$R2_test_MFVB = cor(y_pred_test$MFVB, y_test)^2

result_R2_MSE$mse_y_pred_train_MFVB = mean((y_pred_train$MFVB - y_train)^2)
result_R2_MSE$mse_y_pred_test_MFVB = mean((y_pred_test$MFVB - y_test)^2)

result_left_cortex$y_pred_train = y_pred_train
result_left_cortex$y_pred_test = y_pred_test
result_left_cortex$result_R2_MSE = result_R2_MSE

# Save results of predictive evaluation in a list and table.
result_eval = matrix(NA, 3, 4)
result_eval[1,1] = result_R2_MSE$R2_train_NM 
result_eval[1,2] = result_R2_MSE$mse_y_pred_train_NM 
result_eval[1,3] = result_R2_MSE$R2_test_NM 
result_eval[1,4] = result_R2_MSE$mse_y_pred_test_NM
result_eval[2,1] = result_R2_MSE$R2_train_HS
result_eval[2,2] = result_R2_MSE$mse_y_pred_train_HS
result_eval[2,3] = result_R2_MSE$R2_test_HS
result_eval[2,4] = result_R2_MSE$mse_y_pred_test_HS
result_eval[3,1] = result_R2_MSE$R2_train_MFVB
result_eval[3,2] = result_R2_MSE$mse_y_pred_train_MFVB
result_eval[3,3] = result_R2_MSE$R2_test_MFVB
result_eval[3,4] = result_R2_MSE$mse_y_pred_test_MFVB
col_names = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')
row_names = c('NM', 'HS', 'MFVB')
colnames(result_eval) = col_names
rownames(result_eval) = row_names

xtable(result_eval, digits = 4)
print(xtable(result_eval, digits = 4), file=paste0(path_results, "table_result_left_cortex.txt"))
write.csv(result_eval, paste0(path_results, 'result_matrix.csv'))

save(result_left_cortex, file = sprintf("%sresult_left_cortex.RData", path_results))

### Right hemisphere ###

# Read in data matrix with imaging data of right hemisphere (training data).
X_cortex_right_train = fread(paste0(path_data, 'X_right_cortex_train.csv'))
X_cortex_right_train = data.matrix(X_cortex_right_train)[,2:ncol(X_cortex_right_train)]
# Number of vertices
M = ncol(X_cortex_right_train)

# Read in data matrix with imaging data of right hemisphere (validation data).
X_right_test = fread(paste0(path_data, "X_right_cortex_val.csv"))
X_right_test = X_right_test[,2:dim(X_right_test)[2]]
X_right_test = data.matrix(X_right_test)

# Read in eigenvalues and eigenfunctions of eigendecomposition of right hemisphere.
Psi = fread(paste0(path_eigendecomp, 'eigenvectors_right_cortex.csv'))
Psi = t(data.matrix(Psi)[,2:ncol(Psi)])
lambda_l = as.numeric(unlist(fread(paste0(path_eigendecomp, 'eigenvalues_right_cortex.csv'))[,2]))
# Square root eigenvalues.
sqrt_lambda = sqrt(lambda_l)
# Create bases through eigenvalues and eigenfunctions.
Bases = Psi*sqrt_lambda
# Number of basis functions
L = nrow(Bases)
# Transformed input data with intercept, confounding variables, and transformed input images (multiplied with bases)
Z = cbind(1, X_confounds_train, t(Bases%*%t(X_cortex_right_train)))
# Number of confounding variables
P = ncol(X_confounds_train)

# Start time of analysis.
time.start = Sys.time()

# Perform model estimation of GPR + Normal, GPR + Horseshoe, and GPR + MFVB.
hs_fit_SOI = fast_horseshoe_lm(y_train,Z, mcmc_sample = 5000, burnin = 2000)
nm_fit_SOI = fast_normal_lm(y_train,Z, mcmc_sample = 5000, burnin = 2000)
mfvb_fit_SOI = fast_mfvb_normal_lm(y_train,Z)

# End time of analysis.
time.end = Sys.time()
time = time.end - time.start

# Concatenate theta results of GPR + Normal, GPR + Horseshoe, and GPR + MFVB.
theta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-c(1:(P+1))],
                       NM = nm_fit_SOI$post_mean$betacoef[-c(1:(P+1))],
                       MFVB = mfvb_fit_SOI$post_mean$betacoef[-c(1:(P+1))])

# Concatenate beta results of GPR + Normal, GPR + Horseshoe, and GPR + MFVB.
beta_fit = data.frame(HS = crossprod(Bases,hs_fit_SOI$post_mean$betacoef[-c(1:(P+1))]),
                      NM = crossprod(Bases,nm_fit_SOI$post_mean$betacoef[-c(1:(P+1))]),
                      MFVB = crossprod(Bases,mfvb_fit_SOI$post_mean$betacoef[-c(1:(P+1))]))

# Inference results: GPR + Normal
# Determine binary significance through 95% equidistant credible interval (1: if interval does not contain 0, 0: if it does contain 0)
beta_NM = t(crossprod(Bases, nm_fit_SOI$mcmc$betacoef[(P+2):(L+P+1),]))
x=eti(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_NM = ifelse(x2==T,0,1)
# Determine binary significance through 95% highest posterior density credible interval (1: if interval does not contain 0, 0: if it does contain 0)
x=hdi(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_NM = ifelse(x2==T,0,1)
# Determine binary significance via test statistics (posterior mean / posterior standard deviation) and FDR correct p-values at 5%
t_beta_mean_NM = apply(beta_NM, 2, mean)
t_beta_sd_NM = apply(beta_NM, 2, sd)
t = t_beta_mean_NM / t_beta_sd_NM
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_NM = ifelse(p_fdr < 0.05, 1, 0)

# Inference results: GPR + Horseshoe
# Determine binary significance through 95% equidistant credible interval (1: if interval does not contain 0, 0: if it does contain 0)
beta_HS = t(crossprod(Bases, hs_fit_SOI$mcmc$betacoef[(P+2):(L+P+1),]))
x=eti(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_HS = ifelse(x2==T,0,1)
# Determine binary significance through 95% highest posterior density credible interval (1: if interval does not contain 0, 0: if it does contain 0)
x=hdi(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_HS = ifelse(x2==T,0,1)
# Determine binary significance via test statistics (posterior mean / posterior standard deviation) and FDR correct p-values at 5%
t_beta_mean_HS = apply(beta_HS, 2, mean)
t_beta_sd_HS = apply(beta_HS, 2, sd)
t = t_beta_mean_HS / t_beta_sd_HS
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_HS = ifelse(p_fdr < 0.05, 1, 0)

# Save results of right hemisphere in list.
result_right_cortex = list()
result_right_cortex$beta_HS_mcmc = beta_HS
result_right_cortex$beta_NM_mcmc = beta_NM
result_right_cortex$beta_NM_mean = t_beta_mean_NM
result_right_cortex$beta_HS_mean = t_beta_mean_HS
result_right_cortex$significance_test_statistic_NM = t_bin_fdr_NM
result_right_cortex$significance_test_statistic_HS = t_bin_fdr_HS
result_right_cortex$significance_CI_ETI_NM = CI_ETI_binary_NM
result_right_cortex$significance_CI_ETI_HS = CI_ETI_binary_HS
result_right_cortex$significance_CI_HPDI_NM = CI_HPDI_binary_NM
result_right_cortex$significance_CI_HPDI_HS = CI_HPDI_binary_HS
result_right_cortex$hs_fit_SOI = hs_fit_SOI
result_right_cortex$nm_fit_SOI = nm_fit_SOI
result_right_cortex$mfvb_fit_SOI = mfvb_fit_SOI
result_right_cortex$beta_fit = beta_fit
result_right_cortex$theta_fit = theta_fit
result_right_cortex$time = time

# Evaluate predictive performance of GPR + Normal and GPR + Horseshoe.
y_pred_train = list()
y_pred_test = list()

y_pred_train$NM = as.vector(X_confounds_train%*% matrix(result_right_cortex$nm_fit$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_cortex_right_train%*%matrix(result_right_cortex$beta_fit$NM,M,1)) + result_right_cortex$nm_fit$post_mean$betacoef[1,1]
y_pred_test$NM = as.vector(X_confounds_test%*% matrix(result_right_cortex$nm_fit$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_right_test%*%matrix(result_right_cortex$beta_fit$NM,M,1)) + result_right_cortex$nm_fit$post_mean$betacoef[1,1]
y_pred_train$HS = as.vector(X_confounds_train%*% matrix(result_right_cortex$hs_fit_SOI$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_cortex_right_train%*%matrix(result_right_cortex$beta_fit$HS,M,1)) + result_right_cortex$hs_fit_SOI$post_mean$betacoef[1,1]
y_pred_test$HS = as.vector(X_confounds_test%*% matrix(result_right_cortex$hs_fit_SOI$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_right_test%*%matrix(result_right_cortex$beta_fit$HS,M,1)) + result_right_cortex$hs_fit_SOI$post_mean$betacoef[1,1]
y_pred_train$MFVB = as.vector(X_confounds_train%*% matrix(result_right_cortex$mfvb_fit$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_cortex_right_train%*%matrix(result_right_cortex$beta_fit$MFVB,M,1)) + result_right_cortex$mfvb_fit$post_mean$betacoef[1,1]
y_pred_test$MFVB = as.vector(X_confounds_test%*% matrix(result_right_cortex$mfvb_fit$post_mean$betacoef[2:33,1], 32, 1))  + as.vector(X_right_test%*%matrix(result_right_cortex$beta_fit$MFVB,M,1)) + result_right_cortex$mfvb_fit$post_mean$betacoef[1,1]

result_R2_MSE = list()
result_R2_MSE$R2_train_NM = cor(y_pred_train$NM, y_train)^2
result_R2_MSE$R2_test_NM = cor(y_pred_test$NM, y_test)^2

result_R2_MSE$mse_y_pred_train_NM = mean((y_pred_train$NM - y_train)^2)
result_R2_MSE$mse_y_pred_test_NM = mean((y_pred_test$NM - y_test)^2)

result_R2_MSE$R2_train_HS = cor(y_pred_train$HS, y_train)^2
result_R2_MSE$R2_test_HS = cor(y_pred_test$HS, y_test)^2

result_R2_MSE$mse_y_pred_train_HS = mean((y_pred_train$HS - y_train)^2)
result_R2_MSE$mse_y_pred_test_HS = mean((y_pred_test$HS - y_test)^2)

result_R2_MSE$R2_train_MFVB = cor(y_pred_train$MFVB, y_train)^2
result_R2_MSE$R2_test_MFVB = cor(y_pred_test$MFVB, y_test)^2

result_R2_MSE$mse_y_pred_train_MFVB = mean((y_pred_train$MFVB - y_train)^2)
result_R2_MSE$mse_y_pred_test_MFVB = mean((y_pred_test$MFVB - y_test)^2)

result_right_cortex$y_pred_train = y_pred_train
result_right_cortex$y_pred_test = y_pred_test
result_right_cortex$result_R2_MSE = result_R2_MSE

# Save predictive results and store in list and tables.
result_eval = matrix(NA, 3, 4)
result_eval[1,1] = result_R2_MSE$R2_train_NM 
result_eval[1,2] = result_R2_MSE$mse_y_pred_train_NM 
result_eval[1,3] = result_R2_MSE$R2_test_NM 
result_eval[1,4] = result_R2_MSE$mse_y_pred_test_NM
result_eval[2,1] = result_R2_MSE$R2_train_HS
result_eval[2,2] = result_R2_MSE$mse_y_pred_train_HS
result_eval[2,3] = result_R2_MSE$R2_test_HS
result_eval[2,4] = result_R2_MSE$mse_y_pred_test_HS
result_eval[3,1] = result_R2_MSE$R2_train_MFVB
result_eval[3,2] = result_R2_MSE$mse_y_pred_train_MFVB
result_eval[3,3] = result_R2_MSE$R2_test_MFVB
result_eval[3,4] = result_R2_MSE$mse_y_pred_test_MFVB
col_names = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')
row_names = c('NM', 'HS', 'MFVB')
colnames(result_eval) = col_names
rownames(result_eval) = row_names

xtable(result_eval, digits = 4)
print(xtable(result_eval, digits = 4), file=paste0(path_results, "table_result_right_cortex.txt"))
write.csv(result_eval, paste0(path_results, 'result_matrix_right.csv'))

save(result_right_cortex, file = sprintf("%sresult_right_cortex.RData", path_results))


### Visualisation ###

# Read in results (left and right cortex)
load(paste0(path_results, 'result_left_cortex.RData'))
load(paste0(path_results, 'result_right_cortex.RData'))

# Save beta maps of GPR + Normal model. 
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$beta_fit$NM, length(result_left_cortex$beta_fit$NM), 1)
xii$data$cortex_right = matrix(result_right_cortex$beta_fit$NM, length(result_right_cortex$beta_fit$NM), 1)

write_xifti(
  xii, 
  file.path(path_results, "beta_NM.dtseries.nii"), 
  file.path(path_results, "beta_NM_cortex_left.surf.gii"), file.path(path_results, "beta_NM_cortex_right.surf.gii")
)

# Save beta maps of GPR + Horseshoe model. 
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$beta_fit$HS, length(result_left_cortex$beta_fit$HS), 1)
xii$data$cortex_right = matrix(result_right_cortex$beta_fit$HS, length(result_right_cortex$beta_fit$HS), 1)

write_xifti(
  xii, 
  file.path(path_results, "beta_HS.dtseries.nii"), 
  file.path(path_results, "beta_HS_cortex_left.surf.gii"), file.path(path_results, "beta_HS_cortex_right.surf.gii")
)

# Save binary significance maps (test statistic-based evaluation) of GPR + Normal model. 
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_test_statistic_NM, length(result_left_cortex$significance_test_statistic_NM), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_test_statistic_NM, length(result_right_cortex$significance_test_statistic_NM), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_test_statistic_NM.dtseries.nii"), 
  file.path(path_results, "significance_test_statistic_NM_cortex_left.surf.gii"), file.path(path_results, "significance_test_statistic_NM_cortex_right.surf.gii")
)

# Save binary significance maps (test statistic-based evaluation) of GPR + Horseshoe model. 
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_test_statistic_HS, length(result_left_cortex$significance_test_statistic_HS), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_test_statistic_HS, length(result_right_cortex$significance_test_statistic_HS), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_test_statistic_HS.dtseries.nii"), 
  file.path(path_results, "significance_test_statistic_HS_cortex_left.surf.gii"), file.path(path_results, "significance_test_statistic_HS_cortex_right.surf.gii")
)

# Save binary significance maps (ETI-based evaluation) of GPR + Normal model. 
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_CI_ETI_NM, length(result_left_cortex$significance_CI_ETI_NM), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_CI_ETI_NM, length(result_right_cortex$significance_CI_ETI_NM), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_CI_ETI_NM.dtseries.nii"), 
  file.path(path_results, "significance_CI_ETI_NM_cortex_left.surf.gii"), file.path(path_results, "significance_CI_ETI_NM_cortex_right.surf.gii")
)

# Save binary significance maps (ETI-based evaluation) of GPR + Horseshoe model. 
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_CI_ETI_HS, length(result_left_cortex$significance_CI_ETI_HS), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_CI_ETI_HS, length(result_right_cortex$significance_CI_ETI_HS), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_CI_ETI_HS.dtseries.nii"), 
  file.path(path_results, "significance_CI_ETI_HS_cortex_left.surf.gii"), file.path(path_results, "significance_CI_ETI_HS_cortex_right.surf.gii")
)

# Save binary significance maps (HPDI-based evaluation) of GPR + Normal model. 
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_CI_HPDI_NM, length(result_left_cortex$significance_CI_HPDI_NM), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_CI_HPDI_NM, length(result_right_cortex$significance_CI_HPDI_NM), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_CI_HPDI_NM.dtseries.nii"), 
  file.path(path_results, "significance_CI_HPDI_NM_cortex_left.surf.gii"), file.path(path_results, "significance_CI_HPDI_NM_cortex_right.surf.gii")
)

# Save binary significance maps (HPDI-based evaluation) of GPR + Horseshoe model. 
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_CI_HPDI_HS, length(result_left_cortex$significance_CI_HPDI_HS), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_CI_HPDI_HS, length(result_right_cortex$significance_CI_HPDI_HS), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_CI_HPDI_HS.dtseries.nii"), 
  file.path(path_results, "significance_CI_HPDI_HS_cortex_left.surf.gii"), file.path(path_results, "significance_CI_HPDI_HS_cortex_right.surf.gii")
)


####################
### no confounds ###
####################

# Set current path to store results in for the analysis without confounding variables.
path_results = paste0(path_core, 'L800_phi', toString(phi), '_nu', toString(nu), '/no_confounds/')
dir.create(file.path(path_results), showWarnings = FALSE)

### Left hemisphere ###

# Read in image data of left hemisphere (training data).
X_cortex_left_train = fread(paste0(path_data, 'X_left_cortex_train.csv'))
X_cortex_left_train = data.matrix(X_cortex_left_train)[,2:ncol(X_cortex_left_train)]
# Number of vertices
M = ncol(X_cortex_left_train)
# Read in output vector of intelligence scores (training data).
y_train = as.numeric(unlist(fread(paste0(path_data, 'y_train.csv'))[,2]))
# Read in output vector of intelligence scores (validation data).
y_test = read.csv(paste0(path_data, "y_val.csv"))[,2]

# Read in image data of left hemisphere (validation data).
X_left_test = fread(paste0(path_data, "X_left_cortex_val.csv"))
X_left_test = X_left_test[,2:(M+1)]
X_left_test = data.matrix(X_left_test)

# Read in eigenvalues and eigenfunctions of eigendecomposition of left hemisphere.
Psi = fread(paste0(path_eigendecomp, 'eigenvectors_left_cortex.csv'))
Psi = t(data.matrix(Psi)[,2:ncol(Psi)])
lambda_l = as.numeric(unlist(fread(paste0(path_eigendecomp, 'eigenvalues_left_cortex.csv'))[,2]))
# Square root eigenvalues.
sqrt_lambda = sqrt(lambda_l)
# Create bases through eigenvalues and eigenfunctions.
Bases = Psi*sqrt_lambda
# Number of basis functions.
L = nrow(Bases)
# Transformed input containing intercept and transformed input image (multiplied with bases).
Z = cbind(1, t(Bases%*%t(X_cortex_left_train)))

# Start time of the analysis.
time.start = Sys.time()

# Perform model estimation of GPR + Normal, GPR + Horseshoe, and GPR + MFVB.
hs_fit_SOI = fast_horseshoe_lm(y_train,Z, mcmc_sample = 5000, burnin = 2000)
nm_fit_SOI = fast_normal_lm(y_train,Z, mcmc_sample = 5000, burnin = 2000)
mfvb_fit_SOI = fast_mfvb_normal_lm(y_train,Z)

# End time of the analysis.
time.end = Sys.time()
time = time.end - time.start

# Concatenate theta results of GPR + Normal, GPR + Horseshoe, and GPR + MFVB into dataframe.
theta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-1],
                       NM = nm_fit_SOI$post_mean$betacoef[-1],
                       MFVB = mfvb_fit_SOI$post_mean$betacoef[-1])

# Concatenate beta results of GPR + Normal, GPR + Horseshoe, and GPR + MFVB into dataframe.
beta_fit = data.frame(HS = crossprod(Bases,hs_fit_SOI$post_mean$betacoef[-1]),
                      NM = crossprod(Bases,nm_fit_SOI$post_mean$betacoef[-1]),
                      MFVB = crossprod(Bases,mfvb_fit_SOI$post_mean$betacoef[-1]))


# Inference results: GPR + Normal
# Determine binary significance through 95% equidistant credible interval (1: if interval does not contain 0, 0: if it does contain 0)
beta_NM = t(crossprod(Bases, nm_fit_SOI$mcmc$betacoef[(P+2):(L+P+1),]))
x=eti(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_NM = ifelse(x2==T,0,1)
# Determine binary significance through 95% highest posterior density credible interval (1: if interval does not contain 0, 0: if it does contain 0)
x=hdi(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_NM = ifelse(x2==T,0,1)
# Determine binary significance via test statistics (posterior mean / posterior standard deviation) and FDR correct p-values at 5%
t_beta_mean_NM = apply(beta_NM, 2, mean)
t_beta_sd_NM = apply(beta_NM, 2, sd)
t = t_beta_mean_NM / t_beta_sd_NM
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_NM = ifelse(p_fdr < 0.05, 1, 0)

# Inference results: GPR + Horseshoe
# Determine binary significance through 95% equidistant credible interval (1: if interval does not contain 0, 0: if it does contain 0)
beta_HS = t(crossprod(Bases, hs_fit_SOI$mcmc$betacoef[(P+2):(L+P+1),]))
x=eti(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_HS = ifelse(x2==T,0,1)
# Determine binary significance through 95% highest posterior density credible interval (1: if interval does not contain 0, 0: if it does contain 0)
x=hdi(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_HS = ifelse(x2==T,0,1)
# Determine binary significance via test statistics (posterior mean / posterior standard deviation) and FDR correct p-values at 5%
t_beta_mean_HS = apply(beta_HS, 2, mean)
t_beta_sd_HS = apply(beta_HS, 2, sd)
t = t_beta_mean_HS / t_beta_sd_HS
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_HS = ifelse(p_fdr < 0.05, 1, 0)

# Save results of left hemisphere in list.
result_left_cortex = list()
result_left_cortex$beta_HS_mcmc = beta_HS
result_left_cortex$beta_NM_mcmc = beta_NM
result_left_cortex$beta_NM_mean = t_beta_mean_NM
result_left_cortex$beta_HS_mean = t_beta_mean_HS
result_left_cortex$significance_test_statistic_NM = t_bin_fdr_NM
result_left_cortex$significance_test_statistic_HS = t_bin_fdr_HS
result_left_cortex$significance_CI_ETI_NM = CI_ETI_binary_NM
result_left_cortex$significance_CI_ETI_HS = CI_ETI_binary_HS
result_left_cortex$significance_CI_HPDI_NM = CI_HPDI_binary_NM
result_left_cortex$significance_CI_HPDI_HS = CI_HPDI_binary_HS
result_left_cortex$hs_fit_SOI = hs_fit_SOI
result_left_cortex$nm_fit_SOI = nm_fit_SOI
result_left_cortex$mfvb_fit_SOI = mfvb_fit_SOI
result_left_cortex$beta_fit = beta_fit
result_left_cortex$theta_fit = theta_fit
result_left_cortex$time = time

# Evaluate predictive results.
y_pred_train = list()
y_pred_test = list()

y_pred_train$NM =as.vector(X_cortex_left_train%*%matrix(result_left_cortex$beta_fit$NM,M,1)) + result_left_cortex$nm_fit$post_mean$betacoef[1,1]
y_pred_test$NM =  as.vector(X_left_test%*%matrix(result_left_cortex$beta_fit$NM,M,1)) + result_left_cortex$nm_fit$post_mean$betacoef[1,1]
y_pred_train$HS = as.vector(X_cortex_left_train%*%matrix(result_left_cortex$beta_fit$HS,M,1)) + result_left_cortex$hs_fit_SOI$post_mean$betacoef[1,1]
y_pred_test$HS =  as.vector(X_left_test%*%matrix(result_left_cortex$beta_fit$HS,M,1)) + result_left_cortex$hs_fit_SOI$post_mean$betacoef[1,1]
y_pred_train$MFVB = as.vector(X_cortex_left_train%*%matrix(result_left_cortex$beta_fit$MFVB,M,1)) + result_left_cortex$mfvb_fit$post_mean$betacoef[1,1]
y_pred_test$MFVB = as.vector(X_left_test%*%matrix(result_left_cortex$beta_fit$MFVB,M,1)) + result_left_cortex$mfvb_fit$post_mean$betacoef[1,1]

result_R2_MSE = list()
result_R2_MSE$R2_train_NM = cor(y_pred_train$NM, y_train)^2
result_R2_MSE$R2_test_NM = cor(y_pred_test$NM, y_test)^2

result_R2_MSE$mse_y_pred_train_NM = mean((y_pred_train$NM - y_train)^2)
result_R2_MSE$mse_y_pred_test_NM = mean((y_pred_test$NM - y_test)^2)

result_R2_MSE$R2_train_HS = cor(y_pred_train$HS, y_train)^2
result_R2_MSE$R2_test_HS = cor(y_pred_test$HS, y_test)^2

result_R2_MSE$mse_y_pred_train_HS = mean((y_pred_train$HS - y_train)^2)
result_R2_MSE$mse_y_pred_test_HS = mean((y_pred_test$HS - y_test)^2)

result_R2_MSE$R2_train_MFVB = cor(y_pred_train$MFVB, y_train)^2
result_R2_MSE$R2_test_MFVB = cor(y_pred_test$MFVB, y_test)^2

result_R2_MSE$mse_y_pred_train_MFVB = mean((y_pred_train$MFVB - y_train)^2)
result_R2_MSE$mse_y_pred_test_MFVB = mean((y_pred_test$MFVB - y_test)^2)

result_left_cortex$y_pred_train = y_pred_train
result_left_cortex$y_pred_test = y_pred_test
result_left_cortex$result_R2_MSE = result_R2_MSE

# Save predictive results in list and table.
result_eval = matrix(NA, 3, 4)
result_eval[1,1] = result_R2_MSE$R2_train_NM 
result_eval[1,2] = result_R2_MSE$mse_y_pred_train_NM 
result_eval[1,3] = result_R2_MSE$R2_test_NM 
result_eval[1,4] = result_R2_MSE$mse_y_pred_test_NM
result_eval[2,1] = result_R2_MSE$R2_train_HS
result_eval[2,2] = result_R2_MSE$mse_y_pred_train_HS
result_eval[2,3] = result_R2_MSE$R2_test_HS
result_eval[2,4] = result_R2_MSE$mse_y_pred_test_HS
result_eval[3,1] = result_R2_MSE$R2_train_MFVB
result_eval[3,2] = result_R2_MSE$mse_y_pred_train_MFVB
result_eval[3,3] = result_R2_MSE$R2_test_MFVB
result_eval[3,4] = result_R2_MSE$mse_y_pred_test_MFVB
col_names = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')
row_names = c('NM', 'HS', 'MFVB')
colnames(result_eval) = col_names
rownames(result_eval) = row_names

xtable(result_eval, digits = 4)
print(xtable(result_eval, digits = 4), file=paste0(path_results, "table_result_left_cortex.txt"))
write.csv(result_eval, paste0(path_results, 'result_matrix.csv'))

save(result_left_cortex, file = sprintf("%sresult_left_cortex.RData", path_results))

### Right cortex ###

# Read in imaging data of right hemisphere (training data).
X_cortex_right_train = fread(paste0(path_data, 'X_right_cortex_train.csv'))
X_cortex_right_train = data.matrix(X_cortex_right_train)[,2:ncol(X_cortex_right_train)]
# Number of vertices
M = ncol(X_cortex_right_train)

# Read in imaging data of right hemisphere (validation data)
X_right_test = fread(paste0(path_data, "X_right_cortex_val.csv"))
X_right_test = X_right_test[,2:dim(X_right_test)[2]]
X_right_test = data.matrix(X_right_test)

# Read in eigenvalues and eigenfunctions of eigendecomposition of right hemisphere.
Psi = fread(paste0(path_eigendecomp, 'eigenvectors_right_cortex.csv'))
Psi = t(data.matrix(Psi)[,2:ncol(Psi)])
lambda_l = as.numeric(unlist(fread(paste0(path_eigendecomp, 'eigenvalues_right_cortex.csv'))[,2]))
# Square root eigenvalues.
sqrt_lambda = sqrt(lambda_l)
# Create bases through eigenvalues and eigenfunctions.
Bases = Psi*sqrt_lambda
# Number of basis functions
L = nrow(Bases)
# Transform input containing intercept and transformed input images (multiplied by bases)
Z = cbind(1, t(Bases%*%t(X_cortex_right_train)))

# Start time of analysis.
time.start = Sys.time()

# Perform model estimation of GPR + Normal, GPR + Horseshoe, and GPR + MFVB
hs_fit_SOI = fast_horseshoe_lm(y_train,Z, mcmc_sample = 5000, burnin = 2000)
nm_fit_SOI = fast_normal_lm(y_train,Z, mcmc_sample = 5000, burnin = 2000)
mfvb_fit_SOI = fast_mfvb_normal_lm(y_train,Z)

# End time of analysis
time.end = Sys.time()
time = time.end - time.start

# Concatenate theta results in dataframe.
theta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-1],
                       NM = nm_fit_SOI$post_mean$betacoef[-1],
                       MFVB = mfvb_fit_SOI$post_mean$betacoef[-1])

# Concatenate beta results in dataframe.
beta_fit = data.frame(HS = crossprod(Bases,hs_fit_SOI$post_mean$betacoef[-1]),
                      NM = crossprod(Bases,nm_fit_SOI$post_mean$betacoef[-1]),
                      MFVB = crossprod(Bases,mfvb_fit_SOI$post_mean$betacoef[-1]))

# Inference results: GPR + Normal
# Determine binary significance through 95% equidistant credible interval (1: if interval does not contain 0, 0: if it does contain 0)
beta_NM = t(crossprod(Bases, nm_fit_SOI$mcmc$betacoef[(P+2):(L+P+1),]))
x=eti(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_NM = ifelse(x2==T,0,1)
# Determine binary significance through 95% highest posterior density credible interval (1: if interval does not contain 0, 0: if it does contain 0)
x=hdi(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_NM = ifelse(x2==T,0,1)
# Determine binary significance via test statistics (posterior mean / posterior standard deviation) and FDR correct p-values at 5%
t_beta_mean_NM = apply(beta_NM, 2, mean)
t_beta_sd_NM = apply(beta_NM, 2, sd)
t = t_beta_mean_NM / t_beta_sd_NM
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_NM = ifelse(p_fdr < 0.05, 1, 0)

# Inference results: GPR + Horseshoe
# Determine binary significance through 95% equidistant credible interval (1: if interval does not contain 0, 0: if it does contain 0)
beta_HS = t(crossprod(Bases, hs_fit_SOI$mcmc$betacoef[(P+2):(L+P+1),]))
x=eti(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_HS = ifelse(x2==T,0,1)
# Determine binary significance through 95% highest posterior density credible interval (1: if interval does not contain 0, 0: if it does contain 0)
x=hdi(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_HS = ifelse(x2==T,0,1)
# Determine binary significance via test statistics (posterior mean / posterior standard deviation) and FDR correct p-values at 5%
t_beta_mean_HS = apply(beta_HS, 2, mean)
t_beta_sd_HS = apply(beta_HS, 2, sd)
t = t_beta_mean_HS / t_beta_sd_HS
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_HS = ifelse(p_fdr < 0.05, 1, 0)

# Save results of right hemisphere in list.
result_right_cortex = list()
result_right_cortex$beta_HS_mcmc = beta_HS
result_right_cortex$beta_NM_mcmc = beta_NM
result_right_cortex$beta_NM_mean = t_beta_mean_NM
result_right_cortex$beta_HS_mean = t_beta_mean_HS
result_right_cortex$significance_test_statistic_NM = t_bin_fdr_NM
result_right_cortex$significance_test_statistic_HS = t_bin_fdr_HS
result_right_cortex$significance_CI_ETI_NM = CI_ETI_binary_NM
result_right_cortex$significance_CI_ETI_HS = CI_ETI_binary_HS
result_right_cortex$significance_CI_HPDI_NM = CI_HPDI_binary_NM
result_right_cortex$significance_CI_HPDI_HS = CI_HPDI_binary_HS
result_right_cortex$hs_fit_SOI = hs_fit_SOI
result_right_cortex$nm_fit_SOI = nm_fit_SOI
result_right_cortex$mfvb_fit_SOI = mfvb_fit_SOI
result_right_cortex$beta_fit = beta_fit
result_right_cortex$theta_fit = theta_fit
result_right_cortex$time = time

# Evaluate predictive performance of GPR + Normal and GPR + Horseshoe.
y_pred_train = list()
y_pred_test = list()

y_pred_train$NM = as.vector(X_cortex_right_train%*%matrix(result_right_cortex$beta_fit$NM,M,1)) + result_right_cortex$nm_fit$post_mean$betacoef[1,1]
y_pred_test$NM =  as.vector(X_right_test%*%matrix(result_right_cortex$beta_fit$NM,M,1)) + result_right_cortex$nm_fit$post_mean$betacoef[1,1]
y_pred_train$HS =  as.vector(X_cortex_right_train%*%matrix(result_right_cortex$beta_fit$HS,M,1)) + result_right_cortex$hs_fit_SOI$post_mean$betacoef[1,1]
y_pred_test$HS = as.vector(X_right_test%*%matrix(result_right_cortex$beta_fit$HS,M,1)) + result_right_cortex$hs_fit_SOI$post_mean$betacoef[1,1]
y_pred_train$MFVB = as.vector(X_cortex_right_train%*%matrix(result_right_cortex$beta_fit$MFVB,M,1)) + result_right_cortex$mfvb_fit$post_mean$betacoef[1,1]
y_pred_test$MFVB =as.vector(X_right_test%*%matrix(result_right_cortex$beta_fit$MFVB,M,1)) + result_right_cortex$mfvb_fit$post_mean$betacoef[1,1]

result_R2_MSE = list()
result_R2_MSE$R2_train_NM = cor(y_pred_train$NM, y_train)^2
result_R2_MSE$R2_test_NM = cor(y_pred_test$NM, y_test)^2

result_R2_MSE$mse_y_pred_train_NM = mean((y_pred_train$NM - y_train)^2)
result_R2_MSE$mse_y_pred_test_NM = mean((y_pred_test$NM - y_test)^2)

result_R2_MSE$R2_train_HS = cor(y_pred_train$HS, y_train)^2
result_R2_MSE$R2_test_HS = cor(y_pred_test$HS, y_test)^2

result_R2_MSE$mse_y_pred_train_HS = mean((y_pred_train$HS - y_train)^2)
result_R2_MSE$mse_y_pred_test_HS = mean((y_pred_test$HS - y_test)^2)

result_R2_MSE$R2_train_MFVB = cor(y_pred_train$MFVB, y_train)^2
result_R2_MSE$R2_test_MFVB = cor(y_pred_test$MFVB, y_test)^2

result_R2_MSE$mse_y_pred_train_MFVB = mean((y_pred_train$MFVB - y_train)^2)
result_R2_MSE$mse_y_pred_test_MFVB = mean((y_pred_test$MFVB - y_test)^2)

result_right_cortex$y_pred_train = y_pred_train
result_right_cortex$y_pred_test = y_pred_test
result_right_cortex$result_R2_MSE = result_R2_MSE

# Save evaluation of predictive results in a list and create table.
result_eval = matrix(NA, 3, 4)
result_eval[1,1] = result_R2_MSE$R2_train_NM 
result_eval[1,2] = result_R2_MSE$mse_y_pred_train_NM 
result_eval[1,3] = result_R2_MSE$R2_test_NM 
result_eval[1,4] = result_R2_MSE$mse_y_pred_test_NM
result_eval[2,1] = result_R2_MSE$R2_train_HS
result_eval[2,2] = result_R2_MSE$mse_y_pred_train_HS
result_eval[2,3] = result_R2_MSE$R2_test_HS
result_eval[2,4] = result_R2_MSE$mse_y_pred_test_HS
result_eval[3,1] = result_R2_MSE$R2_train_MFVB
result_eval[3,2] = result_R2_MSE$mse_y_pred_train_MFVB
result_eval[3,3] = result_R2_MSE$R2_test_MFVB
result_eval[3,4] = result_R2_MSE$mse_y_pred_test_MFVB

col_names = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')
row_names = c('NM', 'HS', 'MFVB')
colnames(result_eval) = col_names
rownames(result_eval) = row_names

xtable(result_eval, digits = 4)
print(xtable(result_eval, digits = 4), file=paste0(path_results, "table_result_right_cortex.txt"))
write.csv(result_eval, paste0(path_results, 'result_matrix_right.csv'))

save(result_right_cortex, file = sprintf("%sresult_right_cortex.RData", path_results))

### Visualisation ###

# Read in results (left and right cortex)
load(paste0(path_results, 'result_left_cortex.RData'))
load(paste0(path_results, 'result_right_cortex.RData'))

# Save beta map of GPR + Normal model.
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$beta_fit$NM, length(result_left_cortex$beta_fit$NM), 1)
xii$data$cortex_right = matrix(result_right_cortex$beta_fit$NM, length(result_right_cortex$beta_fit$NM), 1)

write_xifti(
  xii, 
  file.path(path_results, "beta_NM.dtseries.nii"), 
  file.path(path_results, "beta_NM_cortex_left.surf.gii"), file.path(path_results, "beta_NM_cortex_right.surf.gii")
)

# Save beta map of GPR + Horseshoe model.
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$beta_fit$HS, length(result_left_cortex$beta_fit$HS), 1)
xii$data$cortex_right = matrix(result_right_cortex$beta_fit$HS, length(result_right_cortex$beta_fit$HS), 1)

write_xifti(
  xii, 
  file.path(path_results, "beta_HS.dtseries.nii"), 
  file.path(path_results, "beta_HS_cortex_left.surf.gii"), file.path(path_results, "beta_HS_cortex_right.surf.gii")
)



# Save binary significance (evaluated via test statistics) map of GPR + Normal model.
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_test_statistic_NM, length(result_left_cortex$significance_test_statistic_NM), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_test_statistic_NM, length(result_right_cortex$significance_test_statistic_NM), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_test_statistic_NM.dtseries.nii"), 
  file.path(path_results, "significance_test_statistic_NM_cortex_left.surf.gii"), file.path(path_results, "significance_test_statistic_NM_cortex_right.surf.gii")
)

# Save binary significance (evaluated via test statistics) map of GPR + Horseshoe model.
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_test_statistic_HS, length(result_left_cortex$significance_test_statistic_HS), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_test_statistic_HS, length(result_right_cortex$significance_test_statistic_HS), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_test_statistic_HS.dtseries.nii"), 
  file.path(path_results, "significance_test_statistic_HS_cortex_left.surf.gii"), file.path(path_results, "significance_test_statistic_HS_cortex_right.surf.gii")
)

# Save binary significance (evaluated via ETI CI) map of GPR + Normal model.
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_CI_ETI_NM, length(result_left_cortex$significance_CI_ETI_NM), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_CI_ETI_NM, length(result_right_cortex$significance_CI_ETI_NM), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_CI_ETI_NM.dtseries.nii"), 
  file.path(path_results, "significance_CI_ETI_NM_cortex_left.surf.gii"), file.path(path_results, "significance_CI_ETI_NM_cortex_right.surf.gii")
)

# Save binary significance (evaluated via ETI CI) map of GPR + Horseshoe model.
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_CI_ETI_HS, length(result_left_cortex$significance_CI_ETI_HS), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_CI_ETI_HS, length(result_right_cortex$significance_CI_ETI_HS), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_CI_ETI_HS.dtseries.nii"), 
  file.path(path_results, "significance_CI_ETI_HS_cortex_left.surf.gii"), file.path(path_results, "significance_CI_ETI_HS_cortex_right.surf.gii")
)

# Save binary significance (evaluated via HPDI CI) map of GPR + Normal model.
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_CI_HPDI_NM, length(result_left_cortex$significance_CI_HPDI_NM), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_CI_HPDI_NM, length(result_right_cortex$significance_CI_HPDI_NM), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_CI_HPDI_NM.dtseries.nii"), 
  file.path(path_results, "significance_CI_HPDI_NM_cortex_left.surf.gii"), file.path(path_results, "significance_CI_HPDI_NM_cortex_right.surf.gii")
)

# Save binary significance (evaluated via HPDI CI) map of GPR + Horseshoe model.
xii = read_xifti(path_xii_file)
xii$data$cortex_left = matrix(result_left_cortex$significance_CI_HPDI_HS, length(result_left_cortex$significance_CI_HPDI_HS), 1)
xii$data$cortex_right = matrix(result_right_cortex$significance_CI_HPDI_HS, length(result_right_cortex$significance_CI_HPDI_HS), 1)

write_xifti(
  xii, 
  file.path(path_results, "significance_CI_HPDI_HS.dtseries.nii"), 
  file.path(path_results, "significance_CI_HPDI_HS_cortex_left.surf.gii"), file.path(path_results, "significance_CI_HPDI_HS_cortex_right.surf.gii")
)