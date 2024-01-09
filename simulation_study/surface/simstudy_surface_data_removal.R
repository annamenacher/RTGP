# Script: Model estimation of GPR + Normal, GPR + Horseshoe, BR + Normal, BR + Horseshoe, Ridge, LASSO, and RTGP-VI
# with splitting training and test data in 10 folds.

options(bitmapType='cairo-png')
dataset = 1
sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# 10-fold cross validation
percent_removal = 10
n_remove = 200

#############
### Paths ###
#############

# Path to folder containing the Human Connectome workbench software
path_wb = ""
# Path to folder to save results in.
path_core = ""
# Path to folder containing the simulated dataset.
path_data = ""
# Path to folder where eigendecompositions are stored.
path_eigen_core = ""
# Path to folder where true beta file is stored. 
path_beta = ''

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
library(bayestestR)
library(dplyr)
library(data.table)

library(ciftiTools)
ciftiTools.setOption('wb_path', path_wb)

#################
### Constants ###
#################

# Kernel hyperparameters
phi = 0.01
nu = 1.4
# Number of basis functions
L = 100

############
### Data ###
############

# Load simulated dataset
load(paste0(path_data, "dataset", toString(dataset), "/data.RData"))

# Input images
X = data$X 
# Output scalars
y = data$y 
# True beta values
beta_true = data$beta_true
# True binary significance values
alpha_true = data$alpha_true 
# True residual variance
true_sigma2 = data$true_sigma2 
# Scalar outputs (test data)
y_test = data$y_test 
# Input images (test data)
img_test = data$img_test 
# Number of vertices
p = M = data$p 
# Number of subjects (training dataset)
n_train = data$n_train 
# Number of subjects (test dataset)
n_test = data$n_test  
# Number of subjects (total dataset)
n = data$n
beta0_true = data$beta0_true 

# Path to eigendecomposition
path_eigen = sprintf("%seigendecomposition_L%s_phi%s_nu%s/", path_eigen_core, L, phi, nu)

# Eigenvalues
lambda_l = read.csv(paste0(path_eigen, "eigenvalues_left_cortex.csv"))[,2]
# Square rooted eigenvalues
sqrt_lambda = sqrt(lambda_l)
# Eigenfunctions
Psi = data.matrix(fread(paste0(path_eigen, "eigenvectors_left_cortex.csv")))
Psi = Psi[,2:ncol(Psi)]
Psi = t(Psi)
# Bases
Bases = Psi*sqrt_lambda
# Transformed input images(multiplied with bases) 
Z = cbind(1,t(Bases%*%t(X)))
Z_test = cbind(1,t(Bases%*%img_test))

# Split data in CV partitions
idx = c(((sim - 1)*n_remove + 1) : (n_remove*sim))
# Training data
X = X[-idx, ]
y = y[-idx]
Z = Z[-idx, ]
# Test data
X_test = X[idx,]
data$img_test = img_test = t(X[idx,])
y_test = data$y_test = y[idx]
Z_test = Z[idx,]

#############
### Paths ###
#############

path_results = path =  sprintf("%s/N%s_L%s_phi%s_nu%s/", path_core, n_train, L, phi, nu)
dir.create(file.path(path_results), showWarnings = FALSE)

path_results = path =  sprintf("%s/N%s_L%s_phi%s_nu%s/dataset%s/", path_core, n_train, L, phi, nu, dataset)
dir.create(file.path(path_results), showWarnings = FALSE)

path_results_BR =  sprintf("%s/N%s_L%s_phi%s_nu%s/dataset%s/BR_CV/", path_core, n_train, L, phi, nu, dataset)
dir.create(file.path(path_results_BR), showWarnings = FALSE)

path_results_LASSO =  sprintf("%s/N%s_L%s_phi%s_nu%s/dataset%s/Ridge_LASSO_CV/", path_core, n_train, L, phi, nu, dataset)
dir.create(file.path(path_results_LASSO), showWarnings = FALSE)

path_results_NM =  sprintf("%s/N%s_L%s_phi%s_nu%s/dataset%s/GPR_CV/", path_core, n_train, L, phi, nu, dataset)
dir.create(file.path(path_results_NM), showWarnings = FALSE)

path_results_RTGP_VI =  sprintf("%s/N%s_L%s_phi%s_nu%s/dataset%s/RTGP_VI_CV/", path_core, n_train, L, phi, nu, dataset)
dir.create(file.path(path_results_RTGP_VI), showWarnings = FALSE)

path_results_all =  sprintf("%s/N%s_L%s_phi%s_nu%s/dataset%s/all/", path_core, n_train, L, phi, nu, dataset)
dir.create(file.path(path_results_all), showWarnings = FALSE)


###########
### GPR ###
###########

# Set current path.
path = path_results_NM
# Create list to save results in.
result_time = list()
result_nm = list()

# Start time for model estimation.
time.start = Sys.time()

# Fit GPR + Normal, GPR + Horseshoe and GPR + MFVB model.
hs_fit_SOI = fast_horseshoe_lm(y,Z, mcmc_sample = 5000, burnin = 2000)
nm_fit_SOI = fast_normal_lm(y,Z, mcmc_sample = 5000, burnin = 2000)
mfvb_fit_SOI = fast_mfvb_normal_lm(y,Z)

# End time for model estimation.
time.end = Sys.time()
time = time.end - time.start

# Save GPR model estimation time.
result_time$time_nm = time

# Concatenate estimated thetas in a dataframe.
theta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-1],
                       NM = nm_fit_SOI$post_mean$betacoef[-1],
                       MFVB = mfvb_fit_SOI$post_mean$betacoef[-1])

# Concatenate estimated betas in a dataframe.
beta_fit = data.frame(HS = crossprod(Bases,hs_fit_SOI$post_mean$betacoef[-1]),
                      NM = crossprod(Bases,nm_fit_SOI$post_mean$betacoef[-1]),
                      MFVB = crossprod(Bases,mfvb_fit_SOI$post_mean$betacoef[-1]))

# Save results in list.
result_nm$beta_fit = beta_fit
result_nm$theta_fit = theta_fit
result_nm$hs_fit_SOI = hs_fit_SOI
result_nm$nm_fit_SOI = nm_fit_SOI
result_nm$time = time
save(result_nm, file = sprintf("%sresult_nm.RData", path_results_NM))

# Save estimated beta of GPR + Normal.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(beta_fit$NM, length(beta_fit$NM), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "beta_NM.dtseries.nii"), 
  file.path(path_results_NM, "beta_NM_left.surf.gii"), file.path(path_results_NM, "beta_NM_cortex_right.surf.gii")
)

# Save estimated beta of GPR + Horseshoe
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(beta_fit$HS, length(beta_fit$HS), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "beta_HS.dtseries.nii"), 
  file.path(path_results_NM, "beta_HS_left.surf.gii"), file.path(path_results_NM, "beta_HS_cortex_right.surf.gii")
)

# Plot scatterplot of true beta vs estimated beta of GPR + Normal.
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

# Plot scatterplot of true beta vs estimated beta of GPR + Horseshoe
df = data.frame(beta_true = beta_true, beta_method = beta_fit$HS)
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta True') +
  ylab('Beta HS') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
p_plot
ggsave(p_plot, filename=paste0(path_results_NM, "scatterplot_beta_true_vs_beta_HS.png"), type="cairo-png")

# Plot scatterplot of estimated beta of GPR + Normal vs estimated beta of GPR + Horseshoe
df = data.frame(beta_true = beta_fit$NM, beta_method = beta_fit$HS)
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta NM') +
  ylab('Beta HS') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
p_plot
ggsave(p_plot, filename=paste0(path_results_NM, "scatterplot_beta_NM_vs_beta_HS.png"), type="cairo-png")


# Determine binary significance of GPR + Normal.
# Binary significance: 95%-Equidistant credible interval vertices that do not contain 0
beta_NM = t(crossprod(Bases, nm_fit_SOI$mcmc$betacoef[(2):(L+1),]))
x=eti(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_NM = ifelse(x2==T,0,1)
# Binary significance: 95%-Highest posterior density credible interval vertices that do not contain 0
x=hdi(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_NM = ifelse(x2==T,0,1)
# Binary significance: Test statistics = posterior beta mean / posterior beta variance -> FDR correction -> Significant vertices determined by p-values lower
# than 5%
t_beta_mean_NM = apply(beta_NM, 2, mean)
t_beta_sd_NM = apply(beta_NM, 2, sd)
t = t_beta_mean_NM / t_beta_sd_NM
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_NM = ifelse(p_fdr < 0.05, 1, 0)

# Determine binary significance of GPR + Horseshoe.
# Binary significance: 95%-Equidistant credible interval vertices that do not contain 0
beta_HS = t(crossprod(Bases, hs_fit_SOI$mcmc$betacoef[(2):(L+1),]))
x=eti(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_HS = ifelse(x2==T,0,1)
# Binary significance: 95%-Highest posterior density credible interval vertices that do not contain 0
x=hdi(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_HS = ifelse(x2==T,0,1)
# Binary significance: Test statistics = posterior beta mean / posterior beta variance -> FDR correction -> Significant vertices determined by p-values lower
# than 5%
t_beta_mean_HS = apply(beta_HS, 2, mean)
t_beta_sd_HS = apply(beta_HS, 2, sd)
t = t_beta_mean_HS / t_beta_sd_HS
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_HS = ifelse(p_fdr < 0.05, 1, 0)


# Save binary significance results of HPDI evaluation of GPR + Normal.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(CI_HPDI_binary_NM, length(beta_fit$NM), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "significance_NM.dtseries.nii"), 
  file.path(path_results_NM, "significance_NM_left.surf.gii"), file.path(path_results_NM, "significance_NM_cortex_right.surf.gii")
)

# Save binary significance results of HPDI evaluation of GPR + Horseshoe.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(CI_HPDI_binary_HS, length(beta_fit$HS), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "significance_HS.dtseries.nii"), 
  file.path(path_results_NM, "significance_HS_left.surf.gii"), file.path(path_results_NM, "significance_HS_cortex_right.surf.gii")
)

# True threshold
threshold_true = 0
# True binary significance mask 
mask = data$alpha_true

# Evaluation of bias of parameter estimates (GPR + Normal)
print('Bias: NM')
(bias_nm_no_effect = mean((beta_fit$NM - beta_true)[mask == 0]))
(bias_nm_effect = mean((beta_fit$NM - beta_true)[mask == 1]))
(bias_nm_total = mean((beta_fit$NM - beta_true)))

# Evaluation of bias of parameter estimates (GPR + Horseshoe)
print('Bias: HS')
(bias_hs_no_effect = mean((beta_fit$HS - beta_true)[mask == 0]))
(bias_hs_effect = mean((beta_fit$HS - beta_true)[mask == 1]))
(bias_hs_total = mean((beta_fit$HS - beta_true)))

# Evaluation of MSE of parameter estimates (GPR + Normal)
print('MSE: NM')
(mse_nm_no_effect = mean((beta_fit$NM - beta_true)[mask == 0]^2))
(mse_nm_effect = mean((beta_fit$NM - beta_true)[mask == 1]^2))
(mse_nm_total = mean((beta_fit$NM - beta_true)^2))

# Evaluation of MSE of parameter estimates (GPR + Horseshoe)
print('MSE: HS')
(mse_hs_no_effect = mean((beta_fit$HS - beta_true)[mask == 0]^2))
(mse_hs_effect = mean((beta_fit$HS - beta_true)[mask == 1]^2))
(mse_hs_total = mean((beta_fit$HS - beta_true)^2))

# Evaluation of inference results (GPR + Normal)
print('Inference results: NM')

# Binary significance: 95%-Equidistant credible interval vertices that do not contain 0
posterior = t(crossprod(Bases,nm_fit_SOI$mcmc$betacoef[-1,]))
x=eti(data.frame(posterior))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary = ifelse(x2==T,0,1)
# Binary significance: 95%-Highest posterior density credible interval vertices that do not contain 0
x=hdi(data.frame(posterior))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary = ifelse(x2==T,0,1)
# Binary significance: Test statistics = posterior beta mean / posterior beta variance -> FDR correction -> Significant vertices determined by p-values lower
# than 5%
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

TP = CI_ETI_binary + mask
TP = length(TP[TP == 2])
TN = CI_ETI_binary + mask
TN = length(TN[TN == 0])
FP = CI_ETI_binary - mask
FP = length(FP[FP == 1])
FN = CI_ETI_binary - mask
FN = length(FN[FN == -1])
(TPR_CI_ETI_nm = TP/(TP + FN))
(TDR_CI_ETI_nm = TP/(TP + FP))
(FPR_CI_ETI_nm = FP/(FP + TN))
(FDR_CI_ETI_nm = FP/(FP + TP))

TP = CI_HPDI_binary + mask
TP = length(TP[TP == 2])
TN = CI_HPDI_binary + mask
TN = length(TN[TN == 0])
FP = CI_HPDI_binary - mask
FP = length(FP[FP == 1])
FN = CI_HPDI_binary - mask
FN = length(FN[FN == -1])
(TPR_CI_HPDI_nm = TP/(TP + FN))
(TDR_CI_HPDI_nm = TP/(TP + FP))
(FPR_CI_HPDI_nm = FP/(FP + TN))
(FDR_CI_HPDI_nm = FP/(FP + TP))

# Evaluation of inference results (GPR + Horseshoe)
print('Inference results: HS')

# Binary significance: 95%-Equidistant credible interval vertices that do not contain 0
posterior = t(crossprod(Bases,hs_fit_SOI$mcmc$betacoef[-1,]))
x=eti(data.frame(posterior))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary = ifelse(x2==T,0,1)
# Binary significance: 95%-Highest posterior density credible interval vertices that do not contain 0
x=hdi(data.frame(posterior))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary = ifelse(x2==T,0,1)
# Binary significance: Test statistics = posterior beta mean / posterior beta variance -> FDR correction -> Significant vertices determined by p-values lower
# than 5%
t_beta_mean = apply(crossprod(Bases,hs_fit_SOI$mcmc$betacoef[-1,]), 1, mean)
t_beta_sd = apply(crossprod(Bases,hs_fit_SOI$mcmc$betacoef[-1,]), 1, sd)
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
(TPR_hs = TP/(TP + FN))
(TDR_hs = TP/(TP + FP))
(FPR_hs = FP/(FP + TN))
(FDR_hs = FP/(FP + TP))

TP = CI_ETI_binary + mask
TP = length(TP[TP == 2])
TN = CI_ETI_binary + mask
TN = length(TN[TN == 0])
FP = CI_ETI_binary - mask
FP = length(FP[FP == 1])
FN = CI_ETI_binary - mask
FN = length(FN[FN == -1])
(TPR_CI_ETI_hs = TP/(TP + FN))
(TDR_CI_ETI_hs = TP/(TP + FP))
(FPR_CI_ETI_hs = FP/(FP + TN))
(FDR_CI_ETI_hs = FP/(FP + TP))

TP = CI_HPDI_binary + mask
TP = length(TP[TP == 2])
TN = CI_HPDI_binary + mask
TN = length(TN[TN == 0])
FP = CI_HPDI_binary - mask
FP = length(FP[FP == 1])
FN = CI_HPDI_binary - mask
FN = length(FN[FN == -1])
(TPR_CI_HPDI_hs = TP/(TP + FN))
(TDR_CI_HPDI_hs = TP/(TP + FP))
(FPR_CI_HPDI_hs = FP/(FP + TN))
(FDR_CI_HPDI_hs = FP/(FP + TP))


# Evaluation of R2 (train) prediction results (GPR + Normal)
print('R2 Train: NM')
pred_nm = nm_fit_SOI$post_mean$betacoef[1] + crossprod(t(X),beta_fit$NM)
(R2_nm = cor(pred_nm,y)^2)

# Evaluation of R2 (train) prediction results (GPR + Horseshoe)
print('R2 Train: HS')
pred_hs = hs_fit_SOI$post_mean$betacoef[1] + crossprod(t(X),beta_fit$HS)
(R2_hs = cor(pred_hs,y)^2)

# Evaluation of bias (train) prediction results (GPR + Normal)
print('Bias Prediction Train: NM')
(bias_y_pred_nm = mean(pred_nm - y))

# Evaluation of bias (train) prediction results (GPR + Horseshoe)
print('Bias Prediction Train: HS')
(bias_y_pred_hs = mean(pred_hs - y))

# Evaluation of MSE (train) prediction results (GPR + Normal)
print('MSE Prediction Train: NM')
(mse_y_pred_nm = mean((pred_nm - y)^2))

# Evaluation of MSE (train) prediction results (GPR + Horseshoe)
print('MSE Prediction Train: HS')
(mse_y_pred_hs = mean((pred_hs - y)^2))

# Evaluation of R2 (test) prediction results (GPR + Normal)
print('R2 Test: NM')
pred_nm_test = nm_fit_SOI$post_mean$betacoef[1] + crossprod(data$img_test,beta_fit$NM)
(R2_nm_test = cor(pred_nm_test,data$y_test)^2)

# Evaluation of R2 (test) prediction results (GPR + Horseshoe)
print('R2 Test: HS')
pred_hs_test = hs_fit_SOI$post_mean$betacoef[1] + crossprod(data$img_test,beta_fit$HS)
(R2_hs_test = cor(pred_hs_test,data$y_test)^2)

# Evaluation of bias (test) prediction results (GPR + Normal)
print('Bias Prediction Test: NM')
(bias_y_pred_nm_test = mean(pred_nm_test - data$y_test))

# Evaluation of bias (test) prediction results (GPR + Horseshoe)
print('Bias Prediction Test: HS')
(bias_y_pred_hs_test = mean(pred_hs_test - data$y_test))

# Evaluation of MSE (test) prediction results (GPR + Normal)
print('MSE Prediction Test: NM')
(mse_y_pred_nm_test = mean((pred_nm_test - data$y_test)^2))

# Evaluation of MSE (test) prediction results (GPR + Horseshoe)
print('MSE Prediction Test: HS')
(mse_y_pred_hs_test = mean((pred_hs_test - data$y_test)^2))



###############
### RTGP-VI ###
###############

# Set to current path.
path = path_results_RTGP_VI

# Maximum number of iterations. 
n_iter = 10000
# Convergence criterion: Stop if difference in ELBO from the previous iteration is smaller than convergence threshold.
eps = 0.00001

# Standard normal density function.
phi_function = function(x){
  x = (1/sqrt(2*pi)) * exp((-0.5)*x^2)
  return(x)
}

# Log-likelihood of RTGP model 
log_likelihood = function(y, Psi, lambda_l, theta, beta0, X, sigma_epsilon, alpha, threshold){
  beta = t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L) %*% theta
  beta_thresh = matrix(beta * ifelse(abs(alpha) > threshold, 1, 0), p, 1)
  mu = beta0 + X%*%beta_thresh
  ll = sum(dnorm(y, mu, sigma_epsilon, log=TRUE))
  return(ll)
}

# Joint density (consisting of prior and likelihood) of RTGP model
joint_density = function(){
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
# Hyperparameters of variance priors
a_beta = b_beta = 0.001
a_alpha = b_alpha = 0.001
a_epsilon = b_epsilon = 0.001
# Bounds of Uniform prior on threshold lambda (lower bound: 25% quantile of estimated beta parameters of GPR + Normal model,
# upper bound: 90% quantile of estimated beta parameters of GPR + Normal)
t_min = max(abs(beta_fit$NM)) * 0.25
t_max =  max(abs(beta_fit$NM)) * 0.9
# Threshold options for discrete Uniform prior on threshold.
threshold_options = seq(t_min , t_max, length.out = 30)
# Hyperparameter of prior on intercept beta0
sigma_beta0 = 10

# Start time for RTGP model estimation
time.start = Sys.time()

# Initialization (with GPR + Normal or fixed quantities)
threshold = (t_min+t_max)/2
theta = matrix(theta_fit$NM,L,1)
beta0 =  nm_fit_SOI$post_mean$betacoef[1]
beta = beta_fit$NM
sigma_epsilon = sigma_epsilon_prior = sqrt(nm_fit_SOI$post_mean$sigma2_eps)
sigma_alpha = sigma_alpha_prior = 1
sigma_beta = sigma_beta_prior = 1
alpha = rnorm(M, beta_fit$NM, sigma_alpha)

# Additional quantities to initialise
# Inverse of variance parameters
exp_1_sigma_epsilon2 = 1/sigma_epsilon^2
exp_1_sigma_alpha2 = 1/sigma_alpha^2
exp_1_sigma_beta2 = 1/sigma_beta^2

# Posterior variance of variational density on the intercept beta0
Sigma_beta0 = 1/(n_train*exp_1_sigma_epsilon2 + (1/sigma_beta0^2))
# Posterior variance of variational density on basis coefficients theta
Sigma_theta = solve(exp_1_sigma_epsilon2 * (t(X%*%diag(ifelse(abs(as.vector(alpha))>threshold,1,0),p)%*%t(Psi)%*%diag(sqrt(lambda_l),L))%*%X%*%diag(ifelse(abs(as.vector(alpha))>threshold,1,0),p)%*%t(Psi)%*%diag(sqrt(lambda_l),L) + diag(exp_1_sigma_beta2,L) + exp_1_sigma_alpha2*t(t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)) %*% (t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L))))
# Expected value of alpha
exp_I_abs_alpha_greater_lamda = ifelse(abs(alpha)>threshold, 1, 0)
# Expected value of beta^2
exp_beta2 = diag((t(Psi) %*% diag(sqrt(lambda_l),L)) %*% (Sigma_theta + theta %*% t(theta)) %*% t(t(Psi) %*% diag(sqrt(lambda_l),L)))
# Re-weighting of the weights of the mixture distribution of alpha
c = exp((-0.5) * exp_1_sigma_epsilon2 * 2 * t(t(y - beta0) %*% (-X)) * beta  - 0.5 * exp_1_sigma_epsilon2 * colSums(X^2) * exp_beta2-0.5 * exp_1_sigma_epsilon2 * colSums(sweep(X, 2, beta, "*")*(sweep(-sweep(X, 2, (beta * exp_I_abs_alpha_greater_lamda), "*") , 1, X %*% (beta * exp_I_abs_alpha_greater_lamda), "+"))))
# Weights of the mixture distribution on alpha
w = matrix(NA, 3, p)
w[1,] = as.vector(pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) * c
w[2,] = as.vector(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) 
w[3,] = as.vector(1 - pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2))) * c
w_sum = apply(w, 2, sum)
w = sweep(w, 2, w_sum, FUN="/")

# Expected value of truncated distribution on alpha
exp_alpha = matrix(NA, 3, p)
exp_alpha[1,] = etruncnorm(a = -Inf, b = -threshold, mean = beta, sd = sigma_alpha)
exp_alpha[2,] = etruncnorm(a = -threshold, b = threshold, mean = beta, sd = sigma_alpha)
exp_alpha[3,] = etruncnorm(a = threshold, b = Inf, mean = beta, sd = sigma_alpha)

# Initialise difference between ELBOs
diff = 100
# Start counter of iterations
iter = 0
# Initialise ELBO 
ELBO = -100000000000000000

# Create vectors and matrices to monitor quantities, such as log-likelihood, joint, difference in ELBOs, ... across iterations
log_likelihood_chain = numeric(n_iter)
joint_chain = numeric(n_iter)
diff_chain = numeric(n_iter)
threshold_chain = numeric(n_iter)
ELBO_chain = numeric(n_iter)
beta0_chain = numeric(n_iter)
theta_chain = matrix(NA, n_iter, L)
alpha_chain = matrix(NA, n_iter, M)
beta_chain = matrix(NA, n_iter, M)
beta_threshold_chain = matrix(NA, n_iter, M)
alpha_threshold_chain = matrix(NA, n_iter, M)
sigma_epsilon_chain = numeric(n_iter)
sigma_alpha_chain = numeric(n_iter)
sigma_beta_chain = numeric(n_iter)
exp_I_abs_alpha_greater_lamda_chain = matrix(NA, n_iter, M)

while((diff > eps) && (iter != n_iter)){
  tryCatch({
    
    gc()
    # Update iteration counter.
    iter = iter + 1
    print(iter)
    # Update the quantities of the previous iteration (beta, alpha, ELBO)
    beta_old = beta
    alpha_old = alpha
    ELBO_old = ELBO
    # Calculate log-likelihood
    ll = log_likelihood(y, Psi, lambda_l, theta, beta0, X, sigma_epsilon, alpha, threshold)
    # Calculate joint density
    joint = joint_density()
    
    # Update beta0 
    exp_I_abs_alpha_greater_lamda = w[1,] + w[3,]
    Sigma_beta0 = 1/(n_train*exp_1_sigma_epsilon2 + (1/(sigma_beta0^2)))
    beta0 = as.vector( Sigma_beta0 * (exp_1_sigma_epsilon2 * sum(y - X%*%(exp_I_abs_alpha_greater_lamda * beta))))
    
    # Update theta
    Sigma_theta = solve(exp_1_sigma_epsilon2 * t(X %*% diag(as.vector(exp_I_abs_alpha_greater_lamda), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) %*% (X %*% diag(as.vector(exp_I_abs_alpha_greater_lamda), p) %*% t(Psi) %*% diag(as.vector(sqrt(lambda_l)), L)) + diag(exp_1_sigma_beta2, L) + exp_1_sigma_alpha2*t(t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)) %*% (t(Psi)%*%diag(as.vector(sqrt(lambda_l)), L)))
    theta = Sigma_theta %*% t(exp_1_sigma_epsilon2 * t(y - beta0) %*% (X%*%diag(as.vector(exp_I_abs_alpha_greater_lamda),p)%*%t(Psi)%*%diag(as.vector(sqrt(lambda_l)),L)) + exp_1_sigma_alpha2 * t(alpha)%*%t(Psi)%*%diag(as.vector(sqrt(lambda_l)),L))  
    # Calculate beta from the updated thetas.
    beta = t(Psi)%*%(sqrt(lambda_l) * theta)
    gc()
    
    # Update alpha
    exp_beta2 = diag((t(Psi) %*% diag(sqrt(lambda_l),L)) %*% (Sigma_theta + theta %*% t(theta)) %*% t(t(Psi) %*% diag(sqrt(lambda_l),L)))
    c = exp((-0.5) * exp_1_sigma_epsilon2 * 2 * t(t(y - beta0 ) %*% (-X)) * beta  - 0.5 * exp_1_sigma_epsilon2 * colSums(X^2) * exp_beta2-0.5 * exp_1_sigma_epsilon2 * colSums(sweep(X, 2, beta, "*")*(sweep(-sweep(X, 2, (beta * exp_I_abs_alpha_greater_lamda), "*") , 1, X %*% (beta * exp_I_abs_alpha_greater_lamda), "+"))))
    
    w = matrix(NA, 3, p)
    w[1,] = as.vector(pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2))) * c
    w[2,] = as.vector(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm(-(threshold + beta)*sqrt(exp_1_sigma_alpha2)))
    w[3,] = as.vector(1 - pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2))) * c
    w_sum = apply(w, 2, sum)
    w = sweep(w, 2, w_sum, FUN="/")
    
    exp_alpha = matrix(NA, 3, p)
    exp_alpha[1,] = beta - sigma_alpha * (phi_function((-threshold - beta) / sigma_alpha) / (pnorm((-threshold - beta) / sigma_alpha) + 10^(-10)))
    exp_alpha[2,] = beta - sigma_alpha * (phi_function((threshold - beta) / sigma_alpha) - phi_function((-threshold - beta) / sigma_alpha)) / (pnorm((threshold - beta) / sigma_alpha) - pnorm((-threshold - beta) / sigma_alpha) + 10^(-10))
    exp_alpha[3,] = beta + sigma_alpha * phi_function((threshold - beta) / sigma_alpha) / (1 - pnorm((threshold - beta) / sigma_alpha) + 10^(-10))
    
    alpha = matrix(w[1,] * exp_alpha[1,] + w[2,] * exp_alpha[2,] + w[3,] * exp_alpha[3,], p, 1)   
    
    # Update threshold
    threshold_t = numeric(length(threshold_options))
    for(t in 1:length(threshold_options)){
      threshold_t[t] = (-0.5) * exp_1_sigma_epsilon2 * sum((y - beta0 - X%*%(beta * ifelse(abs(alpha) > threshold_options[t], 1, 0)))^2)
    }
    
    logProb = threshold_t - max(threshold_t)
    Prob = exp(logProb)/sum(exp(logProb))
    threshold = sum(Prob*threshold_options)
    
    # Update sigma_epsilon
    var_beta_alpha = exp_beta2 * exp_I_abs_alpha_greater_lamda - beta^2 * exp_I_abs_alpha_greater_lamda^2
    exp_beta_alpha_X = sum(sweep(X^2, 2, var_beta_alpha, "*")) + sum(diag((X%*%(beta * exp_I_abs_alpha_greater_lamda)) %*% t(X%*%(beta * exp_I_abs_alpha_greater_lamda))))
    
    # Update sigma_epsilon
    exp_1_sigma_epsilon2 = as.vector((n_train + a_epsilon) / (b_epsilon + 0.5 * sum((y - beta0  - X%*%(beta * exp_I_abs_alpha_greater_lamda))^2)))
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
    
    # Calculate thresholded beta and alpha
    beta_thresholded =  as.vector( beta * ifelse(abs(as.vector(alpha))>threshold, 1, 0) )
    alpha_thresholded = as.vector( ifelse(abs(alpha)>threshold, 1, 0) )
    
    # Print the number of significant vertices
    print(table(ifelse(exp_I_abs_alpha_greater_lamda>0.5,1,0)))
    
    gc()
    
    # Log quantities, such as log-likelihood, joint, threshold, ...
    log_likelihood_chain[(iter)] = ll
    joint_chain[iter] = joint
    threshold_chain[iter] = threshold
    
    # Calculate ELBO
    exp_ln_sigma_alpha2 = log(0.5 * sum((alpha-beta)^2) + b_alpha) - digamma(p/2 + a_alpha)
    exp_ln_sigma_epsilon2 = log(0.5 * (t(y)%*%y + n_train * (Sigma_beta0 + beta0^2) + exp_beta_alpha_X - 2 * t(y-beta0)%*%(X%*%(beta * exp_I_abs_alpha_greater_lamda)) - sum(2* beta0 * y)) + b_epsilon) - digamma(n_train/2 + a_epsilon)
    exp_ln_sigma_beta2 = log(0.5 * sum(diag(Sigma_theta) + theta^2) + b_beta) - digamma(L/2 + a_beta)
    part = log(pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part1 = sum(w[1,] * part)
    part = -(-threshold - beta)*sqrt(exp_1_sigma_alpha2) * phi_function((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) / (2*pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part1 = alpha_part1 + sum(w[1,] * part)
    part = log(pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm((-threshold-beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part2 = sum(w[2,] * part)
    part = ((-threshold - beta)*sqrt(exp_1_sigma_alpha2) * phi_function((-threshold - beta)*sqrt(exp_1_sigma_alpha2)) - (threshold - beta)*sqrt(exp_1_sigma_alpha2)*phi_function((threshold - beta)*sqrt(exp_1_sigma_alpha2))) / (2 * (pnorm((threshold - beta)*sqrt(exp_1_sigma_alpha2)) - pnorm((-threshold - beta)*sqrt(exp_1_sigma_alpha2))) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part2 = alpha_part2 + sum(w[2,] * part)
    part = log(1 - pnorm((threshold-beta)*sqrt(exp_1_sigma_alpha2)) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part3 = sum(w[3,] * part)
    part = ((threshold-beta)*sqrt(exp_1_sigma_alpha2)*phi_function((threshold-beta)*sqrt(exp_1_sigma_alpha2))) / (2*(1-pnorm((threshold-beta)*sqrt(exp_1_sigma_alpha2))) + 10^(-10))
    part[part == -Inf] = 0
    part[part == Inf] = 0
    alpha_part3 = alpha_part3 + sum(w[3,] * part)
    
    ELBO = - (n_train/2) * log(2*pi) - (n_train/2) * exp_ln_sigma_epsilon2 - 0.5 * exp_1_sigma_epsilon2 * sum((y - beta0  - X%*%beta_thresholded )^2) -
      0.5 * log(2*pi) - 0.5 * log(sigma_beta0^2) - (1/sigma_beta0^2) * (Sigma_beta0 + beta0^2)  -
      (L/2) * log(2*pi) - (L/2) * exp_ln_sigma_beta2 -
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
    
    # Calculate difference in ELBO
    diff = as.vector(ELBO) - as.vector(ELBO_old)
    diff_chain[iter] = diff
    
    # Log more quantities 
    ELBO_chain[iter] = ELBO
    beta0_chain[iter] = beta0
    theta_chain[iter,] = theta
    beta_chain[iter,] = beta
    alpha_chain[iter,] = alpha
    beta_threshold_chain[iter,] = beta_thresholded
    alpha_threshold_chain[iter,] = alpha_thresholded
    sigma_epsilon_chain[iter] = sigma_epsilon
    sigma_alpha_chain[iter] = sigma_alpha
    sigma_beta_chain[iter] = sigma_beta
    exp_I_abs_alpha_greater_lamda_chain[iter,] = exp_I_abs_alpha_greater_lamda
    
    # Save results every 10 iterations.
    if(iter %% 10 == 0){
      
      # Save trace plot of ELBO.
      df = data.frame('ELBO' = ELBO_chain[2:iter],
                      'Iteration' = 2:iter)
      p_plot = ggplot(df, aes(x = Iteration, y = ELBO)) + 
        geom_line(aes()) +
        theme_classic() + 
        theme(panel.background = element_rect(fill = "white",
                                              colour = "white",
                                              size = 1, linetype = "solid"),
              plot.background = element_rect(fill = "white")) +
        ggtitle('ELBO')
      ggsave(p_plot, filename=paste0(path, "ELBO_chain_left.png"), type="cairo-png")
      
      path_out = path 
      
      # Save scatterplot comparing true beta values to estimated beta values.
      df = data.frame(beta_true = beta_true, beta_method = beta)
      p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
        geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
        geom_abline(intercept = 0) +
        xlim(min(df), max(df)) +
        ylim(min(df), max(df)) +
        xlab('Beta True') +
        ylab('Beta RTGP') +
        theme(text = element_text(size=35)) +
        guides(fill = guide_legend(title=''))
      ggsave(p_plot, filename=paste0(path_out, "RTGP_true_beta_rtgp.png"), type="cairo-png")
      
      # Save scatterplot comparing true beta values to estimated thresholded beta values.
      df = data.frame(beta_true = beta_true, beta_method = beta_thresholded)
      p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
        geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
        geom_abline(intercept = 0) +
        xlim(min(df), max(df)) +
        ylim(min(df), max(df)) +
        xlab('Beta True') +
        ylab('Beta RTGP') +
        theme(text = element_text(size=35)) +
        guides(fill = guide_legend(title=''))
      ggsave(p_plot, filename=paste0(path_out, "RTGP_true_beta_rtgp_thresholded.png"), type="cairo-png")
      
      # Save scatterplot comparing estimated beta values of GPR + Normal to estimated beta values of RTGP model.
      df = data.frame(beta_true = beta_fit$NM, beta_method = beta)
      p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
        geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
        geom_abline(intercept = 0) +
        xlim(min(df), max(df)) +
        ylim(min(df), max(df)) +
        xlab('Beta NM') +
        ylab('Beta RTGP') +
        theme(text = element_text(size=35)) +
        guides(fill = guide_legend(title=''))
      ggsave(p_plot, filename=paste0(path_out, "RTGP_NM_beta_rtgp.png"), type="cairo-png")
      
      # Save scatterplot comparing estimated beta values of GPR + Hprseshoe to estimated beta values of RTGP model.
      df = data.frame(beta_true = beta_fit$HS, beta_method = beta)
      p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
        geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
        geom_abline(intercept = 0) +
        xlim(min(df), max(df)) +
        ylim(min(df), max(df)) +
        xlab('Beta HS') +
        ylab('Beta RTGP') +
        theme(text = element_text(size=35)) +
        guides(fill = guide_legend(title=''))
      ggsave(p_plot, filename=paste0(path_out, "RTGP_HS_beta_rtgp.png"), type="cairo-png")
      
     
      # Save expected alpha map 
      xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
      xii$data$cortex_left = matrix(exp_I_abs_alpha_greater_lamda, M, 1)
      xii$data$cortex_right = matrix(0, 1803, 1)
      
      write_xifti(
        xii, 
        file.path(path_out, "alpha_thresholded_RTGP_VI.dtseries.nii"), 
        file.path(path_out, "alpha_thresholded_RTGP_VI_left.surf.gii"), file.path(path_out, "alpha_thresholded_RTGP_VI_right.surf.gii")
      )
      
      # Save beta map 
      xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
      xii$data$cortex_left = matrix(beta, M, 1)
      xii$data$cortex_right = matrix(0, 1803, 1)
      
      write_xifti(
        xii, 
        file.path(path_out, "beta_RTGP_VI.dtseries.nii"), 
        file.path(path_out, "beta_RTGP_VI_left.surf.gii"), file.path(path_out, "beta_RTGP_VI_right.surf.gii")
      )
      
      # Save thresholded beta map 
      xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
      xii$data$cortex_left = matrix(beta_thresholded, M, 1)
      xii$data$cortex_right = matrix(0, 1803, 1)
      
      write_xifti(
        xii, 
        file.path(path_out, "beta_thresholded_RTGP_VI.dtseries.nii"), 
        file.path(path_out, "beta_thresholded_RTGP_VI_left.surf.gii"), file.path(path_out, "beta_thresholded_RTGP_VI_right.surf.gii")
      )
      
      # Save the other results in list
      result = list()
      result$beta0 = beta0_chain
      result$theta = theta_chain
      result$beta = beta_chain
      result$alpha = alpha_chain
      result$sigma_epsilon = sigma_epsilon_chain
      result$sigma_beta = sigma_beta_chain
      result$sigma_alpha = sigma_alpha_chain
      result$beta_thresholded = beta_threshold_chain
      result$alpha_thresholded = alpha_threshold_chain
      result$log_likelihood = log_likelihood_chain
      result$joint = joint_chain
      result$diff_chain = diff_chain
      result$threshold = threshold_chain
      result$ELBO = ELBO_chain
      result$iter = iter
      result$exp_I_abs_alpha_greater_lamda = exp_I_abs_alpha_greater_lamda_chain
      
      save(result, file = sprintf("%sresult_rtgp_vi_left_cortex.RData", path_out))
      rm(result)
      gc()
      
    }
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  time.end = Sys.time()
  time = time.end - time.start
}

# Index of iteration with maximum ELBO value (iteration of convergence)
iter_ELBO = which(ELBO_chain == max(ELBO_chain[ELBO_chain!=0]))

# Evaluation of bias of parameter estimates
print('Bias: RTGP')
(bias_rtgp_no_effect = mean((beta_threshold_chain[iter_ELBO,] - beta_true)[mask == 0]))
(bias_rtgp_effect = mean((beta_threshold_chain[iter_ELBO,] - beta_true)[mask == 1]))
(bias_rtgp_total = mean((beta_threshold_chain[iter_ELBO,] - beta_true)))

# Evaluation of MSE of parameter estimates
print('MSE: RTGP')
(mse_rtgp_no_effect = mean((beta_threshold_chain[iter_ELBO,] - beta_true)[mask == 0]^2))
(mse_rtgp_effect = mean((beta_threshold_chain[iter_ELBO,] - beta_true)[mask == 1]^2))
(mse_rtgp_total = mean((beta_threshold_chain[iter_ELBO,] - beta_true)^2))

# Evaluation of inference results 
print('Inference results: RTGP')

t_bin_fdr = ifelse(exp_I_abs_alpha_greater_lamda_chain[iter_ELBO,] > 0.5, 1, 0)
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

# Evaluation of R2 (train) of prediction results
print('R2 Train: RTGP')
pred_rtgp = beta0_chain[iter_ELBO] + crossprod(t(X),beta_threshold_chain[iter_ELBO,])
(R2_rtgp = cor(pred_rtgp,y)^2)

# Evaluation of bias (train) of prediction results
print('Bias Prediction Train: RTGP')
(bias_y_pred_rtgp = mean(pred_rtgp - y))

# Evaluation of MSE (train) of prediction results
print('MSE Prediction Train: RTGP')
(mse_y_pred_rtgp = mean((pred_rtgp - y)^2))

# Evaluation of R2 (test) of prediction results
print('R2 Test: RTGP')
pred_rtgp_test = beta0_chain[iter_ELBO] + crossprod(data$img_test,beta_threshold_chain[iter_ELBO,])
(R2_rtgp_test = cor(pred_rtgp_test,data$y_test)^2)

# Evaluation of bias (test) of prediction results
print('Bias Prediction Test: RTGP')
(bias_y_pred_rtgp_test = mean(pred_rtgp_test - data$y_test))

# Evaluation of MSE (test) of prediction results
print('MSE Prediction Test: RTGP')
(mse_y_pred_rtgp_test = mean((pred_rtgp_test - data$y_test)^2))

#####################
### Ridge + LASSO ###
#####################

# Create list to save results in.
result_ridge_lasso = list()

# Start time of model estimation.
time.start = Sys.time()

# Fit Ridge Regression with regularisation hyperparameter identified via 100-fold CV.
lambdas = 10^seq(2, -3, by = -.1)
cv_ridge = cv.glmnet(x = X, y =  y, alpha = 0, lambda = lambdas)
optimal_lambda = cv_ridge$lambda.min
ridge_reg = glmnet(x = X, y = y, lambda = optimal_lambda, alpha = 0, family = 'gaussian')
beta_ridge = as.numeric(unlist(coefficients(ridge_reg)))[-1]
beta0_ridge = as.numeric(unlist(coefficients(ridge_reg)))[1]

# Fit LASSO Regression with regularisation hyperparameter identified via 100-fold CV.
lambdas = 10^seq(2, -3, by = -.1)
cv_lasso = cv.glmnet(x = X, y =  y, alpha = 1, lambda = lambdas)
optimal_lambda_lasso = cv_lasso$lambda.min
lasso_model = glmnet(x = X, y = y, alpha = 1, lambda = optimal_lambda_lasso)
beta_lasso = as.numeric(unlist(coefficients(lasso_model)))[-1]
beta0_lasso = as.numeric(unlist(coefficients(lasso_model)))[1]

# End time of model estimation.
time.end = Sys.time()
time = time.end - time.start

# Concatenate beta coefficients of Ridge and LASSO.
beta_fit = data.frame(LASSO = beta_lasso,
                      Ridge = beta_ridge)

# Save results in list.
result_ridge_lasso$beta_fit = beta_fit
result_ridge_lasso$time = time
save(result_ridge_lasso, file = sprintf("%sresult_ridge_lasso.RData", path_results_LASSO))

# Save beta map from Ridge model.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(beta_fit$Ridge, length(beta_fit$Ridge), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_LASSO, "beta_Ridge.dtseries.nii"), 
  file.path(path_results_LASSO, "beta_Ridge_left.surf.gii"), file.path(path_results_LASSO, "beta_Ridge_cortex_right.surf.gii")
)

# Save beta map from LASSO model.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(beta_fit$LASSO, length(beta_fit$LASSO), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_LASSO, "beta_LASSO.dtseries.nii"), 
  file.path(path_results_LASSO, "beta_LASSO_left.surf.gii"), file.path(path_results_LASSO, "beta_LASSO_cortex_right.surf.gii")
)

# Plot scatterplot between true beta values and estimated beta values from Ridge regression.
df = data.frame(beta_true = beta_true, beta_method = beta_fit$Ridge)
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta True') +
  ylab('Beta Ridge') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
p_plot
ggsave(p_plot, filename=paste0(path_results_LASSO, "scatterplot_beta_true_vs_beta_Ridge.png"), type="cairo-png")

# Plot scatterplot between true beta values and estimated beta values from LASSO regression.
df = data.frame(beta_true = beta_true, beta_method = beta_fit$LASSO)
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta True') +
  ylab('Beta LASSO') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
p_plot
ggsave(p_plot, filename=paste0(path_results_LASSO, "scatterplot_beta_true_vs_beta_LASSO.png"), type="cairo-png")

# Plot scatterplot between estimated beta values from Ridge regression and estimated beta values from LASSO regression.
df = data.frame(beta_true = beta_fit$Ridge, beta_method = beta_fit$LASSO)
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta Ridge') +
  ylab('Beta LASSO') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
p_plot
ggsave(p_plot, filename=paste0(path_results_LASSO, "scatterplot_beta_Ridge_vs_beta_LASSO.png"), type="cairo-png")

# Set true threshold parameter.
threshold_true = 0
# Set mask of true significance.
mask = data$alpha_true

# Evaluate bias of parameter estimates of Ridge model
print('Bias: Ridge')
(bias_ridge_no_effect = mean((beta_fit$Ridge - beta_true)[mask == 0]))
(bias_ridge_effect = mean((beta_fit$Ridge - beta_true)[mask == 1]))
(bias_ridge_total = mean((beta_fit$Ridge - beta_true)))

# Evaluate bias of parameter estimates of LASSO model
print('Bias: LASSO')
(bias_lasso_no_effect = mean((beta_fit$LASSO - beta_true)[mask == 0]))
(bias_lasso_effect = mean((beta_fit$LASSO - beta_true)[mask == 1]))
(bias_lasso_total = mean((beta_fit$LASSO - beta_true)))

# Evaluate MSE of parameter estimates of Ridge model
print('MSE: Ridge')
(mse_ridge_no_effect = mean((beta_fit$Ridge - beta_true)[mask == 0]^2))
(mse_ridge_effect = mean((beta_fit$Ridge - beta_true)[mask == 1]^2))
(mse_ridge_total = mean((beta_fit$Ridge - beta_true)^2))

# Evaluate MSE of parameter estimates of LASSO model
print('MSE: LASSO')
(mse_lasso_no_effect = mean((beta_fit$LASSO - beta_true)[mask == 0]^2))
(mse_lasso_effect = mean((beta_fit$LASSO - beta_true)[mask == 1]^2))
(mse_lasso_total = mean((beta_fit$LASSO - beta_true)^2))

# Evaluate R2 (train) of prediction results of Ridge model
print('R2 Train: Ridge')
pred_ridge = predict(ridge_reg, newx = X)
(R2_ridge = cor(pred_ridge,y)^2)

# Evaluate R2 (train) of prediction results of LASSO model
print('R2 Train: LASSO')
pred_lasso = predict(lasso_model, newx = X)
(R2_lasso = cor(pred_lasso,y)^2)

# Evaluate bias (train) of prediction results of Ridge model
print('Bias Prediction Train: Ridge')
(bias_y_pred_ridge = mean(pred_ridge - y))

# Evaluate bias (train) of prediction results of LASSO model
print('Bias Prediction Train: LASSO')
(bias_y_pred_lasso = mean(pred_lasso - y))

# Evaluate MSE (train) of prediction results of Ridge model
print('MSE Prediction Train: Ridge')
(mse_y_pred_ridge = mean((pred_ridge - y)^2))

# Evaluate MSE (train) of prediction results of LASSO model
print('MSE Prediction Train: LASSO')
(mse_y_pred_lasso = mean((pred_lasso - y)^2))

# Evaluate R2 (test) of prediction results of Ridge model
print('R2 Test: Ridge')
pred_ridge_test = predict(ridge_reg, newx = t(img_test))
(R2_ridge_test = cor(pred_ridge_test,data$y_test)^2)

# Evaluate R2 (test) of prediction results of LASSO model
print('R2 Test: LASSO')
pred_lasso_test = predict(lasso_model, newx = t(img_test))
(R2_lasso_test = cor(pred_lasso_test,data$y_test)^2)

# Evaluate bias (test) of prediction results of Ridge model
print('Bias Prediction Test: RIdge')
(bias_y_pred_ridge_test = mean(pred_ridge_test - data$y_test))

# Evaluate bias (test) of prediction results of LASSO model
print('Bias Prediction Test: LASSO')
(bias_y_pred_lasso_test = mean(pred_lasso_test - data$y_test))

# Evaluate MSE (test) of prediction results of Ridge model
print('MSE Prediction Test: Ridge')
(mse_y_pred_ridge_test = mean((pred_ridge_test - data$y_test)^2))

# Evaluate MSE (test) of prediction results of LASSO model
print('MSE Prediction Test: LASSO')
(mse_y_pred_lasso_test = mean((pred_lasso_test - data$y_test)^2))


######################################
### BR + Normal and BR + Horseshoe ###
######################################

# Create list to save results in.
result_br = list()

# Start time of model estimation.
time.start = Sys.time()

# Fit BR + Normal, BR + Horseshoe, and BR + MFVB.
hs_fit_SOI = fast_horseshoe_lm(y,cbind(1,X), mcmc_sample = 5000, burnin = 2000)
nm_fit_SOI = fast_normal_lm(y,cbind(1,X), mcmc_sample = 5000, burnin = 2000)
mfvb_fit_SOI = fast_mfvb_normal_lm(y,cbind(1,X))

# End time of model estimation.
time.end = Sys.time()
time = time.end - time.start

# Concatenate estimated beta coefficients of models in dataframe.
beta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-1],
                      NM = nm_fit_SOI$post_mean$betacoef[-1],
                      MFVB = mfvb_fit_SOI$post_mean$betacoef[-1])

# Save results in list.
result_br$beta_fit = beta_fit
result_br$time = time
result_br$hs_fit_SOI = hs_fit_SOI
result_br$nm_fit_SOI = nm_fit_SOI
result_br$mfvb_fit_SOI = mfvb_fit_SOI
save(result_br, file = sprintf("%sresult_br.RData", path_results_BR))

# Save estimated beta map of BR + Normal model.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(beta_fit$NM, length(beta_fit$NM), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_BR, "beta_NM.dtseries.nii"), 
  file.path(path_results_BR, "beta_NM_left.surf.gii"), file.path(path_results_BR, "beta_NM_cortex_right.surf.gii")
)

# Save estimated beta map of BR + Horseshoe model.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(beta_fit$HS, length(beta_fit$HS), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_BR, "beta_HS.dtseries.nii"), 
  file.path(path_results_BR, "beta_HS_left.surf.gii"), file.path(path_results_BR, "beta_HS_cortex_right.surf.gii")
)

# Plot scatterplot of true beta values vs estimated beta values of BR + Normal model.
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
ggsave(p_plot, filename=paste0(path_results_BR, "scatterplot_beta_true_vs_beta_NM.png"), type="cairo-png")

# Plot scatterplot of true beta values vs estimated beta values of BR + Horseshoe model.
df = data.frame(beta_true = beta_true, beta_method = beta_fit$HS)
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta True') +
  ylab('Beta HS') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
p_plot
ggsave(p_plot, filename=paste0(path_results_BR, "scatterplot_beta_true_vs_beta_HS.png"), type="cairo-png")

# Plot scatterplot of estimated beta values of BR + Normal model vs estimated beta values of BR + Horseshoe model.
df = data.frame(beta_true = beta_fit$NM, beta_method = beta_fit$HS)
p_plot = ggplot(df, aes(x = beta_true, y = beta_method)) +
  geom_pointdensity(size = 2, adjust = 0.1, show.legend = F) +
  geom_abline(intercept = 0) +
  xlim(min(df), max(df)) +
  ylim(min(df), max(df)) +
  xlab('Beta NM') +
  ylab('Beta HS') +
  theme(text = element_text(size=35)) +
  guides(fill = guide_legend(title=''))
p_plot
ggsave(p_plot, filename=paste0(path_results_BR, "scatterplot_beta_NM_vs_beta_HS.png"), type="cairo-png")

# BR + Normal: Inference evaluation 
# Binary significance: 95%-Equidistant credible interval vertices that do not contain 0
beta_NM = nm_fit_SOI$mcmc$betacoef[-1,]
x=eti(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_NM = ifelse(x2==T,0,1)
# Binary significance: 95%-Highest posterior density credible interval vertices that do not contain 0
x=hdi(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_NM = ifelse(x2==T,0,1)
# Binary significance: Test statistics = posterior beta mean / posterior beta variance -> FDR correction -> Significant vertices determined by p-values lower
# than 5%
t_beta_mean_NM = apply(beta_NM, 2, mean)
t_beta_sd_NM = apply(beta_NM, 2, sd)
t = t_beta_mean_NM / t_beta_sd_NM
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_NM = ifelse(p_fdr < 0.05, 1, 0)

# BR + Horseshoe: Inference evaluation
# Binary significance: 95%-Equidistant credible interval vertices that do not contain 0
beta_HS = hs_fit_SOI$mcmc$betacoef[-1,]
x=eti(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_HS = ifelse(x2==T,0,1)
# Binary significance: 95%-Highest posterior density credible interval vertices that do not contain 0
x=hdi(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_HS = ifelse(x2==T,0,1)
# Binary significance: Test statistics = posterior beta mean / posterior beta variance -> FDR correction -> Significant vertices determined by p-values lower
# than 5%
t_beta_mean_HS = apply(beta_HS, 2, mean)
t_beta_sd_HS = apply(beta_HS, 2, sd)
t = t_beta_mean_HS / t_beta_sd_HS
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_HS = ifelse(p_fdr < 0.05, 1, 0)

# Save binary inference map of BR + Normal model.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(CI_HPDI_binary_NM, length(beta_fit$NM), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_BR, "significance_NM.dtseries.nii"), 
  file.path(path_results_BR, "significance_NM_left.surf.gii"), file.path(path_results_BR, "significance_NM_cortex_right.surf.gii")
)

# Save binary inference map of BR + Horseshoe model.
xii = read_xifti(paste0(path_beta, 'beta_2k.dtseries.nii'))
xii$data$cortex_left = matrix(CI_HPDI_binary_HS, length(beta_fit$HS), 1)
xii$data$cortex_right = matrix(0, 1803, 1)

write_xifti(
  xii, 
  file.path(path_results_BR, "significance_HS.dtseries.nii"), 
  file.path(path_results_BR, "significance_HS_left.surf.gii"), file.path(path_results_BR, "significance_HS_cortex_right.surf.gii")
)

# True significance threshold
threshold_true = 0
# True binary significance mask
mask = data$alpha_true

# Evaluate bias of parameter estimates of BR + Normal model.
print('Bias: NM')
(bias_br_nm_no_effect = mean((beta_fit$NM - beta_true)[mask == 0]))
(bias_br_nm_effect = mean((beta_fit$NM - beta_true)[mask == 1]))
(bias_br_nm_total = mean((beta_fit$NM - beta_true)))

# Evaluate bias of parameter estimates of BR + Horseshoe model.
print('Bias: HS')
(bias_br_hs_no_effect = mean((beta_fit$HS - beta_true)[mask == 0]))
(bias_br_hs_effect = mean((beta_fit$HS - beta_true)[mask == 1]))
(bias_br_hs_total = mean((beta_fit$HS - beta_true)))

# Evaluate MSE of parameter estimates of BR + Normal model.
print('MSE: NM')
(mse_br_nm_no_effect = mean((beta_fit$NM - beta_true)[mask == 0]^2))
(mse_br_nm_effect = mean((beta_fit$NM - beta_true)[mask == 1]^2))
(mse_br_nm_total = mean((beta_fit$NM - beta_true)^2))

# Evaluate MSE of parameter estimates of BR + Horseshoe model.
print('MSE: HS')
(mse_br_hs_no_effect = mean((beta_fit$HS - beta_true)[mask == 0]^2))
(mse_br_hs_effect = mean((beta_fit$HS - beta_true)[mask == 1]^2))
(mse_br_hs_total = mean((beta_fit$HS - beta_true)^2))

# BR + Normal: Inference evaluation 
print('Inference results: NM')
# Binary significance: 95%-Equidistant credible interval vertices that do not contain 0
posterior = nm_fit_SOI$mcmc$betacoef[-1,]
x=eti(data.frame(posterior))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary = ifelse(x2==T,0,1)
# Binary significance: 95%-Highest posterior density credible interval vertices that do not contain 0
x=hdi(data.frame(posterior))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary = ifelse(x2==T,0,1)
# Binary significance: Test statistics = posterior beta mean / posterior beta variance -> FDR correction -> Significant vertices determined by p-values lower
# than 5%
t_beta_mean = apply(nm_fit_SOI$mcmc$betacoef[-1,], 1, mean)
t_beta_sd = apply(nm_fit_SOI$mcmc$betacoef[-1,], 1, sd)
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
(TPR_br_nm = TP/(TP + FN))
(TDR_br_nm = TP/(TP + FP))
(FPR_br_nm = FP/(FP + TN))
(FDR_br_nm = FP/(FP + TP))

TP = CI_ETI_binary + mask
TP = length(TP[TP == 2])
TN = CI_ETI_binary + mask
TN = length(TN[TN == 0])
FP = CI_ETI_binary - mask
FP = length(FP[FP == 1])
FN = CI_ETI_binary - mask
FN = length(FN[FN == -1])
(TPR_CI_ETI_br_nm = TP/(TP + FN))
(TDR_CI_ETI_br_nm = TP/(TP + FP))
(FPR_CI_ETI_br_nm = FP/(FP + TN))
(FDR_CI_ETI_br_nm = FP/(FP + TP))

TP = CI_HPDI_binary + mask
TP = length(TP[TP == 2])
TN = CI_HPDI_binary + mask
TN = length(TN[TN == 0])
FP = CI_HPDI_binary - mask
FP = length(FP[FP == 1])
FN = CI_HPDI_binary - mask
FN = length(FN[FN == -1])
(TPR_CI_HPDI_br_nm = TP/(TP + FN))
(TDR_CI_HPDI_br_nm = TP/(TP + FP))
(FPR_CI_HPDI_br_nm = FP/(FP + TN))
(FDR_CI_HPDI_br_nm = FP/(FP + TP))

# BR + Horseshoe: Inference evaluation 
print('Inference results: HS')
# Binary significance: 95%-Equidistant credible interval vertices that do not contain 0
posterior = hs_fit_SOI$mcmc$betacoef[-1,]
x=eti(data.frame(posterior))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary = ifelse(x2==T,0,1)
# Binary significance: 95%-Highest posterior density credible interval vertices that do not contain 0
x=hdi(data.frame(posterior))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary = ifelse(x2==T,0,1)
# Binary significance: Test statistics = posterior beta mean / posterior beta variance -> FDR correction -> Significant vertices determined by p-values lower
# than 5%
t_beta_mean = apply(hs_fit_SOI$mcmc$betacoef[-1,], 1, mean)
t_beta_sd = apply(hs_fit_SOI$mcmc$betacoef[-1,], 1, sd)
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
(TPR_br_hs = TP/(TP + FN))
(TDR_br_hs = TP/(TP + FP))
(FPR_br_hs = FP/(FP + TN))
(FDR_br_hs = FP/(FP + TP))

TP = CI_ETI_binary + mask
TP = length(TP[TP == 2])
TN = CI_ETI_binary + mask
TN = length(TN[TN == 0])
FP = CI_ETI_binary - mask
FP = length(FP[FP == 1])
FN = CI_ETI_binary - mask
FN = length(FN[FN == -1])
(TPR_CI_ETI_br_hs = TP/(TP + FN))
(TDR_CI_ETI_br_hs = TP/(TP + FP))
(FPR_CI_ETI_br_hs = FP/(FP + TN))
(FDR_CI_ETI_br_hs = FP/(FP + TP))

TP = CI_HPDI_binary + mask
TP = length(TP[TP == 2])
TN = CI_HPDI_binary + mask
TN = length(TN[TN == 0])
FP = CI_HPDI_binary - mask
FP = length(FP[FP == 1])
FN = CI_HPDI_binary - mask
FN = length(FN[FN == -1])
(TPR_CI_HPDI_br_hs = TP/(TP + FN))
(TDR_CI_HPDI_br_hs = TP/(TP + FP))
(FPR_CI_HPDI_br_hs = FP/(FP + TN))
(FDR_CI_HPDI_br_hs = FP/(FP + TP))

# Evaluate R2 (train) of prediction results of BR + Normal model.
print('R2 Train: NM')
pred_br_nm = nm_fit_SOI$post_mean$betacoef[1] + crossprod(t(X),beta_fit$NM)
(R2_br_nm = cor(pred_br_nm,y)^2)

# Evaluate R2 (train) of prediction results of BR + Horseshoe model.
print('R2 Train: HS')
pred_br_hs = hs_fit_SOI$post_mean$betacoef[1] + crossprod(t(X),beta_fit$HS)
(R2_br_hs = cor(pred_br_hs,y)^2)

# Evaluate bias (train) of prediction results of BR + Normal model.
print('Bias Prediction Train: NM')
(bias_y_pred_br_nm = mean(pred_br_nm - y))

# Evaluate bias (train) of prediction results of BR + Horseshoe model.
print('Bias Prediction Train: HS')
(bias_y_pred_br_hs = mean(pred_br_hs - y))

# Evaluate MSE (train) of prediction results of BR + Normal model.
print('MSE Prediction Train: NM')
(mse_y_pred_br_nm = mean((pred_br_nm - y)^2))

# Evaluate MSE (train) of prediction results of BR + Horseshoe model.
print('MSE Prediction Train: HS')
(mse_y_pred_br_hs = mean((pred_br_hs - y)^2))

# Evaluate R2 (test) of prediction results of BR + Normal model.
print('R2 Test: NM')
pred_br_nm_test = nm_fit_SOI$post_mean$betacoef[1] + crossprod(data$img_test,beta_fit$NM)
(R2_br_nm_test = cor(pred_br_nm_test,data$y_test)^2)

# Evaluate R2 (test) of prediction results of BR + Horseshoe model.
print('R2 Test: HS')
pred_br_hs_test = hs_fit_SOI$post_mean$betacoef[1] + crossprod(data$img_test,beta_fit$HS)
(R2_br_hs_test = cor(pred_br_hs_test,data$y_test)^2)

# Evaluate bias (test) of prediction results of BR + Normal model.
print('Bias Prediction Test: NM')
(bias_y_pred_br_nm_test = mean(pred_br_nm_test - data$y_test))

# Evaluate bias (test) of prediction results of BR + Horseshoe model.
print('Bias Prediction Test: HS')
(bias_y_pred_br_hs_test = mean(pred_br_hs_test - data$y_test))

# Evaluate MSE (test) of prediction results of BR + Normal model.
print('MSE Prediction Test: NM')
(mse_y_pred_br_nm_test = mean((pred_br_nm_test - data$y_test)^2))

# Evaluate MSE (test) of prediction results of BR + Horseshoe model.
print('MSE Prediction Test: HS')
(mse_y_pred_br_hs_test = mean((pred_br_hs_test - data$y_test)^2))


##############
### Tables ###
##############

row_name = c('RTGP', 'GPR + Normal', 'GPR + Horseshoe', 'BR + Normal', 'BR + Horseshoe', 'Ridge', 'LASSO')
row_name_inference = c('RTGP', 'GPR + Normal', 'GPR + Horseshoe', 'BR + Normal', 'BR + Horseshoe')

column_name_inference = c('TPR', 'TDR', 'FPR', 'FDR')
column_name_parameters = c('Bias', 'MSE')
column_name_prediction = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')

# Table summarising the inference results
matrix_inference = matrix(NA, 5, 4)
matrix_inference[1,1] = TPR_rtgp
matrix_inference[1,2] = TDR_rtgp
matrix_inference[1,3] = FPR_rtgp
matrix_inference[1,4] = FDR_rtgp
matrix_inference[2,1] = TPR_CI_HPDI_nm
matrix_inference[2,2] = TDR_CI_HPDI_nm
matrix_inference[2,3] = FPR_CI_HPDI_nm
matrix_inference[2,4] = FDR_CI_HPDI_nm
matrix_inference[3,1] = TPR_CI_HPDI_hs
matrix_inference[3,2] = TDR_CI_HPDI_hs
matrix_inference[3,3] = FPR_CI_HPDI_hs
matrix_inference[3,4] = FDR_CI_HPDI_hs
matrix_inference[4,1] = TPR_CI_HPDI_br_nm
matrix_inference[4,2] = TDR_CI_HPDI_br_nm
matrix_inference[4,3] = FPR_CI_HPDI_br_nm
matrix_inference[4,4] = FDR_CI_HPDI_br_nm
matrix_inference[5,1] = TPR_CI_HPDI_br_hs
matrix_inference[5,2] = TDR_CI_HPDI_br_hs
matrix_inference[5,3] = FPR_CI_HPDI_br_hs
matrix_inference[5,4] = FDR_CI_HPDI_br_hs

matrix_inference = data.frame(matrix_inference)
colnames(matrix_inference) = column_name_inference
rownames(matrix_inference) = row_name_inference
xtable(matrix_inference, digits = 4)
print(xtable(matrix_inference, digits = 4), file=paste0(path_results_all, "table_inference.txt"))

# Table summarising the parameter estimate results
matrix_parameters = matrix(NA, 7, 2)

matrix_parameters[1,1] = bias_rtgp_total
matrix_parameters[1,2] = mse_rtgp_total

matrix_parameters[2,1] = bias_nm_total
matrix_parameters[2,2] = mse_nm_total

matrix_parameters[3,1] = bias_hs_total
matrix_parameters[3,2] = mse_hs_total

matrix_parameters[4,1] = bias_br_nm_total
matrix_parameters[4,2] = mse_br_nm_total

matrix_parameters[5,1] = bias_br_hs_total
matrix_parameters[5,2] = mse_br_hs_total

matrix_parameters[6,1] = bias_ridge_total
matrix_parameters[6,2] = mse_ridge_total

matrix_parameters[7,1] = bias_lasso_total
matrix_parameters[7,2] = mse_lasso_total

colnames(matrix_parameters) = column_name_parameters
rownames(matrix_parameters) = row_name
xtable(matrix_parameters, digits = 4)
print(xtable(matrix_parameters, digits = 4), file=paste0(path_results_all, "table_parameters.txt"))

# Table summarising the prediction results
matrix_prediction = matrix(NA, 7, 4)
matrix_prediction[1,1] = R2_rtgp
matrix_prediction[1,2] = mse_y_pred_rtgp
matrix_prediction[1,3] = R2_rtgp_test
matrix_prediction[1,4] = mse_y_pred_rtgp_test

matrix_prediction[2,1] = R2_nm
matrix_prediction[2,2] = mse_y_pred_nm
matrix_prediction[2,3] = R2_nm_test
matrix_prediction[2,4] = mse_y_pred_nm_test

matrix_prediction[3,1] = R2_hs
matrix_prediction[3,2] = mse_y_pred_hs
matrix_prediction[3,3] = R2_hs_test
matrix_prediction[3,4] = mse_y_pred_hs_test

matrix_prediction[4,1] = R2_br_nm
matrix_prediction[4,2] = mse_y_pred_br_nm
matrix_prediction[4,3] = R2_br_nm_test
matrix_prediction[4,4] = mse_y_pred_br_nm_test

matrix_prediction[5,1] = R2_br_hs
matrix_prediction[5,2] = mse_y_pred_br_hs
matrix_prediction[5,3] = R2_br_hs_test
matrix_prediction[5,4] = mse_y_pred_br_hs_test

matrix_prediction[6,1] = R2_ridge
matrix_prediction[6,2] = mse_y_pred_ridge
matrix_prediction[6,3] = R2_ridge_test
matrix_prediction[6,4] = mse_y_pred_ridge_test

matrix_prediction[7,1] = R2_lasso
matrix_prediction[7,2] = mse_y_pred_lasso
matrix_prediction[7,3] = R2_lasso_test
matrix_prediction[7,4] = mse_y_pred_lasso_test

colnames(matrix_prediction) = column_name_prediction
rownames(matrix_prediction) = row_name
xtable(matrix_prediction, digits = 4)
print(xtable(matrix_prediction, digits = 4), file=paste0(path_results_all, "table_prediction.txt"))

# Save all three tables and concatenated results.
result_matrix = list()
result_matrix$matrix_inference = matrix_inference
result_matrix$matrix_parameters = matrix_parameters
result_matrix$matrix_prediction = matrix_prediction
save(result_matrix, file = paste0(path_results_all, 'result_matrix.RData'))

