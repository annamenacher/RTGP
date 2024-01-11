# Script: File to estimate Bayesian regression models with all vertices (BR + Normal and BR + Horseshoe model) with ABCD data.

#############
### Paths ###
#############

# Path to folder where Human Connectome workbench software is stored.
wb_path = ''

# Path to folder where to store results in.
path_core = ""
path_results = path_results_NM = path = paste0(path_core, "/BR/")
dir.create(file.path(path_results), showWarnings = FALSE)

# Path to folder where data is stored in.
path_data = ''

# Path to file with example cifti file (e.g. COPE map of a random subject of ABCD study)
path_xii_image = ''

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

options(bitmapType='cairo-png')
library(ciftiTools)
ciftiTools.setOption('wb_path', wb_path)

############
### Data ###
############

### Left hemisphere ###

# Read in data with confounding variables (training data).
X_confounds_train = fread(paste0(path_data, 'X_confounds_train.csv'))
X_confounds_train = data.matrix(X_confounds_train)[,2:ncol(X_confounds_train)]
# Read in imaging data of left hemisphere (training data).
X_cortex_left_train = fread(paste0(path_data, 'X_left_cortex_train.csv'))
X_cortex_left_train = data.matrix(X_cortex_left_train)[,2:ncol(X_cortex_left_train)]
# Number of vertices
M = ncol(X_cortex_left_train)
# Output vector: intelligence scores (training data)
y_train = as.numeric(unlist(fread(paste0(path_data, 'y_train.csv'))[,2]))

# Read in data with confounding variables (validation data).
X_confounds_test = fread(paste0(path_data, "X_confounds_val.csv"))
X_confounds_test = X_confounds_test[,2:33]
X_confounds_test = data.matrix(X_confounds_test)

# Output vector: intelligence scores (validation data)
y_test = read.csv("/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/not_centered/y_val.csv")[,2]

# Read in imaging data of left hemisphere (validation data).
X_left_test = fread("/well/nichols/users/fxo071/IOI-GP/ABCD/data/test_statistics/not_centered/X_left_cortex_val.csv")
X_left_test = X_left_test[,2:(M+1)]
X_left_test = data.matrix(X_left_test)

# Number of confounding variables
P = ncol(X_confounds_train)

################
### Analysis ###
################

# Create lists to store results in.
result_time = list()
result_nm = list()

# Start time of analysis.
time.start = Sys.time()

# Fit BR + Normal, BR + Horseshoe, and BR + MFVB models.
hs_fit_SOI = fast_horseshoe_lm(y_train,cbind(1,X_confounds_train, X_cortex_left_train), mcmc_sample = 5000, burnin = 2000)
nm_fit_SOI = fast_normal_lm(y_train,cbind(1,X_confounds_train, X_cortex_left_train), mcmc_sample = 5000, burnin = 2000)
mfvb_fit_SOI = fast_mfvb_normal_lm(y_train,cbind(1,X_confounds_train, X_cortex_left_train))

# End time of analysis.
time.end = Sys.time()
time = time.end - time.start

# Save time to estimate models.
result_time$time_nm = time

# Concatenate beta values of models in dataframe.
beta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-c(1:(P+1))],
                       NM = nm_fit_SOI$post_mean$betacoef[-c(1:(P+1))],
                       MFVB = mfvb_fit_SOI$post_mean$betacoef[-c(1:(P+1))])

# Save results in list.
beta_fit_left = beta_fit
result_nm$beta_fit = beta_fit
result_nm$hs_fit_SOI = hs_fit_SOI
result_nm$nm_fit_SOI = nm_fit_SOI
result_nm$mfvb_fit_SOI = mfvb_fit_SOI
result_nm$time = time

save(result_nm, file = sprintf("%sresult_nm_left.RData", path_results_NM))

# Plot scatterplot between BR + Normal beta values and BR + Horseshoe beta values.
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
ggsave(p_plot, filename=paste0(path_results_NM, "scatterplot_beta_NM_vs_beta_HS_left.png"), type="cairo-png")


# Inference results: BR + Normal
# Determine binary significance results via 95% ETI CI (significant if 0 is not included in CI, not significant if 0 is included in CI)
beta_NM = nm_fit_SOI$mcmc$betacoef[-c(1:(P+1)),]
x=eti(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_NM_left = ifelse(x2==T,0,1)
# Determine binary significance results via 95% HPDI CI (significant if 0 is not included in CI, not significant if 0 is included in CI)
x=hdi(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_NM_left = ifelse(x2==T,0,1)
# Determine binary significance results via test statistics (FDR correct p-values at 5%)
t_beta_mean_NM = apply(beta_NM, 2, mean)
t_beta_sd_NM = apply(beta_NM, 2, sd)
t = t_beta_mean_NM / t_beta_sd_NM
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_NM_left = ifelse(p_fdr < 0.05, 1, 0)

# Inference results: BR + Horseshoe
# Determine binary significance results via 95% ETI CI (significant if 0 is not included in CI, not significant if 0 is included in CI)
beta_HS = hs_fit_SOI$mcmc$betacoef[-c(1:(P+1)),]
x=eti(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_HS_left = ifelse(x2==T,0,1)
# Determine binary significance results via 95% HPDI CI (significant if 0 is not included in CI, not significant if 0 is included in CI)
x=hdi(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_HS_left = ifelse(x2==T,0,1)
# Determine binary significance results via test statistics (FDR correct p-values at 5%)
t_beta_mean_HS = apply(beta_HS, 2, mean)
t_beta_sd_HS = apply(beta_HS, 2, sd)
t = t_beta_mean_HS / t_beta_sd_HS
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_HS_left = ifelse(p_fdr < 0.05, 1, 0)


# Evaluation of prediction results
# Evaluate R2 (train) of BR + Normal model.
print('R2 Train: BR + Normal')
pred_nm = nm_fit_SOI$post_mean$betacoef[1] + crossprod(t(X_cortex_left_train),beta_fit$NM) + X_confounds_train%*%nm_fit_SOI$post_mean$betacoef[2:(P+1)]
(R2_nm = cor(pred_nm,y_train)^2)

# Evaluate R2 (train) of BR + Horseshoe model.
print('R2 Train: BR + Horseshoe')
pred_hs = hs_fit_SOI$post_mean$betacoef[1] + crossprod(t(X_cortex_left_train),beta_fit$HS) + X_confounds_train%*%hs_fit_SOI$post_mean$betacoef[2:(P+1)]
(R2_hs = cor(pred_hs,y_train)^2)

# Evaluate bias (train) of BR + Normal model.
print('Bias Prediction Train: BR + Normal')
(bias_y_pred_nm = mean(pred_nm - y_train))

# Evaluate bias (train) of BR + Horseshoe model.
print('Bias Prediction Train: BR + Horseshoe')
(bias_y_pred_hs = mean(pred_hs - y_train))

# Evaluate MSE (train) of BR + Normal model.
print('MSE Prediction Train: BR + Normal')
(mse_y_pred_nm = mean((pred_nm - y_train)^2))

# Evaluate MSE (train) of BR + Horseshoe model.
print('MSE Prediction Train: BR + Horseshoe')
(mse_y_pred_hs = mean((pred_hs - y_train)^2))

# Evaluate R2 (test) of BR + Normal model.
print('R2 Test: BR + Normal')
pred_nm_test = nm_fit_SOI$post_mean$betacoef[1] + crossprod(t(X_left_test),beta_fit$NM) + X_confounds_test%*%nm_fit_SOI$post_mean$betacoef[2:(P+1)]
(R2_nm_test = cor(pred_nm_test,y_test)^2)

# Evaluate R2 (test) of BR + Horseshoe model.
print('R2 Test: BR + Horseshoe')
pred_hs_test = hs_fit_SOI$post_mean$betacoef[1] + crossprod(t(X_left_test),beta_fit$HS) + X_confounds_test%*%hs_fit_SOI$post_mean$betacoef[2:(P+1)]
(R2_hs_test = cor(pred_hs_test,y_test)^2)

# Evaluate bias (test) of BR + Normal model.
print('Bias Prediction Test: BR + Normal')
(bias_y_pred_nm_test = mean(pred_nm_test - y_test))

# Evaluate bias (test) of BR + Horseshoe model.
print('Bias Prediction Test: BR + Horseshoe')
(bias_y_pred_hs_test = mean(pred_hs_test - y_test))

# Evaluate MSE (test) of BR + Normal model.
print('MSE Prediction Test: BR + Normal')
(mse_y_pred_nm_test = mean((pred_nm_test - y_test)^2))

# Evaluate MSE (test) of BR + Horseshoe model.
print('MSE Prediction Test: BR + Horseshoe')
(mse_y_pred_hs_test = mean((pred_hs_test - y_test)^2))

##############
### Tables ###
##############

# Save prediction results in tables.
row_names = c('NM', 'HS')
column_name_prediction = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')

matrix_prediction = matrix(NA, 2, 4)
matrix_prediction[1,1] = R2_nm
matrix_prediction[1,2] = mse_y_pred_nm
matrix_prediction[1,3] = R2_nm_test
matrix_prediction[1,4] = mse_y_pred_nm_test

matrix_prediction[2,1] = R2_hs
matrix_prediction[2,2] = mse_y_pred_hs
matrix_prediction[2,3] = R2_hs_test
matrix_prediction[2,4] = mse_y_pred_hs_test

colnames(matrix_prediction) = column_name_prediction
rownames(matrix_prediction) = row_names
xtable(matrix_prediction, digits = 4)
print(xtable(matrix_prediction, digits = 4), file=paste0(path, "table_prediction_left.txt"))

result_matrix = list()
result_matrix$matrix_prediction = matrix_prediction
save(result_matrix, file = paste0(path, 'result_matrix_left.RData'))


### Right hemisphere ###

# Read in imaging data of right hemisphere (training data).
X_cortex_right_train = fread(paste0(path_data, 'X_right_cortex_train.csv'))
X_cortex_right_train = data.matrix(X_cortex_right_train)[,2:ncol(X_cortex_right_train)]
# Number of vertices 
M = ncol(X_cortex_right_train)

# Read in imaging data of right hemisphere (validation data).
X_right_test = fread(paste0(path_data, "X_right_cortex_val.csv"))
X_right_test = X_right_test[,2:(M+1)]
X_right_test = data.matrix(X_right_test)

# Create lists to save results in.
result_time = list()
result_nm = list()

# Start time of analysis
time.start = Sys.time()

# Fit BR + Normal, BR + Horseshoe, and BR + MFVB models.
hs_fit_SOI = fast_horseshoe_lm(y_train,cbind(1,X_confounds_train, X_cortex_right_train), mcmc_sample = 5000, burnin = 2000)
nm_fit_SOI = fast_normal_lm(y_train,cbind(1,X_confounds_train, X_cortex_right_train), mcmc_sample = 5000, burnin = 2000)
mfvb_fit_SOI = fast_mfvb_normal_lm(y_train,cbind(1,X_confounds_train, X_cortex_right_train))

# End time of analysis.
time.end = Sys.time()
time = time.end - time.start

# Save time needed for analysis.
result_time$time_nm = time

# Concatenate beta values of models in dataframe.
beta_fit = data.frame(HS = hs_fit_SOI$post_mean$betacoef[-c(1:(P+1))],
                       NM = nm_fit_SOI$post_mean$betacoef[-c(1:(P+1))],
                       MFVB = mfvb_fit_SOI$post_mean$betacoef[-c(1:(P+1))])

# Save results of models of right hemisphere in list.
beta_fit_right = beta_fit
result_nm$beta_fit = beta_fit
result_nm$hs_fit_SOI = hs_fit_SOI
result_nm$nm_fit_SOI = nm_fit_SOI
result_nm$mfvb_fit_SOI = mfvb_fit_SOI
result_nm$time = time
save(result_nm, file = sprintf("%sresult_nm_right.RData", path_results_NM))

# Plot scatterplot between estimated beta values of BR + Normal model and BR + Horseshoe model.
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
ggsave(p_plot, filename=paste0(path_results_NM, "scatterplot_beta_NM_vs_beta_HS_right.png"), type="cairo-png")


# Inference results: BR + Normal
# Determine binary significance results via 95% ETI CI (significant if 0 is not included in CI, not significant if 0 is included in CI)
beta_NM = nm_fit_SOI$mcmc$betacoef[-c(1:(P+1)),]
x=eti(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_NM_right = ifelse(x2==T,0,1)
# Determine binary significance results via 95% HPDI CI (significant if 0 is not included in CI, not significant if 0 is included in CI)
x=hdi(data.frame(beta_NM))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_NM_right = ifelse(x2==T,0,1)
# Determine binary significance results via test statistics (FDR correct p-values at 5%)
t_beta_mean_NM = apply(beta_NM, 2, mean)
t_beta_sd_NM = apply(beta_NM, 2, sd)
t = t_beta_mean_NM / t_beta_sd_NM
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_NM_right = ifelse(p_fdr < 0.05, 1, 0)

# Inference results: BR + Horseshoe
# Determine binary significance results via 95% ETI CI (significant if 0 is not included in CI, not significant if 0 is included in CI)
beta_HS = hs_fit_SOI$mcmc$betacoef[-c(1:(P+1)),]
x=eti(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_ETI_binary_HS_right = ifelse(x2==T,0,1)
# Determine binary significance results via 95% HPDI CI (significant if 0 is not included in CI, not significant if 0 is included in CI)
x=hdi(data.frame(beta_HS))
x2 = apply(matrix(1:M,M,1),1,function(j) between(0,x[j,3],x[j,4]))
CI_HPDI_binary_HS_right = ifelse(x2==T,0,1)
# Determine binary significance results via test statistics (FDR correct p-values at 5%)
t_beta_mean_HS = apply(beta_HS, 2, mean)
t_beta_sd_HS = apply(beta_HS, 2, sd)
t = t_beta_mean_HS / t_beta_sd_HS
p_val = 2*pnorm(q=t, lower.tail=FALSE)
p_fdr = p.adjust(p_val, method = 'fdr', n = length(p_val))
t_bin_fdr_HS_right = ifelse(p_fdr < 0.05, 1, 0)

# Evaluation of prediction results
# Evaluate R2 (train) of BR + Normal model.
print('R2 Train: BR + Normal')
pred_nm = nm_fit_SOI$post_mean$betacoef[1] + crossprod(t(X_cortex_right_train),beta_fit$NM) + X_confounds_train%*%nm_fit_SOI$post_mean$betacoef[2:(P+1)]
(R2_nm = cor(pred_nm,y_train)^2)

# Evaluate R2 (train) of BR + Horseshoe model.
print('R2 Train: BR + Horseshoe')
pred_hs = hs_fit_SOI$post_mean$betacoef[1] + crossprod(t(X_cortex_right_train),beta_fit$HS) + X_confounds_train%*%hs_fit_SOI$post_mean$betacoef[2:(P+1)]
(R2_hs = cor(pred_hs,y_train)^2)

# Evaluate bias (train) of BR + Normal model.
print('Bias Prediction Train: BR + Normal')
(bias_y_pred_nm = mean(pred_nm - y_train))

# Evaluate R2 (train) of BR + Horseshoe model.
print('Bias Prediction Train: BR + Horseshoe')
(bias_y_pred_hs = mean(pred_hs - y_train))

# Evaluate MSE (train) of BR + Normal model.
print('MSE Prediction Train: BR + Normal')
(mse_y_pred_nm = mean((pred_nm - y_train)^2))

# Evaluate MSE (train) of BR + Horseshoe model.
print('MSE Prediction Train: BR + Horseshoe')
(mse_y_pred_hs = mean((pred_hs - y_train)^2))

# Evaluate R2 (test) of BR + Normal model.
print('R2 Test: BR + Normal')
pred_nm_test = nm_fit_SOI$post_mean$betacoef[1] + crossprod(t(X_right_test),beta_fit$NM) + X_confounds_test%*%nm_fit_SOI$post_mean$betacoef[2:(P+1)]
(R2_nm_test = cor(pred_nm_test,y_test)^2)

# Evaluate R2 (test) of BR + Horseshoe model.
print('R2 Test: BR + Horseshoe')
pred_hs_test = hs_fit_SOI$post_mean$betacoef[1] + crossprod(t(X_right_test),beta_fit$HS) + X_confounds_test%*%hs_fit_SOI$post_mean$betacoef[2:(P+1)]
(R2_hs_test = cor(pred_hs_test,y_test)^2)

# Evaluate bias (test) of BR + Normal model.
print('Bias Prediction Test: BR + Normal')
(bias_y_pred_nm_test = mean(pred_nm_test - y_test))

# Evaluate bias (test) of BR + Horseshoe model.
print('Bias Prediction Test: BR + Horseshoe')
(bias_y_pred_hs_test = mean(pred_hs_test - y_test))

# Evaluate MSE (test) of BR + Normal model.
print('MSE Prediction Test: BR + Normal')
(mse_y_pred_nm_test = mean((pred_nm_test - y_test)^2))

# Evaluate MSE (test) of BR + Horseshoe model.
print('MSE Prediction Test: BR + Horseshoe')
(mse_y_pred_hs_test = mean((pred_hs_test - y_test)^2))

##############
### Tables ###
##############

# Store predictive results in tables.
row_names = c('NM', 'HS')
column_name_prediction = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')

matrix_prediction = matrix(NA, 2, 4)
matrix_prediction[1,1] = R2_nm
matrix_prediction[1,2] = mse_y_pred_nm
matrix_prediction[1,3] = R2_nm_test
matrix_prediction[1,4] = mse_y_pred_nm_test

matrix_prediction[2,1] = R2_hs
matrix_prediction[2,2] = mse_y_pred_hs
matrix_prediction[2,3] = R2_hs_test
matrix_prediction[2,4] = mse_y_pred_hs_test

colnames(matrix_prediction) = column_name_prediction
rownames(matrix_prediction) = row_names
xtable(matrix_prediction, digits = 4)
print(xtable(matrix_prediction, digits = 4), file=paste0(path, "table_prediction_right.txt"))

result_matrix = list()
result_matrix$matrix_prediction = matrix_prediction
save(result_matrix, file = paste0(path, 'result_matrix_right.RData'))

### Visualisation ###

# Save beta map of BR + Normal model.
xii = read_xifti(path_xii_image)
xii$data$cortex_left = matrix(beta_fit_left$NM, length(beta_fit_left$NM), 1)
xii$data$cortex_right = matrix(beta_fit_right$NM, length(beta_fit_right$NM), 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "beta_NM.dtseries.nii"), 
  file.path(path_results_NM, "beta_NM_left.surf.gii"), file.path(path_results_NM, "beta_NM_cortex_right.surf.gii")
)

# Save beta map of BR + Horseshoe model.
xii = read_xifti(path_xii_image)
xii$data$cortex_left = matrix(beta_fit_left$HS, length(beta_fit_left$HS), 1)
xii$data$cortex_right = matrix(beta_fit_right$HS, length(beta_fit_right$HS), 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "beta_HS.dtseries.nii"), 
  file.path(path_results_NM, "beta_HS_left.surf.gii"), file.path(path_results_NM, "beta_HS_cortex_right.surf.gii")
)


# Save binary significance map (evaluated via HPDI) of BR + Normal model.
xii = read_xifti(path_xii_image)
xii$data$cortex_left = matrix(CI_HPDI_binary_NM_left, 29696, 1)
xii$data$cortex_right = matrix(CI_HPDI_binary_NM_right, 29716, 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "significance_NM.dtseries.nii"), 
  file.path(path_results_NM, "significance_NM_left.surf.gii"), file.path(path_results_NM, "significance_NM_cortex_right.surf.gii")
)

# Save binary significance map (evaluated via HPDI) of BR + Horseshoe model.
xii = read_xifti(path_xii_image)
xii$data$cortex_left = matrix(CI_HPDI_binary_HS_left, 29696, 1)
xii$data$cortex_right = matrix(CI_HPDI_binary_HS_right, 29716, 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "significance_HS.dtseries.nii"), 
  file.path(path_results_NM, "significance_HS_left.surf.gii"), file.path(path_results_NM, "significance_HS_cortex_right.surf.gii")
)
