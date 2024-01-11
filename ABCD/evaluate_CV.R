# Script: Evaluate results of 10-fold cross validation for ABCD data.

#############
### Paths ###
#############

# Path to folder that contains Human connectome workbench software.
wb_path = ""

# Path to results of RTGP-VI model results of left hemisphere (phi = 0.2, nu = 1.2).
path_L = ''
# Path to results of RTGP-VI model results of right hemisphere (phi = 0.05, nu = 1.2).
path_R = ''
# Path to folder where data is stored.
path_data = ""
# Path to folder where to save results in.
path_out = ""

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
library(data.table)
library(ciftiTools)
options(bitmapType='cairo-png')
ciftiTools.setOption('wb_path', wb_path)

############
### Data ###
############

# Data matrix containing confounding variables of training dataset.
X_confounds_all = fread(paste0(path_data, 'X_confounds_train.csv'))
X_confounds_all = data.matrix(X_confounds_all)[,2:ncol(X_confounds_all)]

# Data matrix containing test statistics maps of left hemisphere of training dataset.
X_left_all = fread(paste0(path_data, 'X_left_cortex_train.csv'))
X_left_all = data.matrix(X_left_all)[,2:ncol(X_left_all)]

# Data matrix containing test statistics maps of right hemisphere of training dataset.
X_right_all = fread(paste0(path_data, 'X_right_cortex_train.csv'))
X_right_all = data.matrix(X_right_all)[,2:ncol(X_right_all)]

# Data vector containing intelligence scores of training dataset.
y_all = as.numeric(unlist(fread(paste0(path_data, 'y_train.csv'))[,2]))

# Create vectors to store prediction evaluations (R2 (train, test) of CV partitions in.
r2_left_train = numeric(10)
r2_right_train = numeric(10)
r2_left_test = numeric(10)
r2_right_test = numeric(10)

# Create vectors to store prediction evaluations (MSE (train, test) of CV partitions in.
mse_left_train = numeric(10)
mse_right_train = numeric(10)
mse_left_test = numeric(10)
mse_right_test = numeric(10)

# Evaluate predictive results of each of the 10-fold CV partitions.
for(i in 1:10){
  
# Load results of left hemisphere. 
load(paste0(path_L, toString(i), '/result_rtgp_vi_left_cortex.RData'))
result_left = result
# Read in CV split in training and test data.
idx = read.csv(paste0(path_L, toString(i), '/idx.csv'))[,2]

# Load results of right hemisphere. 
load(paste0(path_R, toString(i), '/result_rtgp_vi_right_cortex.RData'))
result_right = result
rm(result)
gc()

# Determine highest ELBO across iterations of left hemisphere.
iter_ELBO_left = which(result_left$ELBO == max(result_left$ELBO))[1]
# Determine highest ELBO across iterations of right hemisphere.
iter_ELBO_right = which(result_right$ELBO == max(result_right$ELBO))[1]

# Create dataset according to CV split.
# Confounding data
X_confounds = X_confounds_all[-idx,]
X_confounds_test = X_confounds_all[idx,]
# Image data of left hemisphere
X_left = X_left_all[-idx,]
X_left_test = X_left_all[idx,]
# Image data of right hemisphere
X_right = X_right_all[-idx,]
X_right_test = X_right_all[idx,]
# Output: Intelligence scores
y = y_all[-idx]
y_test = y_all[idx]

# Evaluate predictive results of left hemisphere.
pred_left_rtgp = result_left$beta0[iter_ELBO_left] + X_confounds %*% result_left$beta_confound[iter_ELBO_left,] + crossprod(t(X_left),result_left$beta_threshold[iter_ELBO_left,])
r2_left_train[i] = cor(pred_left_rtgp,y)^2
mse_left_train[i] = mean((pred_left_rtgp - y)^2)

pred_left_rtgp_test = result_left$beta0[iter_ELBO_left] + X_confounds_test %*% result_left$beta_confound[iter_ELBO_left, ] + crossprod(t(X_left_test),result_left$beta_thresholded[iter_ELBO_left,])
r2_left_test[i] = cor(pred_left_rtgp_test, y_test)^2
mse_left_test[i] = mean((pred_left_rtgp_test - y_test)^2)

# Evaluate predictive results of right hemisphere.
pred_right_rtgp = result_right$beta0[iter_ELBO_right] + X_confounds %*% result_right$beta_confound[iter_ELBO_right,] + crossprod(t(X_right),result_right$beta_threshold[iter_ELBO_right,])
r2_right_train[i] = cor(pred_right_rtgp,y)^2
mse_right_train[i] = mean((pred_right_rtgp - y)^2)

pred_right_rtgp_test = result_right$beta0[iter_ELBO_right] + X_confounds_test %*% result_right$beta_confound[iter_ELBO_right, ] + crossprod(t(X_right_test),result_right$beta_thresholded[iter_ELBO_right,])
r2_right_test[i] = cor(pred_right_rtgp_test, y_test)^2
mse_right_test[i] = mean((pred_right_rtgp_test - y_test)^2)

}

# Concatenate and save results in tables.
matrix_results = matrix(NA, 4, 4)
row_name = c('Mean (left)', 'Sd (left)', 'Mean (right)', 'Sd (right)')
col_name = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')

matrix_results[1,1] = mean(r2_left_train, na.rm = T)
matrix_results[2,1] = sd(r2_left_train, na.rm = T)

matrix_results[1,2] = mean(mse_left_train, na.rm = T)
matrix_results[2,2] = sd(mse_left_train, na.rm = T)

matrix_results[1,3] = mean(r2_left_test, na.rm = T)
matrix_results[2,3] = sd(r2_left_test, na.rm = T)

matrix_results[1,4] = mean(mse_left_test, na.rm = T)
matrix_results[2,4] = sd(mse_left_test, na.rm = T)

matrix_results[3,1] = mean(r2_right_train, na.rm = T)
matrix_results[4,1] = sd(r2_right_train, na.rm = T)

matrix_results[3,2] = mean(mse_right_train, na.rm = T)
matrix_results[4,2] = sd(mse_right_train, na.rm = T)

matrix_results[3,3] = mean(r2_right_test, na.rm = T)
matrix_results[4,3] = sd(r2_right_test, na.rm = T)

matrix_results[3,4] = mean(mse_right_test, na.rm = T)
matrix_results[4,4] = sd(mse_right_test, na.rm = T)


colnames(matrix_results) = col_name
rownames(matrix_results) = row_name
xtable(matrix_results, digits = 4)
print(xtable(matrix_results, digits = 4), file=paste0(path_out, "table_cv.txt"))

write.csv(matrix_results, paste0(path_out, 'result_matrix.csv'))


