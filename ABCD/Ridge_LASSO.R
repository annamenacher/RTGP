# Script: File to estimate Ridge regression and LASSO models with ABCD data.

#############
### Paths ###
#############

# Path to folder which contains Human connectome workbench software.
wb_path = ''

# Path to folder where to store results in.
path_core = ""
path_results = path_results_NM = path = paste0(path_core, "/glmnet/")
dir.create(file.path(path_results), showWarnings = FALSE)

# Path to folder where data is stored in.
path_data = ''

# Path to file with example cifti file (e.g. COPE cifti file of a random subject of ABCD study)
path_xii_image = ""

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
library(glmnet)

options(bitmapType='cairo-png')
library(ciftiTools)
ciftiTools.setOption('wb_path', wb_path)

################
### Analysis ###
################

### Left hemisphere ###

# Read in data matrix which contains confounding variables (training data)
X_confounds_train = fread(paste0(path_data, 'X_confounds_train.csv'))
X_confounds_train = data.matrix(X_confounds_train)[,2:ncol(X_confounds_train)]
# Read in imaging data of left hemisphere (training data)
X_cortex_left_train = fread(paste0(path_data, 'X_left_cortex_train.csv'))
X_cortex_left_train = data.matrix(X_cortex_left_train)[,2:ncol(X_cortex_left_train)]
# Number of vertices 
M = ncol(X_cortex_left_train)
# Output vector: Intelligence scores (training data)
y_train = as.numeric(unlist(fread(paste0(path_data, 'y_train.csv'))[,2]))

# Read in data matrix which contains confounding variables (validation data)
X_confounds_test = fread(paste0(path_data, "X_confounds_val.csv"))
X_confounds_test = X_confounds_test[,2:33]
X_confounds_test = data.matrix(X_confounds_test)

# Output vector: Intelligence scores (validation data)
y_test = read.csv(paste0(path_data, "y_val.csv"))[,2]

# Read in imaging data of left hemisphere (validation data)
X_left_test = fread(paste0(path_data, "X_left_cortex_val.csv"))
X_left_test = X_left_test[,2:(M+1)]
X_left_test = data.matrix(X_left_test)

# Number of confounding variables 
P = ncol(X_confounds_train)

# Create list to store results in.
result_time = list()
result_nm = list()

# Start time of analysis.
time.start = Sys.time()

# Use 100-fold CV to identify regularisation parameter of Ridge regression and use lambda with lowest predictive MSE for analysis.
lambdas = 10^seq(2, -3, by = -.1)
cv_ridge = cv.glmnet(x = cbind(X_confounds_train, X_cortex_left_train), y =  y_train, alpha = 0, lambda = lambdas)
optimal_lambda = cv_ridge$lambda.min
ridge_reg = glmnet(x = cbind(X_confounds_train, X_cortex_left_train), y =  y_train, lambda = optimal_lambda, alpha = 0, family = 'gaussian')
beta_ridge = as.numeric(unlist(coefficients(ridge_reg)))[-c(1:(P+1))]
beta0_ridge = as.numeric(unlist(coefficients(ridge_reg)))[1]
beta_confound_ridge = as.numeric(unlist(coefficients(ridge_reg)))[2:(P+1)]

# Use 100-fold CV to identify regularisation parameter of LASSO regression and use lambda with lowest predictive MSE for analysis.
lambdas = 10^seq(2, -3, by = -.1)
cv_lasso = cv.glmnet(x = cbind(X_confounds_train, X_cortex_left_train), y =  y_train, alpha = 1, lambda = lambdas)
optimal_lambda_lasso = cv_lasso$lambda.min
lasso_model = glmnet(x = cbind(X_confounds_train, X_cortex_left_train), y =  y_train, alpha = 1, lambda = optimal_lambda_lasso)
beta_lasso = as.numeric(unlist(coefficients(lasso_model)))[-c(1:(P+1))]
beta0_lasso = as.numeric(unlist(coefficients(lasso_model)))[1]
beta_confound_lasso = as.numeric(unlist(coefficients(lasso_model)))[2:(P+1)]

# End time of analysis.
time.end = Sys.time()
time = time.end - time.start

# Save resulting time in list.
result_time$time_nm = time

# Concatenate estimated beta values of Ridge and LASSO regression in dataframe.
beta_fit = data.frame(LASSO = beta_lasso,
                       Ridge = beta_ridge)

# Save results in list.
beta_fit_left = beta_fit
result_nm$beta_fit = beta_fit
result_nm$ridge_reg = ridge_reg
result_nm$lasso_model = lasso_model
result_nm$time = time
save(result_nm, file = sprintf("%sresult_nm_left.RData", path_results_NM))

# Plot scatterplot of estimated beta values between Ridge and LASSO regression.
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
ggsave(p_plot, filename=paste0(path_results_NM, "scatterplot_beta_Ridge_vs_beta_LASSO_left.png"), type="cairo-png")

# Evaluation of prediction results
# Evaluate R2 (train) of Ridge regression model.
print('R2 Train: Ridge')
pred_nm = predict(ridge_reg, newx = cbind(X_confounds_train, X_cortex_left_train))
(R2_nm = cor(pred_nm,y_train)^2)

# Evaluate R2 (train) of LASSO regression model.
print('R2 Train: LASSO')
pred_hs = predict(lasso_model, newx = cbind(X_confounds_train, X_cortex_left_train))
(R2_hs = cor(pred_hs,y_train)^2)

# Evaluate bias (train) of Ridge regression model.
print('Bias Prediction Train: Ridge')
(bias_y_pred_nm = mean(pred_nm - y_train))

# Evaluate bias (train) of LASSO regression model.
print('Bias Prediction Train: LASSO')
(bias_y_pred_hs = mean(pred_hs - y_train))

# Evaluate MSE (train) of Ridge regression model.
print('MSE Prediction Train: Ridge')
(mse_y_pred_nm = mean((pred_nm - y_train)^2))

# Evaluate MSE (train) of LASSO regression model.
print('MSE Prediction Train: LASSO')
(mse_y_pred_hs = mean((pred_hs - y_train)^2))

# Evaluate R2 (test) of Ridge regression model.
print('R2 Test: Ridge')
pred_nm_test = predict(ridge_reg, newx = cbind(X_confounds_test, X_left_test))
(R2_nm_test = cor(pred_nm_test,y_test)^2)

# Evaluate R2 (test) of LASSO regression model.
print('R2 Test: LASSO')
pred_hs_test = predict(lasso_model, newx = cbind(X_confounds_test, X_left_test))
(R2_hs_test = cor(pred_hs_test,y_test)^2)

# Evaluate bias (test) of Ridge regression model.
print('Bias Prediction Test: Ridge')
(bias_y_pred_nm_test = mean(pred_nm_test - y_test))

# Evaluate bias (test) of LASSO regression model.
print('Bias Prediction Test: LASSO')
(bias_y_pred_hs_test = mean(pred_hs_test - y_test))

# Evaluate MSE (test) of Ridge regression model.
print('MSE Prediction Test: Ridge')
(mse_y_pred_nm_test = mean((pred_nm_test - y_test)^2))

# Evaluate MSE (test) of LASSO regression model.
print('MSE Prediction Test: LASSO')
(mse_y_pred_hs_test = mean((pred_hs_test - y_test)^2))

# Create tables to store predictive results in.
row_names = c('Ridge', 'LASSO')
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

# Read in imaging data of right hemisphere (test data)
X_right_test = fread(paste0(path_data, "X_right_cortex_val.csv"))
X_right_test = X_right_test[,2:(M+1)]
X_right_test = data.matrix(X_right_test)

# Create lists to store results in.
result_time = list()
result_nm = list()

# Start time of analysis.
time.start = Sys.time()

# Use 100-fold CV to determine best regularisation parameter for Ridge regression analysis.
lambdas = 10^seq(2, -3, by = -.1)
cv_ridge = cv.glmnet(x = cbind(X_confounds_train, X_cortex_right_train), y =  y_train, alpha = 0, lambda = lambdas)
optimal_lambda = cv_ridge$lambda.min
ridge_reg = glmnet(x = cbind(X_confounds_train, X_cortex_right_train), y =  y_train, lambda = optimal_lambda, alpha = 0, family = 'gaussian')
beta_ridge = as.numeric(unlist(coefficients(ridge_reg)))[-c(1:(P+1))]
beta0_ridge = as.numeric(unlist(coefficients(ridge_reg)))[1]
beta_confound_ridge = as.numeric(unlist(coefficients(ridge_reg)))[2:(P+1)]

# Use 100-fold CV to determine best regularisation parameter for LASSO regression analysis.
lambdas = 10^seq(2, -3, by = -.1)
cv_lasso = cv.glmnet(x = cbind(X_confounds_train, X_cortex_right_train), y =  y_train, alpha = 1, lambda = lambdas)
optimal_lambda_lasso = cv_lasso$lambda.min
lasso_model = glmnet(x = cbind(X_confounds_train, X_cortex_right_train), y =  y_train, alpha = 1, lambda = optimal_lambda_lasso)
beta_lasso = as.numeric(unlist(coefficients(lasso_model)))[-c(1:(P+1))]
beta0_lasso = as.numeric(unlist(coefficients(lasso_model)))[1]
beta_confound_lasso = as.numeric(unlist(coefficients(lasso_model)))[2:(P+1)]

# End time of analysis.
time.end = Sys.time()
time = time.end - time.start

# Save resulting time in list.
result_time$time_nm = time

# Concatenate estimated beta values of Ridge and LASSO in dataframe.
beta_fit = data.frame(LASSO = beta_lasso,
                       Ridge = beta_ridge)

# Save results in list.
beta_fit_right = beta_fit
result_nm$beta_fit = beta_fit
result_nm$ridge_reg = ridge_reg
result_nm$lasso_model = lasso_model
result_nm$time = time
save(result_nm, file = sprintf("%sresult_nm_right.RData", path_results_NM))

# Plot scatterplot between estimated beta values of Ridge and LASSO regression. 
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
ggsave(p_plot, filename=paste0(path_results_NM, "scatterplot_beta_Ridge_vs_beta_LASSO_right.png"), type="cairo-png")

# Evaluation of prediction results
# Evaluate R2 (train) of Ridge regression.
print('R2 Train: Ridge')
pred_nm = predict(ridge_reg, newx = cbind(X_confounds_train, X_cortex_right_train))
(R2_nm = cor(pred_nm,y_train)^2)

# Evaluate R2 (train) of LASSO regression.
print('R2 Train: LASSO')
pred_hs = predict(lasso_model, newx = cbind(X_confounds_train, X_cortex_right_train))
(R2_hs = cor(pred_hs,y_train)^2)

# Evaluate bias (train) of Ridge regression.
print('Bias Prediction Train: Ridge')
(bias_y_pred_nm = mean(pred_nm - y_train))

# Evaluate R2 (train) of LASSO regression.
print('Bias Prediction Train: LASSO')
(bias_y_pred_hs = mean(pred_hs - y_train))

# Evaluate MSE (train) of Ridge regression.
print('MSE Prediction Train: Ridge')
(mse_y_pred_nm = mean((pred_nm - y_train)^2))

# Evaluate MSE (train) of LASSO regression.
print('MSE Prediction Train: LASSO')
(mse_y_pred_hs = mean((pred_hs - y_train)^2))

# Evaluate R2 (test) of Ridge regression.
print('R2 Test: Ridge')
pred_nm_test = predict(ridge_reg, newx = cbind(X_confounds_test, X_right_test))
(R2_nm_test = cor(pred_nm_test,y_test)^2)

# Evaluate R2 (test) of LASSO regression.
print('R2 Test: LASSO')
pred_hs_test = predict(lasso_model, newx = cbind(X_confounds_test, X_right_test))
(R2_hs_test = cor(pred_hs_test,y_test)^2)

# Evaluate bias (test) of Ridge regression.
print('Bias Prediction Test: Ridge')
(bias_y_pred_nm_test = mean(pred_nm_test - y_test))

# Evaluate bias (test) of LASSO regression.
print('Bias Prediction Test: LASSO')
(bias_y_pred_hs_test = mean(pred_hs_test - y_test))

# Evaluate MSE (test) of Ridge regression.
print('MSE Prediction Test: Ridge')
(mse_y_pred_nm_test = mean((pred_nm_test - y_test)^2))

# Evaluate MSE (test) of LASSO regression.
print('MSE Prediction Test: LASSO')
(mse_y_pred_hs_test = mean((pred_hs_test - y_test)^2))

# Create tables to store predictive results in.
row_names = c('Ridge', 'LASSO')
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

# Create image map of Ridge regression estimated beta values.
xii = read_xifti(path_xii_image)
xii$data$cortex_left = matrix(beta_fit_left$Ridge, length(beta_fit_left$Ridge), 1)
xii$data$cortex_right = matrix(beta_fit_right$Ridge, length(beta_fit_right$Ridge), 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "beta_Ridge.dtseries.nii"), 
  file.path(path_results_NM, "beta_Ridge_left.surf.gii"), file.path(path_results_NM, "beta_Ridge_cortex_right.surf.gii")
)

# Create image map of LASSO regression estimated beta values.
xii = read_xifti(path_xii_image)
xii$data$cortex_left = matrix(beta_fit_left$LASSO, length(beta_fit_left$LASSO), 1)
xii$data$cortex_right = matrix(beta_fit_right$LASSO, length(beta_fit_right$LASSO), 1)

write_xifti(
  xii, 
  file.path(path_results_NM, "beta_LASSO.dtseries.nii"), 
  file.path(path_results_NM, "beta_LASSO_left.surf.gii"), file.path(path_results_NM, "beta_LASSO_cortex_right.surf.gii")
)
