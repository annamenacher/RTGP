# Script: Evaluate RTGP-VI model results for various kernel hyperparameter settings.

#############
### Paths ###
#############

# Path to folder where Human connectome software is located.
path_wb = ""

# Path to folder where data is stored in.
path_data = ""

# Path to folder where results with model estimates are stored in.
path_results = ""

# Path to folder where to save concatenated results in.
path_out = ""

#################
### Libraries ###
#################

options(bitmapType='cairo-png')
library(ciftiTools)
ciftiTools.setOption('wb_path', wb_path)
library(data.table)
library(xtable)
library(ggplot2)

##################
### Evaluation ###
##################

# Create vectors to store results (predictive MSE, ELBO, number of significant activations).
y_MSE_left = numeric(5)
y_MSE_right = numeric(5)
ELBO_left = numeric(5)
ELBO_right = numeric(5)
alpha_activation_left = matrix(NA, 5, 2)
alpha_activation_right = matrix(NA, 5, 2)

# Read in data matrix which contains confounding variables (training data)
X_confounds_train = fread(paste0(path_data, 'X_confounds_train.csv'))
X_confounds_train = data.matrix(X_confounds_train)[,2:ncol(X_confounds_train)]
# Read in imaging data of left hemisphere (training data)
X_cortex_left_train = fread(paste0(path_data, 'X_left_cortex_train.csv'))
X_cortex_left_train = data.matrix(X_cortex_left_train)[,2:ncol(X_cortex_left_train)]
# Number of vertices (left hemisphere)
M_left = ncol(X_cortex_left_train)
# Read in output vector of intelligence scores (training data)
y_train = as.numeric(unlist(fread(paste0(path_data, 'y_train.csv'))[,2]))

# Read in data matrix which contains confounding variables (validation data)
X_confounds_test = fread(paste0(path_data, "X_confounds_val.csv"))
X_confounds_test = X_confounds_test[,2:33]
X_confounds_test = data.matrix(X_confounds_test)
# Read in output vector of intelligence scores (validation data)
y_test = read.csv(paste0(path_data, "y_val.csv"))[,2]
# Read in imaging data of left hemisphere (validation data)
X_left_test = fread(paste0(path_data, "X_left_cortex_val.csv"))
X_left_test = X_left_test[,2:(M_left+1)]
X_left_test = data.matrix(X_left_test)

# Read in imaging data of right hemisphere (training data)
X_cortex_right_train = fread(paste0(path_data, 'X_right_cortex_train.csv'))
X_cortex_right_train = data.matrix(X_cortex_right_train)[,2:ncol(X_cortex_right_train)]
# Number of vertices (right hemisphere)
M_right = ncol(X_cortex_right_train)

# Read in imaging data of right hemisphere (validation data)
X_right_test = fread(paste0(path_data, "X_right_cortex_val.csv"))
X_right_test = X_right_test[,2:(M_right+1)]
X_right_test = data.matrix(X_right_test)

# Evaluate each of the 5 kernel hyperparameter settings.
for(sim in 1:5){

# Number of basis functions
L = 800
# Kernel hyperparameter: bandwidth parameter
phi_list = c(0.01, 0.05, 0.1, 0.2, 0.3)
phi = phi_list[sim]
# Kernel hyperparameter: exponent parameter
nu = 1.2

# Set current path to results.
path = paste0(path_results, toString(sim), "/")

### Left hemisphere ###

# Load results (left hemisphere). 
load(paste0(path, "result_rtgp_vi_left_cortex.RData"))
result_left = result

# Determine highest ELBO value when convergence was reached.
ELBO = result_left$ELBO
iter = table(ELBO == 0)[1]
log_likelihood = result_left$log_likelihood[2:iter] 
ELBO = ELBO[1:iter]
iter_left = which(ELBO == max(ELBO)) 

# Evaluate predictive results (left hemisphere)
beta0 = result_left$beta0[iter_left] 
beta = result_left$beta[iter_left,] 
beta_thresholded = result_left$beta_thresholded[iter_left,] 
beta_confound = result_left$beta_confound[iter_left,] 
alpha_thresholded = result_left$alpha_thresholded[iter_left,] 
alpha = result_left$alpha[iter_left,] 
beta0_chain = result_left$beta0[1:iter]
sigma_alpha_chain = result_left$sigma_alpha[1:iter]
sigma_beta_chain = result_left$sigma_beta[1:iter]
sigma_epsilon_chain = result_left$sigma_epsilon[1:iter]
threshold_chain = result_left$threshold[1:iter]

y_pred_train = as.vector(X_confounds_train%*% matrix(beta_confound, 32, 1))  + as.vector(X_cortex_left_train%*%matrix(beta_thresholded,M_left,1)) + beta0
y_pred_test = as.vector(X_confounds_test%*% matrix(beta_confound, 32, 1))  + as.vector(X_left_test%*%matrix(beta_thresholded,M_left,1)) + beta0

# Evaluate R2 (train) of RTGP model (left hemisphere)
result_R2_train = cor(y_pred_train, y_train)^2
# Evaluate R2 (test) of RTGP model (left hemisphere)
result_R2_test = cor(y_pred_test, y_test)^2

# Evaluate MSE (train) of RTGP model (left hemisphere)
result_MSE_train = mean((y_pred_train - y_train)^2)
# Evaluate MSE (test) of RTGP model (left hemisphere)
result_MSE_test = mean((y_pred_test - y_test)^2)

# Save results in list for predictive MSE, ELBO, and number of significant activations.
y_MSE_left[sim] = result_MSE_test
ELBO_left[sim] = ELBO[iter_left]
alpha_activation_left[sim,] = table(alpha_thresholded)

result_eval = matrix(NA, 1, 4)
result_eval[1,1] = result_R2_train
result_eval[1,2] = result_MSE_train
result_eval[1,3] = result_R2_test
result_eval[1,4] = result_MSE_test

# Create table for predictive results of left hemisphere.
col_names = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')
row_names = c('RTGP-VI')
colnames(result_eval) = col_names
rownames(result_eval) = row_names

xtable(result_eval, digits = 4)
print(xtable(result_eval, digits = 4), file=paste0(path, "table_result_left_cortex.txt"))
write.csv(result_eval, paste0(path, 'result_matrix_left.csv'))

# Save trace plot of ELBO values across iterations.
df = data.frame('ELBO' = ELBO[2:iter],
                'Iteration' = 2:length(ELBO))
p_plot = ggplot(df, aes(x = Iteration, y = ELBO)) + 
  geom_line(aes()) +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 1, linetype = "solid"),
        plot.background = element_rect(fill = "white")) +
  ggtitle('ELBO')
ggsave(p_plot, filename=paste0(path, "RTGP_VI_ELBO_left.png"), type="cairo-png")

### Right hemisphere ###

# Load results (right hemisphere). 
load(paste0(path, "result_rtgp_vi_right_cortex.RData"))
result_right = result

# Determine highest ELBO value when convergence was reached.
ELBO = result_right$ELBO
iter = table(ELBO == 0)[1]
log_likelihood = result_right$log_likelihood[2:iter] 
ELBO = ELBO[1:iter]
iter_right = which(ELBO == max(ELBO)) 

# Evaluate predictive results (right hemisphere)
beta0 = result_right$beta0[iter_right] 
beta = result_right$beta[iter_right,] 
beta_thresholded = result_right$beta_thresholded[iter_right,] 
beta_confound = result_right$beta_confound[iter_right,] 
alpha_thresholded = result_right$alpha_thresholded[iter_right,] 
alpha = result_right$alpha[iter_right,] 
beta0_chain = result_right$beta0[1:iter]
sigma_alpha_chain = result_right$sigma_alpha[1:iter]
sigma_beta_chain = result_right$sigma_beta[1:iter]
sigma_epsilon_chain = result_right$sigma_epsilon[1:iter]
threshold_chain = result_right$threshold[1:iter]

y_pred_train = as.vector(X_confounds_train%*% matrix(beta_confound, 32, 1))  + as.vector(X_cortex_left_train%*%matrix(beta_thresholded,M_left,1)) + beta0
y_pred_test = as.vector(X_confounds_test%*% matrix(beta_confound, 32, 1))  + as.vector(X_left_test%*%matrix(beta_thresholded,M_left,1)) + beta0

# Evaluate R2 (train) of RTGP model (right hemisphere)
result_R2_train = cor(y_pred_train, y_train)^2
# Evaluate R2 (test) of RTGP model (right hemisphere)
result_R2_test = cor(y_pred_test, y_test)^2

# Evaluate MSE (train) of RTGP model (right hemisphere)
result_MSE_train = mean((y_pred_train - y_train)^2)
# Evaluate MSE (test) of RTGP model (right hemisphere)
result_MSE_test = mean((y_pred_test - y_test)^2)

# Save results in list for predictive MSE, ELBO, and number of significant activations.
y_MSE_right[sim] = result_MSE_test
ELBO_right[sim] = ELBO[iter_left]
alpha_activation_right[sim,] = table(alpha_thresholded)

result_eval = matrix(NA, 1, 4)
result_eval[1,1] = result_R2_train
result_eval[1,2] = result_MSE_train
result_eval[1,3] = result_R2_test
result_eval[1,4] = result_MSE_test

# Create table for predictive results of right hemisphere.
col_names = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')
row_names = c('RTGP-VI')
colnames(result_eval) = col_names
rownames(result_eval) = row_names

xtable(result_eval, digits = 4)
print(xtable(result_eval, digits = 4), file=paste0(path, "table_result_right_cortex.txt"))
write.csv(result_eval, paste0(path, 'result_matrix_right.csv'))

# Save trace plot of ELBO values across iterations.
df = data.frame('ELBO' = ELBO[2:iter],
                'Iteration' = 2:length(ELBO))
p_plot = ggplot(df, aes(x = Iteration, y = ELBO)) + 
  geom_line(aes()) +
  theme_classic() + 
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 1, linetype = "solid"),
        plot.background = element_rect(fill = "white")) +
  ggtitle('ELBO')
ggsave(p_plot, filename=paste0(path, "RTGP_VI_ELBO_left.png"), type="cairo-png")

}

# Save ELBO, MSE, and number of significant results.
write.csv(ELBO_left, paste0(path_out, "ELBO_left.csv"))
write.csv(ELBO_right, paste0(path_out, "ELBO_right.csv"))
write.csv(y_MSE_left, paste0(path_out, "MSE_left.csv"))
write.csv(y_MSE_right, paste0(path_out, "MSE_right.csv"))
write.csv(alpha_activation_left, paste0(path_out, "activation_left.csv"))
write.csv(alpha_activation_right, paste0(path_out, "activation_right.csv"))

# Create plots to plot ELBO values across kernel hyperparameter settings (left hemisphere).
data <- data.frame(
  Eigenvalues = 1:length(ELBO_left),
  Percentage = ELBO_left
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('ELBO') +
  theme_minimal() +
  geom_point()

print(plot)

# Create plots to plot MSE values across kernel hyperparameter settings (left hemisphere).
data <- data.frame(
  Eigenvalues = 1:length(y_MSE_left),
  Percentage = y_MSE_left
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('MSE (Prediction)') +
  theme_minimal() +
  geom_point()

print(plot)

# Create plots to plot ELBO values across kernel hyperparameter settings (right hemisphere).
data <- data.frame(
  Eigenvalues = 1:length(ELBO_right),
  Percentage = ELBO_right
)

# Create a basic line plot
plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('ELBO') +
  theme_minimal() +
  geom_point()

# Display the plot
print(plot)

# Create plots to plot MSE values across kernel hyperparameter settings (right hemisphere).
data <- data.frame(
  Eigenvalues = 1:length(y_MSE_right),
  Percentage = y_MSE_right
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('MSE (Prediction)') +
  theme_minimal() +
  geom_point()

print(plot)
