# Script: Evaluate if RTGP model reliably converges to the same ELBO for different initialisation.

#############
### Paths ###
#############

# Path to Human conncectome workbench software.
path_wb = ""
# Path to RTGP results with different initialisation.
path = ""

#################
### Libraries ###
#################

options(bitmapType='cairo-png')
library(ciftiTools)
ciftiTools.setOption('wb_path', path_wb)
library(data.table)
library(xtable)
library(ggplot2)

#################
### Constants ###
#################

# Sample size of training dataset
n_train = 2000
# Variance of prior on intercept
sigma_beta0 = 10
# True residual noise 
sigma_epsilon = 0.2

##################
### Evaluation ###
##################

# Calculate true log-likelihood value
mu = beta0_true + X%*%beta_true
ll_true = sum(dnorm(y, mu, sigma_epsilon, log=TRUE))

# Create vectors to store results.
ELBO = numeric(5)
sigma_epsilon = numeric(5)
sigma_alpha = numeric(5)
sigma_beta = numeric(5)
beta0 = numeric(5)
threshold  = numeric(5)
ll = numeric(5)
activation = numeric(5)

for(sim in 1:5){
  
  # Set path to folder with results.
  path_out = paste0(path, 'sim', toString(sim), '/RTGP_VI/')
  # Load results.
  load(paste0(path_out, 'result_rtgp_vi_left_cortex.RData'))
  # Determine which iteration has the maximum ELBO.
  iter = which(result$ELBO == max(result$ELBO))
  # Store results in vectors.
  ELBO[sim] = result$ELBO[iter]
  sigma_epsilon[sim] = result$sigma_epsilon[iter]
  sigma_alpha[sim] = result$sigma_alpha[iter]
  sigma_beta[sim] = result$sigma_beta[iter]
  beta0[sim] = result$beta0[iter]
  threshold[sim] = result$threshold[iter]
  ll[sim] = result$log_likelihood[iter]
  activation[sim] = table(result$alpha_thresholded[iter,])[2]
}

# Create table comparing posterior estimates of 5 runs with different initialisations.
matrix_table = matrix(NA, 8, 5)
matrix_table[1,] = ELBO
matrix_table[2,] = ll
matrix_table[3,] = threshold
matrix_table[4,] = beta0
matrix_table[5,] = sigma_beta
matrix_table[6,] = sigma_alpha
matrix_table[7,] = sigma_epsilon
matrix_table[8,] = activation

row_names = c('ELBO', 'log-likelihood', 'threshold', 'beta0', 'sigma_beta', 'sigma_alpha', 'sigma_epsilon', 'activation')
column_names = c('Run 1', 'Run 2', 'Run 3', 'Run 4', 'Run 5')
rownames(matrix_table) = row_names
colnames(matrix_table) = column_names

# Print and store table.
xtable(matrix_table, digits = 4)
print(xtable(matrix_table, digits = 4), file=paste0(path, "matrix_table_model_coefficients.txt"))
write.csv(matrix_table, paste0(path, 'matrix_table_model_coefficients.csv'))


# Create a matrix to store ELBO values across iterations (maximum number of iterations is 1000).
ELBO_df = matrix(NA, 5, 1000)

for(sim in 1:5){
  # Set path.
  path_out = paste0(path, 'sim', toString(sim), '/RTGP_VI/')
  # Load results.
  load(paste0(path_out, 'result_rtgp_vi_left_cortex.RData'))
  # Store ELBO values across iterations for initialisation run.
  ELBO_df[sim,] = result$ELBO
}

# Plot ELBO values across iterations for RTGP model with 5 different initialisations.
data = data.frame(
  x = 2:100,                
  y = t(ELBO_df[,2:100])  
)
data_long <- tidyr::gather(data, key = "group", value = "value", -x)

ggplot(data_long, aes(x = x, y = value, color = group)) +
  geom_line() +
  labs(title = "ELBO (with 5 different initialisations)",
       x = "Iterations",
       y = "ELBO",
       color = "Run") +
  theme_minimal()

# Create matrices to store inference, parameter estimate, and prediction results for 5 runs with different initialisations.
inference_df = matrix(NA, 5, 4)
parameters_df = matrix(NA, 5, 2)
prediction_df = matrix(NA, 5, 4)

for(sim in 1:5){
  # Set path.
  path_out = paste0(path, 'sim', toString(sim), '/RTGP_VI/')
  # Load results.
  load(paste0(path_out, 'result_matrix.RData'))
  # Store inference, parameter estimates, and prediction results.
  inference_df[sim,] = as.vector(as.numeric(unlist(result_matrix$matrix_inference[7,])))
  parameters_df[sim,] = as.vector(as.numeric(unlist(result_matrix$matrix_parameters[3,])))
  prediction_df[sim,] = as.vector(as.numeric(unlist(result_matrix$matrix_prediction[3,])))
  
}

# Create tables with inference, parameter estimates, and prediction results.
row_names = c('Run 1', 'Run 2', 'Run 3', 'Run 4', 'Run 5')
col_names_inference = c('TPR', 'TDR', 'FPR', 'FDR')
col_names_parameters = c('Bias', 'MSE')
col_names_prediction = c('R2 (train)', 'MSE (train)', 'R2 (test)', 'MSE (test)')

colnames(inference_df) = col_names_inference
rownames(inference_df) = row_names
colnames(parameters_df) = col_names_parameters
rownames(parameters_df) = row_names
colnames(prediction_df) = col_names_prediction
rownames(prediction_df) = row_names

xtable(inference_df, digits = 4)
print(xtable(inference_df, digits = 4), file=paste0(path, "inference.txt"))
write.csv(inference_df, paste0(path, 'matrix_inference.csv'))

xtable(parameters_df, digits = 4)
print(xtable(parameters_df, digits = 4), file=paste0(path, "parameters.txt"))
write.csv(parameters_df, paste0(path, 'matrix_parameters.csv'))

xtable(prediction_df, digits = 4)
print(xtable(inference_df, digits = 4), file=paste0(path, "prediction.txt"))
write.csv(prediction_df, paste0(path, 'matrix_prediction.csv'))

# Save tables.
write.csv(inference_df, paste0(path, 'inference.csv'))
write.csv(parameters_df, paste0(path, 'parameters.csv'))
write.csv(prediction_df, paste0(path, 'prediction.csv'))

