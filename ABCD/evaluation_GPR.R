# Script: Evaluate GPR model results for various kernel hyperparameter settings.

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

# Create vectors and matrices to store results in of GPR + Normal and GPR + Horseshoe.
NM_ETI_activation_left = matrix(NA, 60, 2)
NM_ETI_activation_right = matrix(NA, 60, 2)
NM_HDI_activation_left = matrix(NA, 60, 2)
NM_HDI_activation_right = matrix(NA, 60, 2)

HS_ETI_activation_left = matrix(NA, 60, 2)
HS_ETI_activation_right = matrix(NA, 60, 2)
HS_HDI_activation_left = matrix(NA, 60, 2)
HS_HDI_activation_right = matrix(NA, 60, 2)

NM_MSE_test_left = numeric(60)
NM_MSE_test_right = numeric(60)
HS_MSE_test_left = numeric(60)
HS_MSE_test_right = numeric(60)

for(sim in 1:60){
  
  # Number of basis functions
  L = 800
  
  # Kernel hyperparameters (bandwidth parameter)
  phi_list = rep(c(0.01, 0.05, seq(0.1, 1, length.out = 10)), 5)
  phi = phi_list[sim]
  # Kernel hyperparameters (exponent parameter)
  nu_list = c(rep(1.2,12), rep(1.4,12), rep(1.6,12), rep(1.8,12), rep(2,12))
  nu = nu_list[sim]
  
  # Set current path to results
  path = paste0(path_out, "L", toString(L), "_phi", toString(phi), "_nu", toString(nu), "/with_confounds/")
  
  # Read in results of left and right hemisphere.
  load(paste0(path, "result_left_cortex.RData"))
  load(paste0(path, "result_right_cortex.RData"))

  NM_ETI_activation_left[sim,] = table(result_left_cortex$significance_CI_ETI_NM)
  NM_HDI_activation_left[sim,] = table(result_left_cortex$significance_CI_HPDI_NM)
  HS_ETI_activation_left[sim,] = table(result_left_cortex$significance_CI_ETI_HS)
  HS_HDI_activation_left[sim,] = table(result_left_cortex$significance_CI_HPDI_HS)
  
  NM_ETI_activation_right[sim,] = table(result_right_cortex$significance_CI_ETI_NM)
  NM_HDI_activation_right[sim,] = table(result_right_cortex$significance_CI_HPDI_NM)
  HS_ETI_activation_right[sim,] = table(result_right_cortex$significance_CI_ETI_HS)
  HS_HDI_activation_right[sim,] = table(result_right_cortex$significance_CI_HPDI_HS)
  
  NM_MSE_test_left[sim] = result_left_cortex$result_R2_MSE$mse_y_pred_test_NM
  HS_MSE_test_left[sim] = result_left_cortex$result_R2_MSE$mse_y_pred_test_HS
  NM_MSE_test_right[sim] = result_right_cortex$result_R2_MSE$mse_y_pred_test_NM
  HS_MSE_test_right[sim] = result_right_cortex$result_R2_MSE$mse_y_pred_test_HS
  
}

NM_ETI_activation_right[NM_ETI_activation_right == 29716] = 0
NM_HDI_activation_right[NM_HDI_activation_right == 29716] = 0
HS_ETI_activation_right[HS_ETI_activation_right == 29716] = 0
HS_HDI_activation_right[HS_HDI_activation_right == 29716] = 0

NM_ETI_activation_left[NM_ETI_activation_left == 29696] = 0
NM_HDI_activation_left[NM_HDI_activation_left == 29696] = 0
HS_ETI_activation_left[HS_ETI_activation_left == 29696] = 0
HS_HDI_activation_left[HS_HDI_activation_left == 29696] = 0


# Plot the number of significant results (evaluated via ETI) in GPR + Normal model across kernel hyperparameter settings.
data <- data.frame(
  Group = rep(c("Left", "Right"), each = 60),
  Variable = rep(1:60, 2),
  Value = c(NM_ETI_activation_left[,2], NM_ETI_activation_right[,2])
)

plot <- ggplot(data, aes(x = Variable, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Hyperparameter Setting",
       y = "Activation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)


# Plot the number of significant results (evaluated via HPDI) in GPR + Normal model across kernel hyperparameter settings.
set.seed(123)
data <- data.frame(
  Group = rep(c("Left", "Right"), each = 60),
  Variable = rep(1:60, 2),
  Value = c(NM_HDI_activation_left[,2], NM_HDI_activation_right[,2])
)

plot <- ggplot(data, aes(x = Variable, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Hyperparameter Setting",
       y = "Activation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)

# Plot the number of significant results (evaluated via ETI) in GPR + Horseshoe model across kernel hyperparameter settings.
data <- data.frame(
  Group = rep(c("Left", "Right"), each = 60),
  Variable = rep(1:60, 2),
  Value = c(HS_ETI_activation_left[,2], HS_ETI_activation_right[,2])
)

plot <- ggplot(data, aes(x = Variable, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Hyperparameter Setting",
       y = "Activation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)


# Plot the number of significant results (evaluated via HPDI) in GPR + Horseshoe model across kernel hyperparameter settings.
data <- data.frame(
  Group = rep(c("Left", "Right"), each = 60),
  Variable = rep(1:60, 2),
  Value = c(HS_HDI_activation_left[,2], HS_HDI_activation_right[,2])
)

# Create a bar plot
plot <- ggplot(data, aes(x = Variable, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Hyperparameter Setting",
       y = "Activation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)


# Plot predictive MSE results of GPR + Normal model across kernel hyperparameter settings (left hemisphere).
data <- data.frame(
  Eigenvalues = 1:length(NM_MSE_test_left),
  Percentage = NM_MSE_test_left
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('NM Left MSE (Prediction)') +
  theme_minimal() +
  geom_point()

print(plot)


# Plot predictive MSE results of GPR + Normal model across kernel hyperparameter settings (right hemisphere).
data <- data.frame(
  Eigenvalues = 1:length(NM_MSE_test_right),
  Percentage = NM_MSE_test_right
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('NM Right MSE (Prediction)') +
  theme_minimal() +
  geom_point()

print(plot)


# Plot predictive MSE results of GPR + Horseshoe model across kernel hyperparameter settings (left hemisphere).
data <- data.frame(
  Eigenvalues = 1:length(HS_MSE_test_left),
  Percentage = HS_MSE_test_left
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('HS Left MSE (Prediction)') +
  theme_minimal() +
  geom_point()

print(plot)


# Plot predictive MSE results of GPR + Horseshoe model across kernel hyperparameter settings (right hemisphere).
data <- data.frame(
  Eigenvalues = 1:length(HS_MSE_test_right),
  Percentage = HS_MSE_test_right
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('HS Right MSE (Prediction)') +
  theme_minimal() +
  geom_point()

print(plot)



###########################################################################
### Consistency check: GPR model with and without confounding variables ###
###########################################################################

# Create vectors and matrices to store results in.
y_MSE_left = numeric(60)
y_MSE_right = numeric(60)
R2_left = numeric(60)
R2_right = numeric(60)
alpha_activation_left = matrix(NA, 60, 2)
alpha_activation_right = matrix(NA, 60, 2)

for(sim in 1:60){
  
  # Kernel hyperparameters: bandwidth parameter.
  phi_list = rep(c(0.01, 0.05, seq(0.1, 1, length.out = 10)), 5)
  phi = phi_list[sim]
  # Kernel hyperparameters: exponent parameter.
  nu_list = c(rep(1.2,12), rep(1.4,12), rep(1.6,12), rep(1.8,12), rep(2,12))
  nu = nu_list[sim]
  
  # Set current path to results.
  path = paste0(path_out, "L", toString(L), "_phi", toString(phi), "_nu", toString(nu), "/no_confounds/")
  
  # Read in results for left and right hemisphere.
  result_left = read.csv(paste0(path, 'result_matrix.csv'), header = T)[2,2:5]
  y_MSE_left[sim] = result_left[1,4]
  R2_left[sim] = result_left[1,3]
  result_right = read.csv(paste0(path, 'result_matrix_right.csv'), header = T)[2,2:5]
  y_MSE_right[sim] = result_right[1,4]
  R2_right[sim] = result_right[1,3]
  load(paste0(path, "result_left_cortex.RData"))
  alpha_activation_left[sim,] = table(result_left_cortex$significance_CI_HPDI_HS)
  load(paste0(path, "result_right_cortex.RData"))
  alpha_activation_right[sim,] = table(result_right_cortex$significance_CI_HPDI_HS)
  
}

# Plot comparison between number of significant results between left and right hemisphere across kernel hyperparameter settings.
data <- data.frame(
  Group = rep(c("Left", "Right"), each = 60),
  Variable = rep(1:60, 2),
  Value = c(alpha_activation_left[,2], alpha_activation_right[,2])
)

plot <- ggplot(data, aes(x = Variable, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Hyperparameter Setting",
       y = "Activation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)

# Plot comparison between R2 predictive results between left and right hemisphere across kernel hyperparameter settings.
data <- data.frame(
  Group = rep(c("Left", "Right"), each = 60),
  Variable = rep(1:60, 2),
  Value = c(R2_left, R2_right)
)

plot <- ggplot(data, aes(x = Variable, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Hyperparameter Setting",
       y = "R2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)




# Plot comparison between MSE predictive results between left and right hemisphere across kernel hyperparameter settings.
data <- data.frame(
  Group = rep(c("Left", "Right"), each = 60),
  Variable = rep(1:60, 2),
  Value = c(y_MSE_left, y_MSE_right)
)

plot <- ggplot(data, aes(x = Variable, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Hyperparameter Setting",
       y = "MSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)










