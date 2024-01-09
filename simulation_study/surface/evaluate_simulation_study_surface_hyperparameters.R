# Script: Evaluate hyperparameter settings of surface-based simulation study to determine kernel hyperparameters
# with highest yielding ELBO (for RTGP).

#############
### Paths ###
#############

# Path to Human conncectome workbench software.
path_wb = ""
# Path to save results in.
path = ""
# Path to RTGP results with different hyperparameter configurations.
path_sim = ''


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

# Number of training dataset subjects 
N = 2000 # N = 500 # N = 1000 # N = 2000
# Number of basis functions
L = 100 

##################
### Evaluation ###
##################

# Create vectors and matrices to store results of 60 different kernerl hyperparameter settings.
ELBO_list = numeric(60)
MSE_list = numeric(60)
alpha_activation = matrix(NA, 60, 2)

for(sim in 1:60){

# Kernel hyperparameter grid: bandwidth parameter
phi_list = rep(c(0.01, 0.05, seq(0.1, 1, length.out = 10)), 5)
phi = phi_list[sim]
# Kernel hyperparameter grid: exponent parameter
nu_list = c(rep(1.2,12), rep(1.4,12), rep(1.6,12), rep(1.8,12), rep(2,12))
nu = nu_list[sim]

# Path to specific hyperparameter setting.
path_sim = paste0(path_sim, 'N', toString(N), '_L', toString(L), '_phi', toString(phi), '_nu', toString(nu), '/dataset1/RTGP_VI/')

# Load results from RTGP.
load(paste0(path_sim, 'result_rtgp_vi_left_cortex.RData'))
load(paste0(path_sim, 'result_matrix.RData'))

# Determine iteration with maximum ELBO value of RTGP model.
iter = which(result$ELBO == max(result$ELBO))

# Extract ELBO value.
ELBO_list[sim] = result$ELBO[iter]
# Extract out-of-sample test MSE.
MSE_list[sim] = result_matrix$matrix_prediction[3,4]
# Extract number of significant vertices.
alpha_activation[sim,] = table(result$alpha_thresholded[iter,])

}

# Save results.
write.csv(ELBO_list, paste0(path, 'ELBO_list.csv'))
write.csv(MSE_list, paste0(path, 'MSE_list.csv'))
write.csv(alpha_activation, paste0(path, 'activation_list.csv'))


# Plot ELBO values across hyperparameter settings.
data <- data.frame(
  Eigenvalues = 1:length(ELBO_list),
  Percentage = ELBO_list
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('ELBO') +
  theme_minimal() +
  geom_point()

print(plot)

# Plot MSE values across hyperparameter settings.
data <- data.frame(
  Eigenvalues = 1:length(MSE_list),
  Percentage = MSE_list
)

plot <- ggplot(data, aes(Eigenvalues, Percentage)) +
  xlab('Hyperparameter Setting') +
  ylab('MSE (Prediction)') +
  theme_minimal() +
  geom_point()

print(plot)


# Plot number of significant values across hyperparameter settings.
set.seed(123)
data <- data.frame(
  Variable = rep(1:60, 1),
  Value = c(alpha_activation[,2])
)

plot <- ggplot(data, aes(x = Variable, y = Value)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "",
       x = "Hyperparameter Setting",
       y = "Activation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)

# Determine which hyperparameter setting to use for the rest of the analysis.
# Best setting according to maximum ELBO.
which(ELBO_list == max(ELBO_list))
phi_list[which(ELBO_list == max(ELBO_list))]
nu_list[which(ELBO_list == max(ELBO_list))]

# Best setting according to minimum MSE
which(MSE_list == min(MSE_list))
phi_list[which(MSE_list == min(MSE_list))]
nu_list[which(MSE_list == min(MSE_list))]

