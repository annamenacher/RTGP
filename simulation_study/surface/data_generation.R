# Script: Generate datasets for cortical surface-based simulation study.

#############
### Paths ###
#############

# Path with folder here the Human Connectome Workbench software is located.
path_wb_bench = ""
# Path with folder where true beta is located and where the datasets should be saved.
path_core = ""

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
ciftiTools.setOption('wb_path', path_wb_bench)

#################
### Constants ###
#################

# Number of datasets
n_dataset = 10

# Number of vertices
p = 1799 
# Number of subjects in the training dataset
n_train = 500 # n_train = 500 # n_train = 1000 # n_train = 2000
# Number of subjects in the test dataset
n_test = 1000
# Total number of subjects 
n = n_train + n_test

for(dataset in 1:n_dataset){

  # True intercept value
  beta0_true = 2
  # Kernel hyperparameters
  phi = 0.077
  nu = 2
  L = 100
  # Residual noise value
  noise = 0.2
  true_sigma2 = noise^2
  
  # True beta values
  beta_true = read.csv(paste0(path_beta, "beta_left_2k.csv"))[,2]
  
  # Generate dataset specific folders.
  path_data_folder = sprintf("%sN%s_L%s_phi%s_nu%s/", path_core, n_train, L, phi, nu)
  path_data = path_results = sprintf("%sN%s_L%s_phi%s_nu%s/dataset%s/", path_core, n_train, L, phi, nu, dataset)
  
  dir.create(file.path(path_data_folder), showWarnings = FALSE)
  dir.create(file.path(path_data), showWarnings = FALSE)
  
  # Read in eigenvalues and eigenfunctions
  Psi = data.matrix(fread(paste0(path_core,"eigendecomposition/eigendecomposition_L100_phi0.077_nu2/eigenvectors_left_cortex.csv")))
  Psi = Psi[,2:ncol(Psi)]
  lambda = read.csv(paste0(path_core, "eigendecomposition/eigendecomposition_L100_phi0.077_nu2/eigenvalues_left_cortex.csv"))[,2]
  # Square rooted eigenvalues
  sqrt_lambda = sqrt(lambda)
  # Bases defined by eigenvalues and eigenfunctions
  Bases = t(Psi)*sqrt_lambda
  # Generate random normal values to simulate input images from a GP
  theta_X = matrix(rnorm(n*L,sd = 1/sqrt(p)),nrow=L,ncol=n)
  img_total = crossprod(Bases,theta_X)
  # Calculate scalar output values
  mean_Y = beta0_true + crossprod(img_total,beta_true)
  y_total = mean_Y + rnorm(n,sd=sqrt(true_sigma2))
  # Transform input images with bases
  Z_total = cbind(1,t(Bases%*%img_total))
  
  # Split data in training and test dataset
  train_idx_SOI = sample(1:n,round(n_train))
  test_idx_SOI = setdiff(1:n,train_idx_SOI)
  
  y = y_total[train_idx_SOI]
  Z =  Z_total[train_idx_SOI,]
  y_test = y_total[test_idx_SOI,]
  img_test = img_total[,test_idx_SOI]
  Z_test =  Z_total[test_idx_SOI,]
  img = img_total[,train_idx_SOI]
  
  # Save example input image of subject 1
  xii = read_xifti('/well/nichols/users/fxo071/IOI-GP/sim_study/11.21.23/data/beta_2k.dtseries.nii')
  xii$data$cortex_left = matrix(img[,1], length(img[,1]), 1)
  xii$data$cortex_right = matrix(0, 1803, 1)
  
  write_xifti(
    xii, 
    file.path(path_results, "X1.dtseries.nii"), 
    file.path(path_results, "X1_left.surf.gii"), file.path(path_results, "X1_cortex_right.surf.gii")
  )
  
  # Save example input image of subject 2
  xii = read_xifti('/well/nichols/users/fxo071/IOI-GP/sim_study/11.21.23/data/beta_2k.dtseries.nii')
  xii$data$cortex_left = matrix(img[,2], length(img[,2]), 1)
  xii$data$cortex_right = matrix(0, 1803, 1)
  
  write_xifti(
    xii, 
    file.path(path_results, "X2.dtseries.nii"), 
    file.path(path_results, "X2_left.surf.gii"), file.path(path_results, "X2_cortex_right.surf.gii")
  )
  
  # Save example input image of subject 3
  xii = read_xifti('/well/nichols/users/fxo071/IOI-GP/sim_study/11.21.23/data/beta_2k.dtseries.nii')
  xii$data$cortex_left = matrix(img[,3], length(img[,3]), 1)
  xii$data$cortex_right = matrix(0, 1803, 1)
  
  write_xifti(
    xii, 
    file.path(path_results, "X3.dtseries.nii"), 
    file.path(path_results, "X3_left.surf.gii"), file.path(path_results, "X3_cortex_right.surf.gii")
  )
  
  # Determine true binary significance
  alpha_true = ifelse(abs(beta_true) > 0, 1, 0) 
  
  # Save true binary significance map
  xii = read_xifti('/well/nichols/users/fxo071/IOI-GP/sim_study/11.21.23/data/beta_2k.dtseries.nii')
  xii$data$cortex_left = matrix(alpha_true, length(img[,1]), 1)
  xii$data$cortex_right = matrix(0, 1803, 1)
  
  write_xifti(
    xii, 
    file.path(path_results, "significance_true.dtseries.nii"), 
    file.path(path_results, "significance_true_left.surf.gii"), file.path(path_results, "significance_true_cortex_right.surf.gii")
  )
  
  # Save generated dataset in list
  data = list()
  data$X = t(img)
  data$y = y
  data$Psi = t(Psi)
  data$lambda = lambda
  data$beta_true = beta_true
  data$alpha_true = alpha_true
  data$true_sigma2 = true_sigma2
  data$Z = Z
  data$y_test = y_test
  data$img_test = img_test
  data$Z_test = Z_test
  data$p = p
  data$n_train = n_train
  data$n_test = n_test
  data$n = n
  data$beta0_true = beta0_true
  save(data, file = paste0(path_data, 'data.RData'))

}