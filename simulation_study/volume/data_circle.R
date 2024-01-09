# Libraries
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

# Functions

plot_2D_funcs = function(f,grids,names=NULL,nrow=NULL,ncol=NULL){
  
  # Plotting functions from Jian Kang to plot 2D images.
  
  fdat = data.frame(f,v1=grids[,1],v2=grids[,2])
  if(!is.null(names)){
    names(fdat) = c(names,"v1","v2")
  }
  fdat1 = melt(fdat,id.vars=c("v1","v2"))
  return(ggplot(fdat1,aes(x=v1,y=v2,fill=value))+facet_wrap(~variable,nrow=nrow,ncol=ncol)+coord_equal()+geom_tile()+scale_fill_gradientn(colours =GP.create.cols()))
}

# Number of datasets
n_dataset = 10

for(dataset in 1:n_dataset){

# Number of voxels
p = 2500
# Training sample size
n_train = 500 # n_train = 500 # n_train = 1000 # n_train = 2000
# Test sample size
n_test = 1000
# Total number of subjects: training + test dataset
n = n_train + n_test

# True intercept value
beta0_true = 2
# Poly-degree for Hermite polynomials
poly_degree = 10 # poly_degree = 10 # poly_degree = 20 # poly_degree = 30
# Kernel hyperparameters of modified squared exponential kernel
a = 0.01
b = 100
# True signal
signal = 1 # signal = 3 # signal = 1 
# True noise
noise = 0.2 # noise = 1 # noise = 0.2
true_sigma2 = noise^2
# Radius of effect circle
r = 0.07 # r = 0.07 # r = 0.13 # r = 0.27

# 2D grid points
v_list = GP.generate.grids(d = 2, num_grids = sqrt(p), grids_lim=c(-1,1))
# True beta vector
beta_true = as.numeric(v_list[,1]^2+v_list[,2]^2<r)
beta_true = beta_true * signal
# Number of active voxels 
effect = sum(beta_true)

# Path for saving the resulting datasets.
path_core = ""
path_data_folder = sprintf("%sN%s_poly%s_signal%s_noise%s_effect%s/", path_core, n_train, poly_degree, signal, noise, effect)
path_data = sprintf("%sN%s_poly%s_signal%s_noise%s_effect%s/dataset%s/", path_core, n_train, poly_degree, signal, noise, effect,dataset)

dir.create(file.path(path_data_folder), showWarnings = FALSE)
dir.create(file.path(path_data), showWarnings = FALSE)

# Perform eigendecomposition to acquire eigenfunctions Psi and eigenvalues lambda.
Psi = GP.eigen.funcs.fast(v_list, poly_degree = poly_degree, a = a, b = b)
lambda = GP.eigen.value(poly_degree = poly_degree, a = a, b = b, d = ncol(v_list))
# Square rooted eigenvalues.
sqrt_lambda = sqrt(lambda)
# Bases functions.
Bases = t(Psi)*sqrt_lambda
# Number of eigenvalues.
L = nrow(Bases)
# Random Normal values to generate input images with.
theta_X = matrix(rnorm(n*L,sd = 1/sqrt(p)),nrow=L,ncol=n)
# Input images drawn from a GP.
img_total = crossprod(Bases,theta_X)
# Output scalar.
mean_Y = beta0_true + crossprod(img_total,beta_true)
y_total = mean_Y + rnorm(n,sd=sqrt(true_sigma2))
# Transformed input images by multiplying input images with Bases and adding an intercept column.
Z_total = cbind(1,t(Bases%*%img_total))

# Split data in training and test data.
train_idx_SOI = sample(1:n,round(n_train))
test_idx_SOI = setdiff(1:n,train_idx_SOI)

y = y_total[train_idx_SOI]
Z =  Z_total[train_idx_SOI,]
y_test = y_total[test_idx_SOI,]
img_test = img_total[,test_idx_SOI]
Z_test =  Z_total[test_idx_SOI,]
img = img_total[,train_idx_SOI]

# Plot 4 examples of input images.
p_plot = plot_2D_funcs(img[,1:4],v_list)
ggsave(p_plot, filename=paste0(path_data, "X.png"), type="cairo-png")

# Plot true beta image.
p_plot = plot_2D_funcs(beta_true,v_list)
ggsave(p_plot, filename=paste0(path_data, "truth_beta.png"), type="cairo-png")

# Set and plot true significance in the areas of true effect determined by the true beta map.
alpha_true = beta_true

# Plot true significance map.
p_plot = plot_2D_funcs(beta_true / signal ,v_list)
ggsave(p_plot, filename=paste0(path_data, "truth_binary_significance.png"), type="cairo-png")

# Save list of generated quantities for each created dataset.
data = list()
data$X = t(img)
data$y = y
data$v_list = v_list
data$Psi = t(Psi)
data$lambda = lambda
data$beta_true = beta_true
data$alpha_true = alpha_true
data$true_sigma2 = true_sigma2
data$y_test = y_test
data$img_test = img_test
data$Z_test = Z_test
data$p = p
data$n_train = n_train
data$n_test = n_test
data$n = n
data$poly_degree = poly_degree
data$a = a 
data$b = b 
data$signal = signal 
data$noise = noise 
data$r = r
data$effect = effect 
data$beta0_true = beta0_true
save(data, file = paste0(path_data, 'data.RData'))
}
