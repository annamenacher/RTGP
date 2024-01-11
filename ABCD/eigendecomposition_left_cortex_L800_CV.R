# Script: Run eigen decompositon of kernel matrix for various kernel hyperparameter settings on the left hemisphere.

# Determine which kernel hyperparameter setting is chosen for kernel decomposition.
sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#############
### Paths ###
#############

# Path to folder which contains data.
path = ''
# Path to folder where to save eigendecomposition results in.
path_out = ''
dir.create(file.path(path_out), showWarnings = FALSE)
path_out = paste0(path_out, 'eigendecomposition_L', toString(L), '_phi', toString(phi), '_nu', toString(nu), '/')
dir.create(file.path(path_out), showWarnings = FALSE)

#################
### Libraries ###
#################

library(data.table)
library(RSpectra)

#################
### Constants ###
#################

# Create list of hyperparameter configurations to test (bandwidth parameter)
phi_list = rep(c(0.01, 0.05, 0.1, 0.2, 0.3), 2)
phi = phi_list[sim]
# Create list of hyperparameter configurations to test (exponent parameter)
nu_list = c(rep(1.2,5), rep(1.4,5))
nu = nu_list[sim]
# Number of basis coefficients
L = 800

##########################
### Eigendecomposition ###
##########################

# Read in distance matrix generated via package 'brainsmash'.
distance_matrix = data.matrix(fread(paste0(path, 'LeftDenseGeodesicDistmat.txt')))

print(dim(distance_matrix))

# Read in mask of missing vertices in data from template. 
mask = as.numeric(unlist(fread(paste0(path, 'mask/idx_left.csv'))[,2]))

# Mask out missing vertices from distance matrix.
distance_matrix = distance_matrix[mask, mask]

# Two parameter exponential radial basis function
covariance_kernel = function(x, s, phi, nu){
  y = s * exp(-phi*x^nu)
  return(y)
}

# Kernel matrix 
distance_matrix = covariance_kernel(distance_matrix, s = 1, phi = phi, nu = nu)

# Eigendecomposition
eigendecomp = eigs_sym(A = distance_matrix, k = L, which = "LM") 

rm(distance_matrix)

# Save results of eigendecomposition (eigenvalues, eigenfunctions)
write.csv(eigendecomp$values, paste0(path_out, 'eigenvalues_left_cortex.csv'))
write.csv(eigendecomp$vectors, paste0(path_out, 'eigenvectors_left_cortex.csv'))
save(eigendecomp, file = paste0(path_out, 'eigendecomposition_left_cortex.RData'))

